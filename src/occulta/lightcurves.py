"""Light-curve level search: detrending, single-event box search, per-sector BLS.

Runs inside worker processes; everything here is pure numpy/astropy on one
QLP light curve file, and returns plain dicts for the orchestrator to store.
"""
import numpy as np
from scipy.ndimage import median_filter
from astropy.io import fits

WIDTHS_H = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.5, 9.0, 13.0]
DETREND_WINDOW_D = 2.0
GAP_D = 0.5
SNR_EVENT = 7.0
MAX_EVENTS = 5


def mad_std(x):
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


KEPLER_TO_BTJD = 2454833.0 - 2457000.0


def read_lc(path, fmt='qlp'):
    """(t [BTJD], flux, background, centroid x, centroid y, meta) of good cadences for any supported product."""
    with fits.open(path, memmap=False) as h:
        hdr = h[0].header
        if fmt == 'qlp':
            d = h[1].data
            # QLP's DET_FLUX/KSPSAP_FLUX are high-passed and erode transits longer than a few hours,
            # so use the systematics-corrected (or raw) aperture flux and detrend it ourselves.
            fcol = 'SYS_RM_FLUX' if 'SYS_RM_FLUX' in d.columns.names else 'SAP_FLUX'
            t, f, q = d['TIME'], d[fcol], d['QUALITY']
            bkg, x, y = d['SAP_BKG'], d['SAP_X'], d['SAP_Y']
            off, grp = 0.0, (hdr.get('CAMERA'), hdr.get('CCD'))
        elif fmt == 'k2sff':
            d = h['BESTAPER'].data
            t, f = d['T'], d['FCOR']
            q = np.where(np.asarray(d['MOVING']) > 0, 1, 0) if 'MOVING' in d.columns.names else np.zeros(len(t))
            bkg = x = y = np.full(len(t), np.nan)
            off, grp = KEPLER_TO_BTJD, (hdr.get('MODULE'), hdr.get('OUTPUT'))
        else:  # 'spoc' (TESS 2-min / 20-s) and 'kepler' share the SPOC layout
            d = h[1].data
            t, f = d['TIME'], d['PDCSAP_FLUX']
            q = d['QUALITY'] if 'QUALITY' in d.columns.names else d['SAP_QUALITY']
            bkg, x, y = d['SAP_BKG'], d['MOM_CENTR1'], d['MOM_CENTR2']
            if fmt == 'kepler':
                off, grp = KEPLER_TO_BTJD, (hdr.get('MODULE'), hdr.get('OUTPUT'))
            else:
                off, grp = 0.0, (hdr.get('CAMERA'), hdr.get('CCD'))
        t, f, q = np.asarray(t, float), np.asarray(f, float), np.asarray(q, int)
        bkg, x, y = np.asarray(bkg, float), np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    meta = {'cam': grp[0], 'ccd': grp[1]}
    return t[m] + off, f[m], bkg[m], x[m], y[m], meta


def read_qlp(path):
    return read_lc(path, 'qlp')


def detrend(t, f):
    """Normalise and remove residual trends with a running median per continuous segment."""
    f = f / np.nanmedian(f) - 1.0
    out = np.full_like(f, np.nan)
    edges = np.where(np.diff(t) > GAP_D)[0] + 1
    segs = np.split(np.arange(t.size), edges)
    dt = np.median(np.diff(t))
    w = max(5, int(DETREND_WINDOW_D / dt) | 1)
    for s in segs:
        if s.size < 20:
            continue
        trend = median_filter(f[s], size=min(w, s.size - (1 - s.size % 2)), mode='nearest')
        out[s] = f[s] - trend
    sig = mad_std(out)
    out[out > 4 * sig] = np.nan          # clip upward outliers only
    return out, segs


def to_grid(t, *arrays):
    dt = np.median(np.diff(t))
    idx = np.rint((t - t[0]) / dt).astype(int)
    n = idx[-1] + 1
    grids = []
    for a in arrays:
        g = np.full(n, np.nan)
        g[idx] = a
        grids.append(g)
    tg = t[0] + np.arange(n) * dt
    return tg, dt, grids


def box_scan(y, k):
    """In-box mean, local flank mean (k/2 each side) and count for every box start."""
    v = np.isfinite(y)
    cs = np.concatenate([[0], np.cumsum(np.where(v, y, 0.0))])
    cn = np.concatenate([[0], np.cumsum(v)])
    n = y.size
    i = np.arange(0, n - k + 1)
    s_in = cs[i + k] - cs[i]
    n_in = cn[i + k] - cn[i]
    g, h = max(1, k // 2), k
    a0, a1 = np.clip(i - g - h, 0, n), np.clip(i - g, 0, n)
    b0, b1 = np.clip(i + k + g, 0, n), np.clip(i + k + g + h, 0, n)
    s_fl = (cs[a1] - cs[a0]) + (cs[b1] - cs[b0])
    n_fl = (cn[a1] - cn[a0]) + (cn[b1] - cn[b0])
    with np.errstate(invalid='ignore', divide='ignore'):
        m_in = s_in / n_in
        m_fl = s_fl / n_fl
    ok = (n_in >= 0.7 * k) & (n_fl >= 0.5 * h)
    return i, m_in, m_fl, n_in, ok


def single_events(tg, dt, y, sigma, invert=False, interior=None):
    yy = -y if invert else y
    best = []
    for w in WIDTHS_H:
        k = max(3, int(round(w / 24.0 / dt)))
        if k * 4 > np.isfinite(yy).sum():
            continue
        i, m_in, m_fl, n_in, ok = box_scan(yy, k)
        depth = m_fl - m_in
        snr_w = depth * np.sqrt(n_in) / sigma
        red = mad_std(np.where(ok, m_in, np.nan))
        snr_r = depth / red if red and np.isfinite(red) else snr_w
        snr = np.where(ok, np.minimum(snr_w, snr_r), -np.inf)
        best.append((w, k, i, depth, snr, n_in))
    if invert:  # compare like with like: only boxes away from gaps and segment edges
        vals = [float(np.max(np.where(interior[i] & interior[np.minimum(i + k - 1, tg.size - 1)], snr, -np.inf)))
                for w, k, i, depth, snr, n_in in best]
        return max([v for v in vals if np.isfinite(v)], default=0.0)
    events = []
    masked = np.zeros(tg.size, bool)
    for _ in range(MAX_EVENTS):
        top = None
        for w, k, i, depth, snr, n_in in best:
            s = np.where(masked[i] | masked[np.minimum(i + k - 1, tg.size - 1)], -np.inf, snr)
            j = int(np.argmax(s))
            if s[j] >= SNR_EVENT and (top is None or s[j] > top[0]):
                top = (float(s[j]), w, k, i[j], float(depth[j]), int(n_in[j]))
        if top is None:
            break
        snr, w, k, j, depth, n_in = top
        events.append({'snr': snr, 'dur_h': w, 'k': k, 'j': j, 'depth': depth, 'n_in': n_in})
        pad = int(0.5 / dt) + k
        masked[max(0, j - pad):j + k + pad] = True
    return events


def characterise(ev, tg, dt, y, bkg, cx, cy, segs_t):
    j, k = ev['j'], ev['k']
    tc = tg[j] + (k - 1) * dt / 2
    sl = slice(j, j + k)
    fl = np.r_[max(0, j - 2 * k):j, j + k:min(tg.size, j + 3 * k)]
    flags = []
    # shape: outer quarters vs central half of the box
    q = max(1, k // 4)
    centre = np.nanmean(y[j + q:j + k - q]) if k - 2 * q > 0 else np.nan
    outer = np.nanmean(np.r_[y[j:j + q], y[j + k - q:j + k]])
    base = np.nanmean(y[fl]) if fl.size else 0.0
    dc, do = base - centre, base - outer
    vratio = float(do / dc) if dc > 0 else np.nan
    if np.isfinite(vratio) and vratio < 0.35:
        flags.append('vshape')
    # distance to data gaps / segment edges
    edge = min(min(abs(tc - a), abs(tc - b)) for a, b in segs_t)
    if edge < ev['dur_h'] / 48 + 0.15:
        flags.append('edge')
    # background excursion (scattered light)
    for name, arr in (('bkg', bkg), ('centroid_x', cx), ('centroid_y', cy)):
        a_in, a_fl = arr[sl], arr[fl]
        s = mad_std(arr)
        if np.isfinite(s) and s > 0 and np.isfinite(a_fl).sum() > 3:
            z = (np.nanmean(a_in) - np.nanmean(a_fl)) / (s / np.sqrt(max(1, np.isfinite(a_in).sum())))
            if abs(z) > 6:
                flags.append(name)
    return {'t0': float(tc), 'dur_h': ev['dur_h'], 'depth_ppm': ev['depth'] * 1e6, 'snr': ev['snr'],
            'n_in': ev['n_in'], 'vratio': vratio, 'flags': ','.join(flags)}


def run_bls(t, y, span):
    from astropy.timeseries import BoxLeastSquares
    m = np.isfinite(y)
    t, y = t[m], y[m]
    # 20-minute bins
    b = np.floor((t - t[0]) / (20 / 1440)).astype(int)
    nb = np.bincount(b)
    keep = nb > 0
    tb = (np.bincount(b, t) / np.maximum(nb, 1))[keep]
    yb = (np.bincount(b, y) / np.maximum(nb, 1))[keep]
    if tb.size < 200:
        return None
    pmax = min(30.0, span / 3.0) if span > 40 else min(15.0, span / 2.0)
    if pmax <= 0.6:
        return None
    durs = np.array([1, 2, 3, 4, 6]) / 24.0
    model = BoxLeastSquares(tb, yb, dy=np.full(yb.size, mad_std(yb)))
    periods = model.autoperiod(durs, minimum_period=0.5, maximum_period=pmax, frequency_factor=3.0)
    res = model.power(periods, durs, objective='snr')
    p = res.power
    i = int(np.nanargmax(p))
    sde = (p[i] - np.nanmedian(p)) / np.nanstd(p)
    P, t0, dur = res.period[i], res.transit_time[i], res.duration[i]
    st = model.compute_stats(P, dur, t0)
    ntr = int(np.sum(st['per_transit_count'] > 0))
    do, de = st['depth_odd'], st['depth_even']
    oe = abs(do[0] - de[0]) / np.hypot(do[1], de[1]) if do[1] > 0 and de[1] > 0 else 0.0
    return {'period': float(P), 't0': float(t0), 'dur_h': float(dur * 24), 'depth_ppm': float(res.depth[i] * 1e6),
            'snr': float(res.depth_snr[i]), 'sde': float(sde), 'ntr': ntr, 'oddeven_sig': float(oe),
            'sec_depth_ppm': float(st['depth_phased'][0] * 1e6)}


def process(path, tic, sector, fmt='qlp', do_bls=True, brightenings=False):
    """Full per-light-curve analysis. Returns dict or {'error': ...}. With `brightenings`, the light curve is also searched
    for microlensing events (rogues.py) and the result carries them under 'brightenings'."""
    try:
        t, f, bkg, cx, cy, meta = read_lc(path, fmt)
        o = np.argsort(t); t, f, bkg, cx, cy = t[o], f[o], bkg[o], cx[o], cy[o]
        if t.size < 300:
            return {'tic': tic, 'sector': sector, 'error': 'too few points'}
        y, segs = detrend(t, f)
        segs_t = [(t[s[0]], t[s[-1]]) for s in segs if s.size]
        tg, dt, (yg, bg, xg, yg2) = to_grid(t, y, bkg, cx, cy)
        sigma = mad_std(yg)
        evs = single_events(tg, dt, yg, sigma)
        interior = np.ones(tg.size, bool)
        for a, b in segs_t:
            interior &= ~(((tg > a - 0.01) & (tg < a + 0.3)) | ((tg > b - 0.3) & (tg < b + 0.01)))
        inv = single_events(tg, dt, yg, sigma, invert=True, interior=interior)
        events = [characterise(e, tg, dt, yg, bg, xg, yg2, segs_t) for e in evs]
        span = t[-1] - t[0]
        bls = run_bls(t, y, span) if do_bls else None
        # 1-hour noise
        k1 = max(1, int(round(1 / 24 / dt)))
        _, m_in, _, _, ok = box_scan(yg, k1)
        noise1h = mad_std(np.where(ok, m_in, np.nan))
        res = None
        if brightenings:
            from . import rogues
            try:
                res = rogues.search_lc(t, f, bkg, cx, cy, rogues.flagged_times(path, fmt))
            except Exception:
                res = []
        out = {'tic': tic, 'sector': sector, 'cam': meta['cam'], 'ccd': meta['ccd'],
                'npts': int(t.size), 't_start': float(t[0]), 't_end': float(t[-1]),
                'noise_ppm': float(sigma * 1e6), 'noise1h_ppm': float(noise1h * 1e6), 'inv_snr': float(inv),
                'events': events, 'bls': bls}
        if res is not None:
            out['brightenings'] = res
        return out
    except Exception as e:  # keep the run alive on corrupt files
        return {'tic': tic, 'sector': sector, 'error': repr(e)[:200]}
