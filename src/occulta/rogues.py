"""Free-floating (rogue) planets: single, symmetric brightenings of a star with the shape of a point-lens
microlensing event (Paczynski 1986), lasting from about an hour to about a day.

    occulta rogues search TIC [--skybot] [--finite-source]
                                  every light curve of one star, downloaded from MAST and searched for brightenings
    occulta rogues scan [--minutes M] [--min-inv-snr X] [--limit N] [--skybot]
                                  re-search the light curves that the survey has already processed and whose strongest
                                  brightening (inv_snr in the database) reaches X (default 7); resumable
    occulta rogues list [--skybot]
                                  rank every brightening found so far (findings/rogues/rogues.csv and report.md)

The survey searches each new light curve for brightenings itself when "search_brightenings": true is set in
config.json. Every brightening is fitted with a flat baseline, a point-lens model (t0, tE, u0 and the blending
fraction fs; a uniform finite source with --finite-source) and the flare template of Davenport et al. (2014), and is
rejected as a flare (the template fits better, or the profile is not mirror-symmetric), an asteroid (the centroid moves
during the event, the background changes, or SkyBoT lists a known minor body), a quality-flagged cadence (momentum dump,
scattered light), a repeating variable, an event at a gap, or a detector systematic (brightenings of many stars on the
same camera and CCD at the same time).
"""
import argparse, json, time
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.stats import chi2 as _chi2, norm

from . import lightcurves as lcs

SNR_TRIGGER = 7.0        # box S/N of a brightening in the trigger series (as for dips)
TRIGGER_WINDOW_D = 4.0   # running median of the trigger series: wider than the 2 d of the dip search, so a 1-d event survives
TE_MIN_D, TE_MAX_D = 1 / 24, 1.5   # Einstein time accepted for a planetary-mass lens (d)
DCHI2_FLAT = 50.0        # the lens must beat a flat baseline by this delta chi2
DCHI2_FLARE = 9.0        # ... and the flare template by this delta chi2 (3 sigma)
ASYM_SIGMA = 3.0         # mirror asymmetry about t0 above this rejects the event (flare: fast rise, slow decay)
MOTION_SIGMA = 5.0       # in-event centroid trend (a moving source) or background change above this rejects the event
BIN_MIN = 10 / 1440      # light curves with shorter cadences are binned to 10 min for the fits
COMMON_DT_D = 0.1        # brightenings of this many stars (config common_mode_n) on one CCD within this time are systematics
SS = (np.arange(5) - 2) / 5   # sub-exposure offsets for the integration over the exposure time


# ---------------------------------------------------------------- models
def magnification(u, rho=0.0):
    """Point-lens magnification A(u) = (u^2 + 2) / (u sqrt(u^2 + 4)) (Paczynski 1986); for rho > 0 the magnification of
    a uniform source disk of radius rho (in Einstein radii; Witt & Mao 1994), integrated over circles about the lens:
    A = (1 / pi rho^2) Int A(r) r theta(r) dr, theta(r) the angle of the circle of radius r that lies on the disk, so that
    the integrand is finite at r = 0."""
    u = np.abs(np.asarray(u, float))
    if rho <= 0:
        u = np.maximum(u, 1e-9)
        return (u * u + 2) / (u * np.sqrt(u * u + 4))
    x, w = np.polynomial.legendre.leggauss(48)
    a, b, c = np.abs(u - rho), u + rho, np.maximum(u - rho, 0.0)
    total = 0.0
    for lo, hi in ((c, a), (a, b)):           # the integrand has a kink at |u - rho|
        r = lo[..., None] + (hi - lo)[..., None] * (x + 1) / 2
        cosang = (r * r + u[..., None] ** 2 - rho * rho) / np.maximum(2 * r * u[..., None], 1e-300)
        ang = 2 * np.arccos(np.clip(np.where(r < rho - u[..., None], -1.0, cosang), -1, 1))
        total = total + ((hi - lo) / 2) * np.sum(w * (r * r + 2) / np.sqrt(r * r + 4) * ang, -1)
    return total / (np.pi * rho * rho)


def lens(t, t0, te, u0, fs, rho=0.0):
    """Relative flux of a microlensing event, F/F_base - 1 = fs (A(u) - 1), u^2 = u0^2 + ((t - t0)/tE)^2; fs is the
    fraction of the baseline flux that comes from the lensed source (1 - fs from blended stars)."""
    return fs * (magnification(np.hypot(u0, (t - t0) / te), rho) - 1)


def flare(t, tpk, fwhm, amp):
    """Flare template of Davenport et al. (2014): fourth-order polynomial rise over one FWHM, double-exponential decay."""
    x = (t - tpk) / fwhm
    rise = 1 + 1.941 * x - 0.175 * x ** 2 - 2.246 * x ** 3 - 1.125 * x ** 4
    out = np.where((x >= -1) & (x <= 0), rise, 0.0)
    xp = np.clip(x, 0, None)
    out = np.where(x > 0, 0.6890 * np.exp(-1.600 * xp) + 0.3030 * np.exp(-0.2783 * xp), out)
    return amp * np.clip(out, 0, None)


def _integrated(f, t, exp, *p):
    return sum(f(t + o * exp, *p) for o in SS) / len(SS)


# ---------------------------------------------------------------- tests of one event
def mirror_asymmetry(t, r, t0, half, sig):
    """Significance (sigma) of the difference between the flux after t0 and the mirrored flux before t0, within `half`
    of t0; r is the flux above the baseline. A symmetric event gives about 0, a flare a large positive value."""
    pre = (t < t0) & (t > t0 - half)
    post = (t > t0) & (t < t0 + half)
    if pre.sum() < 3 or post.sum() < 3:
        return np.nan
    tp, rp = (2 * t0 - t[pre])[::-1], r[pre][::-1]
    m = (t[post] >= tp[0]) & (t[post] <= tp[-1])
    if m.sum() < 3:
        return np.nan
    d = r[post][m] - np.interp(t[post][m], tp, rp)
    c = np.sum(d ** 2) / (2 * sig ** 2)            # each difference carries the noise of two points
    return float(norm.isf(max(_chi2.sf(c, m.sum()), 1e-300)))


def motion(t, a, t0, half):
    """(trend, shift): significance of a linear trend of a centroid (or background) series inside the event, and of its
    mean shift there, both against a straight line through the data outside the event. A moving asteroid drags the
    centroid across the aperture (trend); a brightening on a neighbour displaces it (shift)."""
    ok = np.isfinite(a)
    t, a = t[ok], a[ok]
    inn = np.abs(t - t0) < half
    if inn.sum() < 4 or (~inn).sum() < 6:
        return 0.0, 0.0
    c = np.polyfit(t[~inn] - t0, a[~inn], 1)
    r = a - np.polyval(c, t - t0)
    s = lcs.mad_std(np.diff(a)) / np.sqrt(2)        # point-to-point scatter, insensitive to slow drifts
    if not np.isfinite(s) or s <= 0:
        return 0.0, 0.0
    ti, ri = t[inn] - t0, r[inn]
    shift = np.mean(ri) / (s / np.sqrt(ri.size))
    sxx = np.sum((ti - ti.mean()) ** 2)
    trend = (np.sum((ti - ti.mean()) * ri) / sxx) / (s / np.sqrt(sxx)) if sxx > 0 else 0.0
    return float(abs(trend)), float(abs(shift))


def _bin(t, *arrays, width):
    b = np.floor((t - t[0]) / width).astype(int)
    n = np.bincount(b)
    k = n > 0
    out = [np.bincount(b, t)[k] / n[k]]
    for a in arrays:
        v = np.isfinite(a)
        s, m = np.bincount(b, np.where(v, a, 0.0)), np.bincount(b, v.astype(float))
        with np.errstate(invalid='ignore', divide='ignore'):
            out.append((s / m)[k])
    return out


def fit_event(t, y, sig, tc, w, exp, rho=None):
    """Fit a flat baseline, a point lens and a flare to one window (relative flux y; each model has a linear baseline).
    Returns a dict with the lens parameters and the chi2 of the three models."""
    x = t - tc
    amp = max(np.nanmax(median_filter(y, 3)) - np.median(y), 3 * sig)
    A = np.vstack([np.ones_like(x), x]).T
    cf, *_ = np.linalg.lstsq(A, y, rcond=None)
    chi_flat = float(np.sum(((y - A @ cf) / sig) ** 2))
    span = max(w, 0.1)
    rho_on = rho is not None

    def r_lens(p):
        m = _integrated(lens, x, exp, p[0], 10 ** p[1], p[2], p[3], p[6] if rho_on else 0.0)
        return (y - m - p[4] - p[5] * x) / sig
    lo = [-span, np.log10(0.5 / 24), 1e-3, 1e-3, -1, -1] + ([1e-3] if rho_on else [])
    hi = [span, np.log10(3.0), 3.0, 1.0, 1, 1] + ([3.0] if rho_on else [])
    best = None
    for te in (w / 2, w, 2 * w):
        for u0 in (0.1, 0.5, 1.5):
            fs = float(np.clip(amp / (magnification(u0) - 1), 2e-3, 0.999))
            p0 = np.clip([0, np.log10(max(te, 0.6 / 24)), u0, fs, cf[0], 0] + ([0.05] if rho_on else []),
                         np.array(lo) + 1e-6, np.array(hi) - 1e-6)
            f = least_squares(r_lens, p0, bounds=(lo, hi), x_scale='jac')
            if best is None or f.cost < best.cost:
                best = f
    pl = best.x

    def r_flare(p):
        return (y - _integrated(flare, x, exp, p[0], 10 ** p[1], p[2]) - p[3] - p[4] * x) / sig
    lo_f, hi_f = [-span, np.log10(1 / 1440), 0, -1, -1], [span, 0.0, 50 * amp + 1e-3, 1, 1]
    xpk = float(x[np.nanargmax(median_filter(y, 3))])
    bestf = None
    for fw in (w / 8, w / 3, w):
        p0 = np.clip([np.clip(xpk, -span, span), np.log10(max(fw, 2 / 1440)), amp, cf[0], 0],
                     np.array(lo_f) + 1e-6, np.array(hi_f) - 1e-6)
        f = least_squares(r_flare, p0, bounds=(lo_f, hi_f), x_scale='jac')
        if bestf is None or f.cost < bestf.cost:
            bestf = f
    chi_lens, chi_flare = 2 * best.cost, 2 * bestf.cost
    t0, te, u0, fs = tc + pl[0], 10 ** pl[1], pl[2], pl[3]
    rr = pl[6] if rho_on else 0.0
    return {'t0': float(t0), 'te_d': float(te), 'u0': float(u0), 'fs': float(fs), 'rho': float(rr),
            'amp_ppm': float(lens(t0, t0, te, u0, fs, rr) * 1e6), 'base': (float(pl[4]), float(pl[5])),
            'chi2_flat': chi_flat, 'chi2_lens': float(chi_lens), 'chi2_flare': float(chi_flare),
            'dchi2_flat': float(chi_flat - chi_lens), 'dchi2_flare': float(chi_flare - chi_lens),
            'chi2r': float(chi_lens / max(1, y.size - len(pl)))}


def _segments(t):
    edges = np.where(np.diff(t) > lcs.GAP_D)[0] + 1
    return [s for s in np.split(np.arange(t.size), edges) if s.size]


def triggers(t, f):
    """Brightenings in a light curve: the dip search of the survey run on the inverted series, after a 4-day running
    median per segment and without the clipping of upward points. Returns ([(tc, box width d, S/N)], sigma)."""
    y = np.full_like(f, np.nan)
    dt = np.median(np.diff(t))
    w = max(5, int(TRIGGER_WINDOW_D / dt) | 1)
    for s in _segments(t):
        if s.size >= 20:
            y[s] = f[s] - median_filter(f[s], size=min(w, s.size - (1 - s.size % 2)), mode='nearest')
    tg, dtg, (yg,) = lcs.to_grid(t, y)
    sigma = lcs.mad_std(yg)
    ev = lcs.single_events(tg, dtg, -yg, sigma) if np.isfinite(sigma) and sigma > 0 else []
    return [(float(tg[e['j']] + (e['k'] - 1) * dtg / 2), e['dur_h'] / 24, float(e['snr'])) for e in ev], sigma


def search_lc(t, f, bkg=None, cx=None, cy=None, flagged=(), finite_source=False, keep_windows=True):
    """Search one light curve (time in BTJD, flux, background, centroid columns, times of quality-flagged cadences)
    for single, symmetric brightenings. Returns a list of dicts, one per brightening, with the fit and its flags;
    the flags name every test that the brightening fails (empty: a microlensing candidate on this light curve)."""
    t, f = np.asarray(t, float), np.asarray(f, float)
    n = t.size
    bkg, cx, cy = (np.full(n, np.nan) if a is None else np.asarray(a, float) for a in (bkg, cx, cy))
    ok = np.isfinite(t) & np.isfinite(f)
    o = np.argsort(t[ok])
    t, f, bkg, cx, cy = (a[ok][o] for a in (t, f, bkg, cx, cy))
    if t.size < 100:
        return []
    f = f / np.nanmedian(f) - 1.0
    trig, _ = triggers(t, f)
    flagged = np.sort(np.asarray(flagged, float))
    segs = [(t[s[0]], t[s[-1]]) for s in _segments(t)]
    cad = float(np.median(np.diff(t)))
    out = []
    for tc, w, snr in trig:
        seg = next(((a, b) for a, b in segs if a - 1e-6 <= tc <= b + 1e-6), None)
        flags = []
        if seg is None or min(tc - seg[0], seg[1] - tc) < 0.3:
            continue                 # gap and segment edges: the same exclusion as for the false-alarm level of the dips
        W = float(np.clip(5 * w, 0.5, 2.0))
        m = (np.abs(t - tc) < W) & (t >= seg[0]) & (t <= seg[1])
        tw, fw, bw, xw, yw = t[m], f[m], bkg[m], cx[m], cy[m]
        if tw.size < 30:
            continue
        exp = cad
        if cad < BIN_MIN * 0.9:
            tw, fw, bw, xw, yw = _bin(tw, fw, bw, xw, yw, width=BIN_MIN)
            exp = BIN_MIN
        sig = lcs.mad_std(np.diff(fw)) / np.sqrt(2)
        if not np.isfinite(sig) or sig <= 0:
            continue
        r = fit_event(tw, fw, sig, tc, w, exp, rho=0.05 if finite_source else None)
        t0, te = r['t0'], r['te_d']
        model = lens(tw, t0, te, r['u0'], r['fs'], r['rho'])
        pk = model.max() if model.size else 0.0
        ext = tw[model >= 0.1 * pk] if pk > 0 else tw[:0]
        half10 = float(max(np.max(np.abs(ext - t0)) if ext.size else 0.0, 2 * exp))
        fwhm_half = float(max(np.max(np.abs(tw[model >= 0.5 * pk] - t0)) if pk > 0 and np.any(model >= 0.5 * pk) else 0.0, exp))
        resid = fw - r['base'][0] - r['base'][1] * (tw - tc)
        r['asym_sigma'] = mirror_asymmetry(tw, resid, t0, half10, sig)
        mt = [motion(tw, a, t0, half10) for a in (xw, yw)]
        r['motion_sigma'] = max(m_[0] for m_ in mt)
        r['offset_sigma'] = max(m_[1] for m_ in mt)
        r['bkg_sigma'] = max(motion(tw, bw, t0, half10))
        r['snr'] = snr
        if not (TE_MIN_D <= te <= TE_MAX_D): flags.append('duration')
        if r['dchi2_flat'] < DCHI2_FLAT: flags.append('weak')
        if r['dchi2_flare'] < DCHI2_FLARE: flags.append('flare')
        if np.isfinite(r['asym_sigma']) and r['asym_sigma'] > ASYM_SIGMA: flags.append('asymmetric')
        if r['motion_sigma'] > MOTION_SIGMA: flags.append('moving')
        if r['bkg_sigma'] > MOTION_SIGMA: flags.append('bkg')
        if flagged.size and np.any(np.abs(flagged - t0) < fwhm_half + exp): flags.append('quality')
        cov = np.sum(np.abs(tw - t0) < fwhm_half) * exp / (2 * fwhm_half)
        if cov < 0.6 or np.sum((tw > t0 - half10) & (tw < t0)) < 3 or np.sum((tw > t0) & (tw < t0 + half10)) < 3:
            flags.append('gap')
        if any(abs(tc2 - t0) > max(3 * half10, 0.5) and s2 >= max(SNR_TRIGGER, 0.5 * snr) for tc2, _, s2 in trig):
            flags.append('repeats')
        r['flags'] = ','.join(flags)
        if keep_windows and not flags:
            r['window'] = {'t': tw, 'f': resid, 'model': model}
        r.pop('base')
        out.append(r)
    return out


def flagged_times(path, fmt):
    """Times (BTJD) of the cadences that the pipeline flagged (QUALITY != 0): momentum dumps, scattered light,
    Argabrightening and the like. Empty when the file has no quality column."""
    from astropy.io import fits
    try:
        with fits.open(path, memmap=False) as h:
            if fmt == 'k2sff':
                return np.array([])
            d = h[1].data
            q = d['QUALITY'] if 'QUALITY' in d.columns.names else d['SAP_QUALITY']
            t = np.asarray(d['TIME'], float) + (lcs.KEPLER_TO_BTJD if fmt == 'kepler' else 0.0)
            q = np.asarray(q, int)
            return t[np.isfinite(t) & (q != 0)]
    except Exception:
        return np.array([])


def search_file(path, fmt, finite_source=False):
    """(events, cam, ccd) of one light-curve file."""
    t, f, bkg, cx, cy, meta = lcs.read_lc(path, fmt)
    return search_lc(t, f, bkg, cx, cy, flagged_times(path, fmt), finite_source), meta['cam'], meta['ccd']


# ---------------------------------------------------------------- database, ranking and reports
COLS = ('t0', 'te_d', 'u0', 'fs', 'amp_ppm', 'snr', 'dchi2_flat', 'dchi2_flare', 'asym_sigma', 'motion_sigma',
        'offset_sigma', 'bkg_sigma', 'chi2r', 'flags')


def store(con, tic, unit, cam, ccd, events, findings):
    """Keep the brightenings of one light curve; the windows of unflagged ones go to findings/rogues/<star>/."""
    from . import sources
    con.execute('DELETE FROM brightening WHERE tic=? AND sector=?', (tic, unit))
    for e in events:
        con.execute(f'INSERT OR REPLACE INTO brightening VALUES (?,?,?,?,{",".join("?" * len(COLS))},?)',
                    (tic, unit, cam, ccd, *[e[k] if not isinstance(e[k], float) or np.isfinite(e[k]) else None for k in COLS], None))
        if 'window' in e:
            d = Path(findings) / 'rogues' / sources.name_of(tic).replace(' ', '')
            d.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(d / f"u{unit}_{e['t0']:.3f}.npz", **e['window'], te=e['te_d'], u0=e['u0'], fs=e['fs'])


def minor_bodies(ra, dec, t_btjd, radius_arcmin=2.0):
    """Known asteroids and comets within the radius at the time of the event (SkyBoT; Berthier et al. 2006), as a list
    of names; None when the service cannot be reached."""
    try:
        from astropy.coordinates import SkyCoord
        from astropy.time import Time
        import astropy.units as un
        from astroquery.imcce import Skybot
        tab = Skybot.cone_search(SkyCoord(ra * un.deg, dec * un.deg), radius_arcmin * un.arcmin,
                                 Time(t_btjd + 2457000.0, format='jd'), location='C57')
        return [str(x) for x in tab['Name']]
    except Exception as e:
        return [] if 'No solar system object' in str(e) else None


def rank(con, findings, common_n=5, skybot=False, star=None):
    """Mark systematics common to one CCD, rank the brightenings, and write findings/rogues/rogues.csv, report.md and,
    for every candidate, <star>/report.md, summary.json and event.png. Returns the list of candidates."""
    from . import sources
    rows = con.execute('SELECT rowid, sector, cam, ccd, t0 FROM brightening').fetchall()
    by = {}
    for rid, s, cam, ccd, t0 in rows:
        by.setdefault((s, cam, ccd), []).append((t0, rid))
    for lst in by.values():
        ts = np.sort([x[0] for x in lst])
        for t0, rid in lst:
            n = np.searchsorted(ts, t0 + COMMON_DT_D) - np.searchsorted(ts, t0 - COMMON_DT_D)
            if n >= common_n:
                con.execute("UPDATE brightening SET flags = CASE WHEN flags='' THEN 'common' ELSE flags||',common' END "
                            "WHERE rowid=? AND instr(flags, 'common')=0", (rid,))
    cands = []
    for rid, tic, s, t0, te, u0, fs, amp, d1, d2 in con.execute(
            "SELECT rowid, tic, sector, t0, te_d, u0, fs, amp_ppm, dchi2_flat, dchi2_flare FROM brightening "
            "WHERE flags='' ORDER BY dchi2_flat DESC").fetchall():
        status, note = 'lens_candidate', ''
        if skybot:
            st = (star or _star)(con, tic)
            mb = minor_bodies(st[0], st[1], t0) if st and st[0] is not None else None
            if mb:
                status, note = 'rejected', 'known minor body: ' + ', '.join(mb[:5])
                con.execute("UPDATE brightening SET flags='minor_body' WHERE rowid=?", (rid,))
            elif mb is None:
                note = 'SkyBoT could not be reached'
        con.execute('UPDATE brightening SET status=? WHERE rowid=?', (status, rid))
        if status == 'lens_candidate':
            cands.append({'tic': tic, 'unit': s, 'unit_name': sources.unit_name(s), 't0_btjd': t0, 'te_d': te, 'u0': u0,
                          'fs': fs, 'amp_ppm': amp, 'dchi2_flat': d1, 'dchi2_flare': d2, 'note': note})
    con.execute("UPDATE brightening SET status='rejected' WHERE flags!=''")
    con.commit()
    out = Path(findings) / 'rogues'
    out.mkdir(parents=True, exist_ok=True)
    allrows = con.execute(f'SELECT tic, sector, cam, ccd, {", ".join(COLS)}, status FROM brightening '
                          'ORDER BY status="lens_candidate" DESC, dchi2_flat DESC').fetchall()
    head = ['tic', 'unit', 'cam', 'ccd', *COLS, 'status']
    import csv
    with open(out / 'rogues.csv', 'w', newline='', encoding='utf-8') as fh:
        wr = csv.writer(fh); wr.writerow(head); wr.writerows(allrows)
    reasons = {}
    for (fl,) in con.execute("SELECT flags FROM brightening WHERE flags!=''"):
        for k in fl.split(','):
            reasons[k] = reasons.get(k, 0) + 1
    L = ['# Free-floating planets: microlensing candidates', '',
         f'{len(allrows)} brightenings fitted, {len(cands)} pass every test. Rejections by test (an event can fail '
         'several): ' + (', '.join(f'{k} {v}' for k, v in sorted(reasons.items(), key=lambda kv: -kv[1])) or 'none') + '.', '',
         '| Rank | Star | Unit | t0 (BTJD) | tE (h) | u0 | fs | Peak (ppm) | dchi2 flat | dchi2 flare |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for i, c in enumerate(cands, 1):
        L.append(f"| {i} | {sources.name_of(c['tic'])} | {c['unit_name']} | {c['t0_btjd']:.4f} | {c['te_d'] * 24:.2f} | "
                 f"{c['u0']:.3f} | {c['fs']:.3f} | {c['amp_ppm']:.0f} | {c['dchi2_flat']:.0f} | {c['dchi2_flare']:.0f} |")
        _candidate_report(out / sources.name_of(c['tic']).replace(' ', ''), c)
    (out / 'report.md').write_text('\n'.join(L) + '\n', encoding='utf-8')
    return cands


def _star(con, tic):
    from . import planets
    return planets.get_star(con, tic)


def _candidate_report(d, c):
    from . import sources
    d.mkdir(parents=True, exist_ok=True)
    (d / 'summary.json').write_text(json.dumps(c, indent=1, default=float))
    L = [f"# {sources.name_of(c['tic'])} — LENS_CANDIDATE", '',
         f"Single brightening in {c['unit_name']} at BTJD {c['t0_btjd']:.4f}: Einstein time {c['te_d'] * 24:.2f} h, impact "
         f"parameter u0 {c['u0']:.3f}, source fraction fs {c['fs']:.3f}, peak {c['amp_ppm']:.0f} ppm.", '',
         f"The point-lens model fits better than a flat baseline by delta chi2 {c['dchi2_flat']:.0f} and better than the "
         f"flare template by delta chi2 {c['dchi2_flare']:.0f}; the event is mirror-symmetric, the centroid and the "
         "background do not move during it, no flagged cadence falls inside it, and no other star on the same CCD "
         "brightens at the same time.", '']
    if c.get('note'):
        L += ['Note: ' + c['note'], '']
    win = sorted(d.glob(f"u{c['unit']}_{c['t0_btjd']:.3f}.npz"))
    if win:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        z = np.load(win[0])
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.plot((z['t'] - c['t0_btjd']) * 24, z['f'] * 1e6, 'o', ms=2.5, color='0.15')
        ax.plot((z['t'] - c['t0_btjd']) * 24, z['model'] * 1e6, color='C3', lw=2)
        ax.set_xlabel('hours from the peak'); ax.set_ylabel('ppm above the baseline')
        ax.set_title(f"{sources.name_of(c['tic'])}  {c['unit_name']}  tE {c['te_d'] * 24:.2f} h  u0 {c['u0']:.2f}", fontsize=10)
        fig.tight_layout(); fig.savefig(d / 'event.png', dpi=110); plt.close(fig)
        L += ['![event](event.png)', '']
    (d / 'report.md').write_text('\n'.join(L), encoding='utf-8')


# ---------------------------------------------------------------- entry points
UNIT_OFFSET = {'QLP': 0, 'SPOC': 1000, 'TESS-SPOC': 4000, 'K2SFF': 2000, 'Kepler': 3000}


def search_star(con, tic, findings, finite_source=False):
    """Download every light curve of one star and search it; returns the number of brightenings kept."""
    import shutil
    from . import planets
    work = Path(findings).parent / 'data' / 'cache' / f'rogue_{tic}'
    work.mkdir(parents=True, exist_ok=True)
    n = 0
    try:
        for author, num, path in planets.fetch_all(tic, work):
            unit = num % 1000 + UNIT_OFFSET[author]
            ev, cam, ccd = search_file(path, planets.AUTHORS[author], finite_source)
            store(con, tic, unit, cam, ccd, ev, findings)
            con.execute('INSERT OR REPLACE INTO rogue_scan VALUES (?,?,1)', (tic, unit))
            n += len(ev)
        con.commit()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return n


def scan(con, findings, minutes=55, min_inv_snr=SNR_TRIGGER, limit=None, finite_source=False, log=print):
    """Re-search the star-sectors that the survey has processed and whose strongest brightening (inv_snr) reaches
    `min_inv_snr`, strongest first, one light curve at a time; resumable through the table rogue_scan."""
    import os
    from . import sources, survey
    deadline = time.time() + minutes * 60
    q = ('SELECT ss.tic, ss.sector, t.url FROM star_sector ss LEFT JOIN target t ON t.tic=ss.tic AND t.sector=ss.sector '
         'WHERE ss.inv_snr >= ? AND NOT EXISTS (SELECT 1 FROM rogue_scan r WHERE r.tic=ss.tic AND r.sector=ss.sector) '
         'ORDER BY ss.inv_snr DESC') + (f' LIMIT {int(limit)}' if limit else '')
    todo = con.execute(q, (min_inv_snr,)).fetchall()
    log(f'rogues: {len(todo)} processed light curves with inv_snr >= {min_inv_snr} to re-search')
    survey.RAW.mkdir(parents=True, exist_ok=True)
    done = 0
    for tic, u, url in todo:
        if time.time() > deadline:
            break
        path = survey.fetch(u, tic, url)
        status = 3
        if path:
            try:
                ev, cam, ccd = search_file(path, sources.fmt_of_unit(u), finite_source)
                store(con, tic, u, cam, ccd, ev, findings)
                status = 1
            except Exception as e:
                log('rogues: search failed', tic, u, repr(e)[:150]); status = 2
            finally:
                try: os.remove(path)
                except OSError: pass
        con.execute('INSERT OR REPLACE INTO rogue_scan VALUES (?,?,?)', (tic, u, status))
        done += 1
        if done % 50 == 0:
            con.commit(); log('rogues:', done, 'light curves re-searched')
    con.commit()
    return done


def main(argv=None):
    from . import paths, survey
    ap = argparse.ArgumentParser(prog='occulta rogues', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='sub', required=True)
    s = sp.add_parser('search', help='download and search every light curve of one star')
    from .cli import tic_arg
    s.add_argument('tic', type=tic_arg)
    sc = sp.add_parser('scan', help='re-search processed light curves whose inv_snr passes a threshold')
    sc.add_argument('--minutes', type=float, default=55)
    sc.add_argument('--min-inv-snr', type=float, default=SNR_TRIGGER)
    sc.add_argument('--limit', type=int)
    sp.add_parser('list', help='rank every brightening found so far')
    for p in (s, sc, sp.choices['list']):
        p.add_argument('--skybot', action='store_true', help='check every candidate against the known minor bodies (SkyBoT, online)')
    for p in (s, sc):
        p.add_argument('--finite-source', action='store_true', help='fit a uniform finite source as well')
    a = ap.parse_args(argv)
    paths.init()
    con = survey.db()
    if a.sub == 'search':
        n = search_star(con, a.tic, survey.FIND, a.finite_source)
        print(f'{n} brightenings fitted')
    elif a.sub == 'scan':
        scan(con, survey.FIND, a.minutes, a.min_inv_snr, a.limit, a.finite_source, survey.log)
    cands = rank(con, survey.FIND, survey.CFG.get('common_mode_n', 5), a.skybot)
    print(f'{len(cands)} microlensing candidates; see', survey.FIND / 'rogues' / 'report.md')
