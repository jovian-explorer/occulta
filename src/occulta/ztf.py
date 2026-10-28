"""Zwicky Transient Facility (ZTF) light curves: public IRSA download, cleaning, flux conversion and a
sparse-cadence box-least-squares search for transits and eclipses.

A ZTF light curve is a dict of equal-length numpy arrays, the same plain-dict pattern as lightcurves.py:
    t     time in BTJD (HJD - 2457000, the time base used for all other occulta sources)
    flux  relative flux, normalised to a median of 1 within each filter
    err   1-sigma flux error
    filt  filter code per point ('zg', 'zr' or 'zi')
plus the scalars 'oid' (dict filter -> ZTF object id used) and 'n_raw'/'n_good' (point counts).
Errors raise ZtfError; searches return a plain dict (or None when the data are too sparse).
See docs/ztf_notes.md for the design and its limits.
"""
import csv
import io

import numpy as np
import requests

IRSA_URL = 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves'
MAX_RADIUS_DEG = 0.1667            # the API limit for CIRCLE/POS
BTJD_OFFSET = 2457000.0
FILTERS = ('zg', 'zr', 'zi')
MAX_MAGERR = 0.1                   # mag; points with a larger error are dropped
MIN_POINTS = 30                    # per filter, after cleaning
UNIT_OF = {'zg': 5001, 'zr': 5002, 'zi': 5003}   # unit codes of sources.py (5000 + filter id), one per filter


class ZtfError(RuntimeError):
    pass


def query_url(ra=None, dec=None, radius_deg=None, oid=None, bands=None):
    """(url, params) of an IRSA ZTF-LC-API request: by position (degrees, ICRS) or by ZTF object id."""
    if oid is not None:
        ids = [oid] if np.isscalar(oid) else list(oid)
        return IRSA_URL, {'ID': [str(int(i)) for i in ids], 'FORMAT': 'csv'}
    if ra is None or dec is None or radius_deg is None:
        raise ZtfError('give ra, dec and radius_deg, or a ZTF object id')
    if not 0 < radius_deg <= MAX_RADIUS_DEG:
        raise ZtfError(f'radius_deg must be in (0, {MAX_RADIUS_DEG}] for the IRSA ZTF API')
    if not -90 <= dec <= 90:
        raise ZtfError('dec must lie in [-90, 90]')
    p = {'POS': f'CIRCLE {float(ra) % 360.0:.6f} {float(dec):.6f} {float(radius_deg):.6f}', 'FORMAT': 'csv'}
    if bands:
        p['BANDNAME'] = ','.join(bands)
    return IRSA_URL, p


def fetch_table(ra=None, dec=None, radius_deg=None, oid=None, bands=None, timeout=120):
    """Download the raw light-curve table (text) from IRSA. Public, no login needed."""
    url, params = query_url(ra, dec, radius_deg, oid, bands)
    try:
        r = requests.get(url, params=params, timeout=timeout)
    except requests.RequestException as e:
        raise ZtfError(f'IRSA request failed: {e!r}') from e
    if r.status_code != 200:
        raise ZtfError(f'IRSA returned HTTP {r.status_code}: {r.text[:200]!r}')
    return r.text


def parse_table(text):
    """Parse an IRSA CSV or IPAC table into {column: numpy array}. Needs oid, hjd, mag, magerr, catflags, filtercode."""
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith('\\')]
    if not lines:
        raise ZtfError('empty response from IRSA')
    if lines[0].lstrip().startswith('<') or lines[0].upper().startswith('ERROR'):
        raise ZtfError(f'IRSA returned an error: {lines[0][:200]!r}')
    if lines[0].startswith('|'):                       # IPAC table: header rows start with '|'
        hdr = [l for l in lines if l.startswith('|')]
        edges = [i for i, c in enumerate(hdr[0]) if c == '|']
        spans = list(zip(edges[:-1], edges[1:]))
        names = [hdr[0][a + 1:b].strip() for a, b in spans]
        rows = [[l[a + 1:b].strip() for a, b in spans] for l in lines if not l.startswith('|')]
    else:
        parsed = list(csv.reader(io.StringIO('\n'.join(lines))))
        names, rows = [n.strip() for n in parsed[0]], parsed[1:]
    need = {'oid', 'hjd', 'mag', 'magerr', 'catflags', 'filtercode'}
    if not need <= set(names):
        raise ZtfError(f'missing columns {sorted(need - set(names))}; got {names[:12]}')
    col = {n: [r[i] for r in rows] for i, n in enumerate(names)}
    nul = ('', 'null', 'NaN')
    out = {'oid': np.array([int(v) for v in col['oid']], dtype=np.int64),
           'filtercode': np.array(col['filtercode'], dtype='U2')}
    for n in ('hjd', 'mag', 'magerr'):
        out[n] = np.array([float(v) if v not in nul else np.nan for v in col[n]])
    out['catflags'] = np.array([int(float(v)) if v not in nul else -1 for v in col['catflags']], dtype=np.int64)
    return out


def good_mask(tab, max_magerr=MAX_MAGERR):
    """Points that are finite, have catflags == 0 (no processing flag set) and magerr <= max_magerr."""
    return (tab['catflags'] == 0) & np.isfinite(tab['hjd']) & np.isfinite(tab['mag']) & np.isfinite(tab['magerr']) \
        & (tab['magerr'] > 0) & (tab['magerr'] <= max_magerr)


def mag_to_flux(mag, magerr):
    """(relative flux, flux error) of magnitudes, normalised to the median magnitude (flux 1 at the median)."""
    mag, magerr = np.asarray(mag, float), np.asarray(magerr, float)
    flux = 10.0 ** (-0.4 * (mag - np.median(mag)))
    return flux, flux * 0.4 * np.log(10.0) * magerr


def build_lightcurve(tab, max_magerr=MAX_MAGERR, min_points=MIN_POINTS):
    """Occulta-format light curve from a parsed table.

    Each filter is normalised to its own median (zero points and colour terms differ between zg, zr and zi), and
    for each filter only the ZTF object id (field/CCD) with the most good points is used, because the same star in
    two fields is two objects with different zero points. Raises ZtfError if no filter has enough points.
    """
    good = good_mask(tab, max_magerr)
    T, F, E, C, oids = [], [], [], [], {}
    for flt in FILTERS:
        m = good & (tab['filtercode'] == flt)
        if not m.any():
            continue
        ids, n = np.unique(tab['oid'][m], return_counts=True)
        best = ids[np.argmax(n)]
        m &= tab['oid'] == best
        if m.sum() < min_points:
            continue
        f, e = mag_to_flux(tab['mag'][m], tab['magerr'][m])
        T.append(tab['hjd'][m] - BTJD_OFFSET)
        F.append(f)
        E.append(e)
        C.append(np.full(int(m.sum()), flt))
        oids[flt] = int(best)
    if not T:
        raise ZtfError(f'no filter has >= {min_points} good points (catflags == 0, magerr <= {max_magerr})')
    t, f, e, c = np.concatenate(T), np.concatenate(F), np.concatenate(E), np.concatenate(C)
    o = np.argsort(t)
    return {'t': t[o], 'flux': f[o], 'err': e[o], 'filt': c[o], 'oid': oids,
            'n_raw': int(len(tab['hjd'])), 'n_good': int(good.sum())}


def fetch_lightcurve(ra=None, dec=None, radius_deg=2.0 / 3600, oid=None, bands=None, timeout=120, **kw):
    """Download, parse and clean in one call. Default radius is 2 arcsec."""
    return build_lightcurve(parse_table(fetch_table(ra, dec, radius_deg, oid, bands, timeout)), **kw)


def period_grid(baseline, pmin=0.2, pmax=30.0, dur=0.08, oversample=2.0):
    """Periods log-spaced so that the phase slip over the baseline stays below dur/oversample.

    A box of duration dur at period P has phase width dur/P; two trial frequencies f and f + df are
    distinguishable once baseline * df < (dur / P) / oversample, i.e. d(ln f) = dur / (oversample * baseline).
    """
    n = int(np.ceil(np.log(pmax / pmin) * oversample * baseline / dur)) + 1
    return pmin * np.exp(np.linspace(0.0, np.log(pmax / pmin), n))


DURATIONS = np.array([0.04, 0.06, 0.09, 0.13, 0.2, 0.3, 0.45])   # days (about 1 h to 11 h)


def search(lc, pmin=0.2, pmax=30.0, durations=DURATIONS, min_in_transit=5, min_transits=3, ncand=30, band=None):
    """Box-least-squares search of a ZTF light curve with irregular, about nightly sampling.

    Uses astropy BoxLeastSquares on the per-filter-normalised flux (all filters together unless band is given),
    over the period_grid of the baseline. Candidates are the ncand best SNR peaks per duration; those with fewer
    than min_in_transit in-box points or fewer than min_transits distinct transit epochs are rejected, which is
    what keeps single bad nights from producing a periodic signal. Returns the best surviving candidate as
    {'period','t0','duration_d','depth','depth_ppm','snr','n_in_transit','n_transits','flags'} or None.
    """
    from astropy.timeseries import BoxLeastSquares
    t, y, e = lc['t'], lc['flux'], lc['err']
    if band is not None:
        m = lc['filt'] == band
        t, y, e = t[m], y[m], e[m]
    if t.size < 3 * min_in_transit:
        return None
    base = float(t.max() - t.min())
    pmax = min(pmax, base / min_transits)
    if pmax <= pmin:
        return None
    model = BoxLeastSquares(t, y, dy=e)
    best = None
    for d in np.asarray(durations, float):
        lo = max(pmin, 2.0 * d)
        if pmax <= lo:
            continue
        res = model.power(period_grid(base, lo, pmax, dur=d), d, objective='snr')
        top = np.argsort(np.nan_to_num(res.depth_snr, nan=-np.inf))[::-1][:ncand]
        for i in top:
            P, t0 = float(res.period[i]), float(res.transit_time[i])
            if res.depth[i] <= 0:                      # only dimming events
                continue
            ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P
            inb = np.abs(ph) < 0.5 * d
            ep = np.unique(np.round((t[inb] - t0) / P))
            if inb.sum() < min_in_transit or ep.size < min_transits:
                continue
            if best is None or res.depth_snr[i] > best['snr']:
                best = {'period': P, 't0': t0, 'duration_d': float(d), 'depth': float(res.depth[i]),
                        'depth_ppm': float(res.depth[i] * 1e6), 'snr': float(res.depth_snr[i]),
                        'n_in_transit': int(inb.sum()), 'n_transits': int(ep.size)}
    if best is not None:
        best['flags'] = alias_flags(best['period'])
    return best


def alias_flags(P, tol=0.02):
    """Warn when the period sits on a sampling alias of the one-night cadence (P near n/k days, k = 1, 2, 3)."""
    for k in (1, 2, 3):
        x = P * k
        if round(x) >= 1 and abs(x - round(x)) < tol * k:
            return f'alias_{k}d'
    return ''


def search_position(ra, dec, radius_deg=2.0 / 3600, **kw):
    """Fetch and search a position in one call; returns {'lc': light curve, 'result': search dict or None}."""
    lc = fetch_lightcurve(ra, dec, radius_deg)
    return {'lc': lc, 'result': search(lc, **kw)}


def main(args):
    """occulta ztf RA DEC [--radius-arcsec R] [--pmin D] [--pmax D] [--band zg|zr|zi]: fetch the public IRSA ZTF light
    curve at a position (degrees, ICRS) and run the sparse-cadence BLS search. Prints a summary; needs the network."""
    import argparse, json
    ap = argparse.ArgumentParser(prog='occulta ztf', description=main.__doc__)
    ap.add_argument('ra', type=float); ap.add_argument('dec', type=float)
    ap.add_argument('--radius-arcsec', type=float, default=2.0)
    ap.add_argument('--pmin', type=float, default=0.2); ap.add_argument('--pmax', type=float, default=30.0)
    ap.add_argument('--band', choices=FILTERS)
    a = ap.parse_args(args)
    try:
        out = search_position(a.ra, a.dec, a.radius_arcsec / 3600, pmin=a.pmin, pmax=a.pmax, band=a.band)
    except ZtfError as e:
        raise SystemExit(f'ZTF: {e}')
    lc = out['lc']
    print(f"ZTF object ids {lc['oid']} (units {[UNIT_OF[f] for f in lc['oid']]}); {lc['n_good']} good of {lc['n_raw']} points, "
          f"BTJD {lc['t'].min():.1f} to {lc['t'].max():.1f}")
    print(json.dumps(out['result'], indent=1) if out['result'] else 'no candidate passes the cuts (too sparse or no periodic dip)')
