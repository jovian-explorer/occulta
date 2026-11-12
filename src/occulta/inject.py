"""Injection-recovery tests on real light curves.

Synthetic limb-darkened transits are multiplied into the flux of a real TESS light curve, which is then passed through
the same detrending, single-event search and BLS that the survey uses, with the survey thresholds. A transit counts
as recovered when the per-sector BLS returns its period (or half or twice it) with the promotion cuts, or when at least
two injected transits are found as single events above `event_snr`, which is what `survey.promote` needs to make a
'multi' candidate. Running this over a random sample of star-sectors gives the completeness of the pipeline as a
function of period, depth and noise.
"""
import numpy as np
from . import lightcurves as L, paths

CFG = paths.config()


def draw(rng, span):
    """One random transit: log-uniform period 0.5-30 d (capped at the sector span), log-uniform depth 100 ppm-3%,
    duration of a central transit across a Sun-like star scaled by sqrt(1-b^2)."""
    P = float(np.exp(rng.uniform(np.log(0.5), np.log(min(30.0, max(1.0, span))))))
    depth = float(np.exp(rng.uniform(np.log(1e-4), np.log(3e-2))))
    b = rng.uniform(0, 0.9)
    dur_h = float(13.0 * (P / 365.25) ** (1 / 3) * np.sqrt(1 - b ** 2))
    return P, depth, max(dur_h, 0.5)


def model(t, P, t0, depth, dur_h):
    """Trapezoid transit, ingress and egress each 15% of the duration."""
    d = dur_h / 24.0
    x = np.abs(((t - t0) / P + 0.5) % 1 - 0.5) * P
    tin = 0.15 * d
    m = np.clip((d / 2 - x) / tin, 0, 1)
    return 1.0 - depth * m


def search(t, f, bkg, cx, cy):
    """The detection half of lightcurves.process, returning events, BLS result and the inverted-flux SNR."""
    o = np.argsort(t); t, f, bkg, cx, cy = t[o], f[o], bkg[o], cx[o], cy[o]
    y, segs = L.detrend(t, f)
    segs_t = [(t[s[0]], t[s[-1]]) for s in segs if s.size]
    tg, dt, (yg, bg, xg, yg2) = L.to_grid(t, y, bkg, cx, cy)
    sigma = L.mad_std(yg)
    evs = L.single_events(tg, dt, yg, sigma)
    interior = np.ones(tg.size, bool)
    for a, b in segs_t:
        interior &= ~(((tg > a - 0.01) & (tg < a + 0.3)) | ((tg > b - 0.3) & (tg < b + 0.01)))
    inv = L.single_events(tg, dt, yg, sigma, invert=True, interior=interior)
    events = [L.characterise(e, tg, dt, yg, bg, xg, yg2, segs_t) for e in evs]
    bls = L.run_bls(t, y, t[-1] - t[0])
    k1 = max(1, int(round(1 / 24 / dt)))
    _, m_in, _, _, ok = L.box_scan(yg, k1)
    return events, bls, inv, float(L.mad_std(np.where(ok, m_in, np.nan)) * 1e6)


def recovered(P, t0, dur_h, epochs, events, bls, inv):
    """Apply the survey promotion cuts to what the search returned for one injected transit."""
    tol = dur_h / 48 + 0.1
    hit = [e for e in events if e['snr'] >= CFG['event_snr'] and e['snr'] >= CFG['inv_ratio'] * inv
           and e['flags'] in ('', 'edge') and np.min(np.abs(epochs - e['t0'])) < tol]
    nev = len({int(np.argmin(np.abs(epochs - e['t0']))) for e in hit})
    bls_ok = False
    if bls and bls['sde'] >= CFG['bls_sde'] and bls['snr'] >= CFG['bls_snr'] and bls['ntr'] >= 3 \
            and bls['oddeven_sig'] < 3:
        for h in (0.5, 1.0, 2.0):
            if abs(bls['period'] / (P * h) - 1) < 0.01:
                ph = ((bls['t0'] - t0) / (P * h) + 0.5) % 1 - 0.5
                bls_ok |= abs(ph * P * h) < tol
    return nev, bool(bls_ok)


def run(path, fmt='qlp', n=5, seed=0):
    """`n` independent injections into one light curve. Returns a list of result dicts (empty on unreadable files)."""
    rng = np.random.default_rng(seed)
    try:
        t, f, bkg, cx, cy, meta = L.read_lc(path, fmt)
    except Exception:
        return []
    if t.size < 300:
        return []
    span = float(t.max() - t.min())
    out = []
    if n == 0:      # null run: the untouched light curve, for the false-alarm side of any cut change
        events, bls, inv, n1h = search(t, f, bkg, cx, cy)
        return [{'P': 0.0, 't0': 0.0, 'depth_ppm': 0.0, 'dur_h': 0.0, 'ntr': bls['ntr'] if bls else 0,
                 'nev': sum(1 for e in events if e.get('flags', '') in ('', 'edge')), 'bls_ok': False, 'rec': False,
                 'noise1h_ppm': n1h, 'inv_snr': float(inv), 'cam': meta['cam'], 'ccd': meta['ccd'],
                 'bls_P': bls['period'] if bls else np.nan, 'bls_sde': bls['sde'] if bls else np.nan,
                 'bls_snr': bls['snr'] if bls else np.nan, 'bls_oe': bls['oddeven_sig'] if bls else np.nan,
                 'ev_max': max([e['snr'] for e in events], default=0.0)}]
    for _ in range(n):
        P, depth, dur_h = draw(rng, span)
        t0 = float(t.min() + rng.uniform(0, P))
        epochs = t0 + P * np.arange(0, int(span / P) + 2)
        near = np.array([np.any(np.abs(t - e) < dur_h / 48) for e in epochs])
        epochs = epochs[near]
        if epochs.size == 0:
            continue
        try:
            events, bls, inv, n1h = search(t, f * model(t, P, t0, depth, dur_h), bkg, cx, cy)
        except Exception:
            continue
        nev, bls_ok = recovered(P, t0, dur_h, epochs, events, bls, inv)
        out.append({'P': P, 't0': t0, 'depth_ppm': depth * 1e6, 'dur_h': dur_h, 'ntr': int(epochs.size),
                    'nev': nev, 'bls_ok': bls_ok, 'rec': bool(bls_ok or nev >= 2), 'noise1h_ppm': n1h,
                    'inv_snr': float(inv), 'cam': meta['cam'], 'ccd': meta['ccd'],
                    'bls_P': bls['period'] if bls else np.nan, 'bls_sde': bls['sde'] if bls else np.nan,
                    'bls_snr': bls['snr'] if bls else np.nan, 'bls_oe': bls['oddeven_sig'] if bls else np.nan,
                    'ev_max': max([e['snr'] for e in events], default=0.0)})
    return out
