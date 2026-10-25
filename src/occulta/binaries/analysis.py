"""Eclipsing-binary vetting of the uncatalogued EB-like rejects (binaries/new_variables_isolated.csv).

For each star: (1) every sector is fetched again (cached as a compact light curve in binaries/eb/lc/);
(2) each sector is normalised by a 2-day running median computed with every dip masked (dips found iteratively,
so an eccentric secondary at any phase is protected); (3) the mid-time of every covered primary eclipse is
measured by fitting an empirical template, a linear ephemeris is fitted to those times, and the O-C residuals
are kept (eclipse timing variations from a third body or apsidal motion show up there); (4) odd and even
eclipses are compared inside each sector, so that sector-to-sector dilution changes cannot mimic a difference,
and the period is doubled when they differ consistently; (5) the secondary eclipse is searched over all
phases, and its own mid-times are measured; (6) the in-eclipse centroid shift is measured per sector.
QLP Sectors 99-104, whose timestamps await correction,
are excluded from all eclipse times. Every sample is drawn in the diagnostic figure.

Usage: occulta binaries analyse [--workers 6] | [--one TIC]  ->  binaries/eb/<TIC>.json/.png
"""
import json, shutil, sys, traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .. import paths

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, binary_dilation

ROOT = paths.HOME
from .. import lightcurves as lcsearch, planets as vet  # noqa: E402

OUT = paths.BIN / 'eb'
WORK = ROOT / 'data' / 'cache' / 'ebvet'
PIX = 21.0  # TESS pixel, arcsec (Ricker et al. 2015)
KEPLER_PIX = 3.98  # Kepler and K2 pixel, arcsec (Koch et al. 2010); units >= 2000 are K2 campaigns and Kepler quarters


# ------------------------------------------------------------------ data
def load(tic):
    """All sectors of one star as one table (t, f, unit, cx, cy), raw flux; cached in eb/lc/<tic>.npz."""
    cache = OUT / 'lc' / f'{tic}.npz'
    if cache.exists():
        z = np.load(cache)
        return z['t'], z['f'], z['u'], z['cx'], z['cy']
    wd = WORK / str(tic)
    wd.mkdir(parents=True, exist_ok=True)
    parts = []
    try:
        for author, unit, path in vet.fetch_all(tic, wd):
            try:
                t, f, _, cx, cy, _ = lcsearch.read_lc(path, vet.AUTHORS[author])
            except Exception:
                continue
            if t.size < 100:
                continue
            o = np.argsort(t)
            u = unit + (1000 if author in ('SPOC', 'TESS-SPOC') else 0)
            parts.append((t[o], f[o], np.full(t.size, u), cx[o], cy[o]))
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    if not parts:
        return None
    t, f, u, cx, cy = (np.concatenate([p[i] for p in parts]) for i in range(5))
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, t=t, f=f.astype(np.float32), u=u.astype(np.int16), cx=cx.astype(np.float32),
                        cy=cy.astype(np.float32))
    return t, f, u, cx, cy


def normalise(t, f, w):
    """Relative flux with slow trends removed; dips (eclipses) are found iteratively and masked from the trend."""
    dt = np.median(np.diff(t))
    win = max(5, int(2.0 / dt) | 1)
    grow = max(1, int(0.75 * w / dt))
    edges = np.where(np.diff(t) > 0.5)[0] + 1
    segs = [s for s in np.split(np.arange(t.size), edges) if s.size > 20]
    y = np.full(t.size, np.nan)
    mask = np.zeros(t.size, bool)
    for _ in range(3):
        for s in segs:
            g = f[s].astype(float).copy()
            good = ~mask[s]
            if good.sum() < 20:
                continue
            g[~good] = np.interp(t[s][~good], t[s][good], g[good])
            trend = median_filter(g, size=min(win, s.size - (1 - s.size % 2)), mode='nearest')
            y[s] = f[s] / trend - 1.0
        rms = lcsearch.mad_std(y[~mask & np.isfinite(y)])
        mask = binary_dilation(np.nan_to_num(y) < -5 * rms, iterations=grow)
    return y


# ------------------------------------------------------------------ eclipse timing
def template(x, y, half):
    """Binned median eclipse profile on a grid of x (days from centre)."""
    bins = np.linspace(-half, half, 61)
    c = 0.5 * (bins[1:] + bins[:-1])
    idx = np.digitize(x, bins) - 1
    prof = np.array([np.median(y[idx == i]) if np.sum(idx == i) >= 3 else np.nan for i in range(c.size)])
    ok = np.isfinite(prof)
    return c[ok], prof[ok]


def eclipse_times(T, Y, P, t0, w, rms, phase=0.0):
    """Mid-time and error of every well-covered eclipse at the given phase, from a template fit."""
    half = 1.5 * w
    x = (((T - t0 - phase * P) / P + 0.5) % 1 - 0.5) * P
    near = np.abs(x) < half
    if near.sum() < 30:
        return []
    tc, tp = template(x[near], Y[near], half)
    if tc.size < 10 or tp.min() > -3 * rms:
        return []
    # the template inherits any offset of the trial t0 (e.g. a survey t0 from QLP Sectors 99-104); centre it on its
    # axis of mirror symmetry so that the times refer to mid-eclipse
    grid = np.linspace(-0.5 * w, 0.5 * w, 201)
    xs = np.linspace(-0.6 * w, 0.6 * w, 61)
    asym = [np.nansum((np.interp(s + xs, tc, tp, left=np.nan, right=np.nan)
                       - np.interp(s - xs, tc, tp, left=np.nan, right=np.nan)) ** 2) for s in grid]
    s0 = grid[int(np.argmin(asym))]
    tc = tc - s0
    n = np.round((T - t0 - phase * P) / P)
    out = []
    shifts = np.linspace(-0.4 * w, 0.4 * w, 161)
    for k in np.unique(n[near]):
        c = t0 + phase * P + k * P + s0
        m = np.abs(T - c) < half
        inn = np.abs(T - c) < 0.5 * w
        if m.sum() < 15 or inn.sum() < max(5, 0.6 * w / np.median(np.diff(T[m]))):
            continue
        chi = np.array([np.sum((Y[m] - np.interp(T[m] - c - s, tc, tp, left=0, right=0)) ** 2) for s in shifts]) / rms ** 2
        i = int(np.argmin(chi))
        if i in (0, shifts.size - 1):
            continue
        a, b, _ = np.polyfit(shifts[i - 1:i + 2], chi[i - 1:i + 2], 2)
        if a <= 0:
            continue
        red = max(1.0, chi[i] / max(1, m.sum() - 1))
        out.append((int(k), c - b / (2 * a), np.sqrt(red / a)))
    return out


def ephemeris(times):
    """Weighted linear ephemeris T = T0 + E P from [(E, t, err)]; returns T0, P, their errors, O-C, chi2_red."""
    E = np.array([x[0] for x in times], float); t = np.array([x[1] for x in times]); e = np.array([x[2] for x in times])
    W = 1 / e ** 2
    A = np.vstack([np.ones_like(E), E]).T
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    p = cov @ (A.T @ (W * t))
    oc = t - A @ p
    chi = float(np.sum(W * oc ** 2) / max(1, E.size - 2))
    s = np.sqrt(max(1.0, chi))
    return p[0], p[1], np.sqrt(cov[0, 0]) * s, np.sqrt(cov[1, 1]) * s, oc, chi


def eclipse_width(T, Y, P, t0, w, level=0.1):
    """Full width (days) of the folded eclipse at `level` times its depth, from a finely binned median profile with
    linear interpolation of the two crossings (level=0.1: close to first-to-fourth contact; level=0.5: FWHM)."""
    x = (((T - t0) / P + 0.5) % 1 - 0.5) * P
    half = min(3 * w, 0.25 * P)
    near = np.abs(x) < half
    if near.sum() < 30:
        return w
    bins = np.linspace(-half, half, 241)
    c = 0.5 * (bins[1:] + bins[:-1])
    idx = np.digitize(x[near], bins) - 1
    yy = Y[near]
    prof = np.array([np.median(yy[idx == i]) if np.sum(idx == i) >= 3 else np.nan for i in range(c.size)])
    ok = np.isfinite(prof)
    c, prof = c[ok], prof[ok]
    if c.size < 20 or prof.min() >= 0:
        return w
    # smooth lightly (3-bin running median) so that single noisy bins do not end the walk early
    prof = median_filter(prof, size=3, mode='nearest')
    i0 = int(np.argmin(prof)); lev = level * prof[i0]
    lo = i0
    while lo > 0 and prof[lo - 1] < lev:
        lo -= 1
    hi = i0
    while hi < c.size - 1 and prof[hi + 1] < lev:
        hi += 1
    def cross(i, j):   # interpolate between bin i (inside) and j (outside)
        if j < 0 or j >= c.size or prof[j] == prof[i]:
            return c[i]
        return c[i] + (c[j] - c[i]) * (lev - prof[i]) / (prof[j] - prof[i])
    width = cross(hi, hi + 1) - cross(lo, lo - 1)
    return float(np.clip(width, 0.3 * w, 3 * w))


def dip(ph, y, c, half):
    """Mean flux deficit inside |phase - c| < half and its error."""
    d = np.abs(((ph - c + 0.5) % 1) - 0.5)
    m = d < half
    if m.sum() < 5:
        return np.nan, np.nan, 0
    return -np.mean(y[m]), np.std(y[m]) / np.sqrt(m.sum()), int(m.sum())


# ------------------------------------------------------------------ analysis
def analyse(tic, P, t0, t14_h, gaia):
    d = load(tic)
    if d is None:
        return {'tic': tic, 'status': 'no data'}
    T, F, U, CX, CY = d
    P_survey = P
    w = t14_h / 24
    Y = np.full(T.size, np.nan)
    for u in np.unique(U):
        m = U == u
        Y[m] = normalise(T[m], F[m], w)
    ok = np.isfinite(Y)
    T, Y, U, CX, CY = T[ok], Y[ok], U[ok].astype(int), CX[ok], CY[ok]
    o = np.argsort(T); T, Y, U, CX, CY = T[o], Y[o], U[o], CX[o], CY[o]

    # QLP Sectors 99-104 carry timestamps awaiting correction (offsets of minutes): they are used for depths and
    # shapes, never for eclipse times
    tm = ~((U >= 99) & (U <= 104))
    Tt, Yt = T[tm], Y[tm]

    # 1. measure the true eclipse width from the fold, then refine the period from individual eclipse times
    rms = lcsearch.mad_std(Y)
    w = eclipse_width(T, Y, P, t0, w)
    for _ in range(2):
        tt = eclipse_times(Tt, Yt, P, t0, w, rms)
        if len(tt) >= 3 and len({x[0] for x in tt}) >= 3:
            t0, P, *_ = ephemeris(tt)

    # 2. odd/even inside each sector
    z = []
    for u in np.unique(U):
        m = U == u
        n = np.round((T[m] - t0) / P)
        inn = np.abs(T[m] - (t0 + n * P)) < 0.25 * w
        dd = {0: [], 1: []}
        # only eclipses whose core is fully covered: one falling partly into a data gap or onto a segment edge
        # would look shallow and fake an odd/even difference
        need = max(3, 0.8 * 0.5 * w / np.median(np.diff(T[m])))
        for k in np.unique(n[inn]):
            mm = inn & (n == k)
            if mm.sum() >= need:
                dd[int(k) % 2].append(-np.mean(Y[m][mm]))
        if len(dd[0]) >= 2 and len(dd[1]) >= 2:     # one eclipse of each parity cannot show a real difference
            a, b = np.mean(dd[0]), np.mean(dd[1])
            sc = lcsearch.mad_std(Y[m][~(np.abs(T[m] - (t0 + n * P)) < w)]) / np.sqrt(
                max(3, int(0.5 * w / np.median(np.diff(T[m])))))
            # the depths of individual eclipses scatter more than the photometric noise allows (spots, trends):
            # use the scatter between eclipses of the same parity when it is larger
            if len(dd[0]) >= 2 and len(dd[1]) >= 2:
                sc = max(sc, np.sqrt(0.5 * (np.var(dd[0], ddof=1) + np.var(dd[1], ddof=1))))
            z.append(((a - b) / np.hypot(sc / np.sqrt(len(dd[0])), sc / np.sqrt(len(dd[1]))), abs(a - b) / max(a, b)))
    oe_z = np.array([x[0] for x in z]); oe_rel = np.array([x[1] for x in z])
    # alternating depths in (nearly) every sector: one discordant sector (e.g. a noisy one) must not veto the doubling
    sgn = np.sign(np.median(oe_z)) if z else 0
    double = bool(z) and np.median(np.abs(oe_z)) > 5 and np.median(oe_rel) > 0.05 and np.mean(np.sign(oe_z) == sgn) >= 0.75
    Pt = 2 * P if double else P
    if double and sgn < 0:
        t0 = t0 + P

    # 3. ephemeris of the primary at the adopted period
    tt = eclipse_times(Tt, Yt, Pt, t0, w, rms)
    eph = None
    if len(tt) >= 3 and len({x[0] for x in tt}) >= 3:
        t0, Pt, et0, eP, oc, chi = ephemeris(tt)
        eph = {'T0_btjd': t0, 'T0_err_d': et0, 'P_d': Pt, 'P_err_d': eP, 'n_times': len(tt),
               'oc_rms_min': float(np.std(oc) * 1440), 'chi2_red': chi,
               'times': [(int(E), float(t), float(e)) for (E, t, e) in tt]}
        if len(tt) >= 6:
            E = np.array([x[0] for x in tt], float); W = 1 / np.array([x[2] for x in tt]) ** 2
            eph['oc_quadratic_coef_d'] = float(np.polyfit(E, oc, 2, w=np.sqrt(W))[0])
    ph = ((T - t0) / Pt) % 1
    hw = 0.5 * w / Pt
    d1, e1, _ = dip(ph, Y, 0.0, 0.25 * hw)

    # 4. secondary eclipse over all phases, then its own times
    best = (np.nan, np.nan, 0.5)
    for c in np.arange(0.0, 1.0, max(hw / 4, 0.001)):
        if abs(((c + 0.5) % 1) - 0.5) < 2.2 * hw:
            continue
        dd, ee, nn = dip(ph, Y, c, 0.5 * hw)
        # a secondary eclipse is a local minimum: deeper than the flux one half-width to either side (this rejects
        # the slope of the primary's wing or of a detrending ramp next to the excluded zone)
        fl = [dip(ph, Y, c + s_ * hw, 0.5 * hw)[0] for s_ in (-1, 1)]
        if not all(np.isfinite(f_) and dd > f_ for f_ in fl):
            continue
        if nn >= 5 and np.isfinite(dd) and ee > 0 and (not np.isfinite(best[0]) or dd / ee > best[0] / best[1]):
            best = (dd, ee, c)
    d2, e2, ph2 = best
    sec = None
    if np.isfinite(d2) and d2 / e2 > 8:
        ts = eclipse_times(Tt, Yt, Pt, t0, w, rms, phase=ph2)
        if len(ts) >= 3 and len({x[0] for x in ts}) >= 3:
            s0, sP, es0, esP, soc, schi = ephemeris(ts)
            sec = {'phase': float(((s0 - t0) / Pt) % 1), 'phase_err': float(es0 / Pt), 'P_d': sP, 'P_err_d': esP,
                   'width_h': 24 * eclipse_width(T, Y, Pt, t0 + ph2 * Pt, w),
                   'fwhm_h': 24 * eclipse_width(T, Y, Pt, t0 + ph2 * Pt, w, 0.5),
                   'n_times': len(ts), 'oc_rms_min': float(np.std(soc) * 1440), 'chi2_red': schi}
    d05, e05, _ = dip(ph, Y, 0.5, 0.5 * hw)

    # 5. recurrence and centroids
    n = np.round((T - t0) / Pt)
    seen = cov = 0
    for k in np.unique(n):
        m = (n == k) & (np.abs(T - (t0 + k * Pt)) < 0.25 * w)
        if m.sum() >= 3:
            cov += 1
            seen += int(-np.mean(Y[m]) > 0.5 * d1)
    shifts = []
    for u in np.unique(U):
        m = (U == u) & np.isfinite(CX) & np.isfinite(CY)
        p = ph[m]
        inn = np.abs(((p + 0.5) % 1) - 0.5) < 0.25 * hw
        out = np.abs(((p + 0.5) % 1) - 0.5) > 2 * hw
        if inn.sum() >= 5 and out.sum() > 50:
            cx, cy = CX[m], CY[m]
            dx = np.median(cx[inn]) - np.median(cx[out]); dy = np.median(cy[inn]) - np.median(cy[out])
            ex = lcsearch.mad_std(cx[out]) / np.sqrt(inn.sum()); ey = lcsearch.mad_std(cy[out]) / np.sqrt(inn.sum())
            shifts.append((int(u), np.hypot(dx, dy) * (KEPLER_PIX if u >= 2000 else PIX), np.hypot(dx / ex, dy / ey) if ex > 0 and ey > 0 else np.nan))
    oot = lcsearch.mad_std(Y[np.abs(((ph + 0.5) % 1) - 0.5) > 3 * hw])

    res = {'tic': tic, 'period_survey': P_survey, 'period': Pt, 'doubled': bool(double), 't0_btjd': float(t0),
           'ephemeris': eph, 'primary_depth': float(d1), 'primary_err': float(e1),
           'secondary_depth': float(d2), 'secondary_err': float(e2), 'secondary_phase_scan': float(ph2), 'secondary': sec,
           'depth_at_half': float(d05), 'depth_at_half_err': float(e05), 'width_h': w * 24,
           'fwhm_h': 24 * eclipse_width(T, Y, Pt, t0, w, 0.5), 'width_survey_h': t14_h,
           'eclipses_covered': cov, 'eclipses_seen': seen, 'units': sorted(int(u) for u in np.unique(U)), 'oot_rms': float(oot),
           'odd_even_z_per_sector': [round(float(x), 1) for x in oe_z],
           'centroid': [{'unit': u, 'shift_arcsec': round(float(a), 2), 'sigma': round(float(b), 1)} for u, a, b in shifts],
           'gaia': gaia, 'status': 'ok'}

    # figure: every sample, opaque large dark points
    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.8])
    a = fig.add_subplot(gs[0, :])
    big = int(np.bincount(U - U.min()).argmax() + U.min())
    k = U == big
    a.plot(T[k], Y[k], 'o', ms=2.5, color='0.1')
    a.set_xlabel('BTJD'); a.set_ylabel('relative flux'); a.set_title(f'TIC {tic}, unit {big}')
    a = fig.add_subplot(gs[1, 0])
    a.plot(ph, Y, 'o', ms=2, color='0.1'); a.set_xlabel(f'phase (P = {Pt:.6f} d)'); a.set_ylabel('relative flux')
    for i, (c, lab) in enumerate(((0.0, 'primary'), (ph2, f'secondary (phase {ph2:.3f})'))):
        a = fig.add_subplot(gs[1, 1 + i])
        x = (((ph - c + 0.5) % 1) - 0.5) * Pt * 24
        m = np.abs(x) < 2.5 * w * 24
        a.plot(x[m], Y[m], 'o', ms=2.5, color='0.1')
        a.set_xlabel('hours from centre'); a.set_title(lab)
    a = fig.add_subplot(gs[2, :])
    if eph:
        E = np.array([x[0] for x in eph['times']]); tt_ = np.array([x[1] for x in eph['times']])
        ee = np.array([x[2] for x in eph['times']])
        a.errorbar(E, (tt_ - (t0 + E * Pt)) * 1440, ee * 1440, fmt='o', ms=4, color='0.1')
        a.axhline(0, color='0.5', lw=0.8)
        a.set_xlabel('primary eclipse number'); a.set_ylabel('O - C (min)')
        a.set_title(f'primary timing: rms {eph["oc_rms_min"]:.2f} min, reduced chi2 {eph["chi2_red"]:.1f}')
    sec_txt = f'{d2*100:.2f}% ({d2/e2:.0f} sigma) at phase {ph2:.3f}' if np.isfinite(d2) else 'none'
    fig.suptitle(f'TIC {tic}: P {Pt:.6f} d{" (doubled)" if double else ""}, primary {d1*100:.2f}%, secondary {sec_txt}, '
                 f'seen {seen}/{cov}')
    fig.tight_layout()
    fig.savefig(OUT / f'{tic}.png', dpi=80)
    plt.close(fig)
    (OUT / f'{tic}.json').write_text(json.dumps(res, indent=1, default=float))
    return res


def one(row):
    tic = int(row['tic'])
    if (OUT / f'{tic}.json').exists():
        return
    s = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())
    try:
        analyse(tic, s['period_d'], s['t0_bjd'] - 2457000.0, s['t14_h'], s.get('gaia', {}))
        print(tic, 'ok', flush=True)
    except Exception:
        (OUT / f'{tic}.json').write_text(json.dumps({'tic': tic, 'status': 'error', 'trace': traceback.format_exc()[-800:]}))
        print(tic, 'error', flush=True)


def main(argv=None):
    """Analyse every uncatalogued stellar signal (all of them, --workers N in parallel, or --one TIC)."""
    argv = sys.argv[1:] if argv is None else list(argv)
    from astropy.table import Table
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(Table.read(paths.BIN / 'new_variables_isolated.csv'))
    workers = 6
    if len(argv) > 1 and argv[0] == '--one':
        rows = [r for r in rows if int(r['tic']) == int(argv[1])]
    elif len(argv) > 1 and argv[0] == '--workers':
        workers = int(argv[1])
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, rows))


if __name__ == '__main__':
    main()
