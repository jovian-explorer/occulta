"""Per-eclipse depths and times, sector by sector, for one eclipsing binary (cached light curve from eb_vet).
Each eclipse is fitted with a symmetric generalised-Gaussian dip on a local linear baseline,
    f(t) = a + b (t - tc) - D exp(-|(t - tc)/s|^g),
which gives its mid-time tc, depth D, and width; the errors come from the covariance matrix scaled by the
reduced chi-square. Within each sector, the secondary phase is measured relative to the primaries of the same
sector, so that it is independent of the long-term period and of any timestamp offset common to the sector.
Usage: occulta binaries eclipses <TIC> <P_d> <T0_primary_btjd> <phase2> [width_h]
Writes binaries/eb/<TIC>_eclipses.csv and <TIC>_by_sector.png."""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from . import analysis as eb_vet
UNIT_NAME = lambda u: f'SPOC S{u - 1000}' if u > 1000 else f'QLP S{u}'


def model(t, a, b, tc, D, s, g):
    return a + b * (t - tc) - D * np.exp(-np.abs((t - tc) / s) ** g)


def fit_one(t, y, tc0, w):
    # re-centre on the observed minimum within +-0.6 eclipse widths of the prediction (the long-term ephemeris may
    # be off by hours when the orbit precesses), using a 5-point running median to suppress single-point noise
    s = np.abs(t - tc0) < 0.6 * w
    if s.sum() >= 5:
        ys = np.convolve(y[s], np.ones(5) / 5, mode='same')
        tc0 = t[s][np.argmin(ys)]
    m = np.abs(t - tc0) < 1.5 * w
    if m.sum() < 15 or np.sum(np.abs(t[m] - tc0) < 0.3 * w) < 3:
        return None
    tt, yy = t[m], y[m]
    p0 = [np.median(yy[np.abs(tt - tc0) > 0.7 * w]) if np.any(np.abs(tt - tc0) > 0.7 * w) else 0.0, 0.0, tc0,
          max(1e-4, -yy.min()), 0.3 * w, 2.0]
    try:
        p, cov = curve_fit(model, tt, yy, p0=p0, bounds=([-1, -1, tc0 - 0.3 * w, 0, 0.02 * w, 0.8], [1, 1, tc0 + 0.3 * w, 1, 2 * w, 8]),
                           maxfev=20000)
    except Exception:
        return None
    r = yy - model(tt, *p)
    chi = np.sum(r ** 2) / max(1, tt.size - 6)
    e = np.sqrt(np.diag(cov))
    # covered both sides of the centre?
    if np.sum(tt < p[2] - 0.3 * p[4]) < 3 or np.sum(tt > p[2] + 0.3 * p[4]) < 3:
        return None
    return p[2], e[2], p[3], e[3], p[4], np.sqrt(chi)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    tic, P, T0, ph2 = int(argv[0]), float(argv[1]), float(argv[2]), float(argv[3])
    w = float(argv[4]) / 24 if len(argv) > 4 else 0.25
    T, F, U, *_ = eb_vet.load(tic)
    rows = []
    for u in np.unique(U):
        m = U == u
        t = T[m]; y = eb_vet.normalise(t, F[m], w)
        ok = np.isfinite(y); t, y = t[ok], y[ok]
        for kind, phase in (('P', 0.0), ('S', ph2)):
            n0 = np.ceil((t[0] - T0 - phase * P) / P); n1 = np.floor((t[-1] - T0 - phase * P) / P)
            for n in np.arange(n0, n1 + 1):
                r = fit_one(t, y, T0 + (n + phase) * P, w)
                if r:
                    rows.append((int(u), kind, int(n), *r))
    rows = np.array(rows, dtype=object)
    out = Path(eb_vet.OUT) / f'{tic}_eclipses.csv'
    with open(out, 'w') as fh:
        fh.write('unit,type,cycle,t_btjd,t_err_d,depth,depth_err,sigma_d,rms\n')
        for r in rows:
            fh.write(','.join(str(x) for x in r) + '\n')
    # per-sector summary: mean depths and secondary phase from same-sector primaries
    print(f'{"unit":>10} {"mid BTJD":>9} {"D_pri %":>14} {"D_sec %":>14} {"ratio":>6} {"phase_2":>16}')
    summ = []
    for u in sorted(set(r[0] for r in rows)):
        pr = [r for r in rows if r[0] == u and r[1] == 'P']; se = [r for r in rows if r[0] == u and r[1] == 'S']
        if not pr:
            continue
        dp = np.array([r[5] for r in pr]); ds = np.array([r[5] for r in se]) if se else np.array([np.nan])
        # local primary ephemeris in this sector (fixed P), then phase of each secondary
        tp = np.array([r[3] - r[2] * P for r in pr]); t0u = np.median(tp)
        phs = np.array([((r[3] - t0u) / P) % 1 for r in se]) if se else np.array([np.nan])
        mid = np.median([r[3] for r in pr])
        summ.append((u, mid, dp.mean(), dp.std() / np.sqrt(len(dp)), ds.mean(), ds.std() / np.sqrt(max(1, len(ds))),
                     np.nanmean(phs), np.nanstd(phs) / np.sqrt(max(1, len(phs)))))
        s = summ[-1]
        print(f'{UNIT_NAME(u):>10} {mid:9.1f} {100*s[2]:6.2f}±{100*s[3]:4.2f} {100*s[4]:6.2f}±{100*s[5]:4.2f} '
              f'{s[4]/s[2]:6.3f} {s[6]:.4f}±{s[7]:.4f}')
    summ = np.array(summ)
    fig, ax = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for k, (col, lab) in enumerate(((2, 'primary depth (%)'), (4, 'secondary depth (%)'))):
        ax[k].errorbar(summ[:, 1], 100 * summ[:, col], 100 * summ[:, col + 1], fmt='o', ms=6, color='0.1')
        ax[k].set_ylabel(lab)
    ax[2].errorbar(summ[:, 1], summ[:, 6], summ[:, 7], fmt='o', ms=6, color='0.1'); ax[2].set_ylabel('secondary phase')
    ax[2].set_xlabel('BTJD'); fig.suptitle(f'TIC {tic}: eclipses by sector (P = {P:.6f} d)')
    fig.tight_layout(); fig.savefig(Path(eb_vet.OUT) / f'{tic}_by_sector.png', dpi=90)


if __name__ == '__main__':
    main()
