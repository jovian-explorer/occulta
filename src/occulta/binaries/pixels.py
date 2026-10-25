"""Difference-image source assignment for an eclipse in crowded TESS pixels.
For every sector, a 15x15 full-frame-image cutout (TESScut) is downloaded. Around each eclipse, the mean image in
eclipse (|dt| < 0.25 duration) is subtracted from the mean image just outside it (0.75 < |dt|/duration < 1.5, both
sides), which gives a difference image containing only the light lost in the eclipse. The out-of-eclipse image is
modelled as a sum of Gaussian point-spread functions at the positions of the Gaia DR3 sources (G < G_target + 4,
projected with the cutout WCS), with a common width, a common small pointing offset, and a flat background; the
difference image is then fitted with the same point-spread functions, and the fitted amplitude of each source,
divided by its out-of-eclipse amplitude, is the fraction of that star's light lost in the eclipse. The eclipsing
star is the source that alone reproduces the difference image; the errors come from the scatter between eclipses.
The flux-weighted centroid of the difference image is also reported relative to each source.
Usage: occulta binaries pixels <TIC> <P_d> <T0_btjd> <duration_h> <ra> <dec> [phase=0]
Writes binaries/eb/diffimg/<TIC>_s<sector>.png and prints a table per sector."""
import sys
from pathlib import Path

from .. import paths

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = paths.HOME
OUT = paths.BIN / 'eb' / 'diffimg'


def gaia_sources(ra, dec, radius_arcsec=150, gmax=None):
    from astroquery.gaia import Gaia
    q = (f"select source_id, ra, dec, phot_g_mean_mag as g, parallax, ruwe from gaiadr3.gaia_source where "
         f"1=CONTAINS(POINT({ra},{dec}),CIRCLE(ra,dec,{radius_arcsec / 3600})) order by phot_g_mean_mag")
    r = Gaia.launch_job(q).get_results()
    if gmax is not None:
        r = r[r['g'] < gmax]
    return r


def merge_close(src, r_arcsec):
    """Sources closer than r_arcsec cannot be separated in 20-arcsec pixels: merge them into one blend at the
    flux-weighted position, named after the brightest member."""
    from astropy.table import Table
    ra, de, g = np.asarray(src['ra'], float), np.asarray(src['dec'], float), np.asarray(src['g'], float)
    used = np.zeros(len(src), bool); rows = []
    for i in np.argsort(g):
        if used[i]:
            continue
        d = np.hypot((ra - ra[i]) * np.cos(np.radians(de[i])), de - de[i]) * 3600
        m = (d < r_arcsec) & ~used
        used |= m
        f = 10 ** (-0.4 * g[m])
        rows.append((int(src['source_id'][i]), float(np.sum(f * ra[m]) / f.sum()), float(np.sum(f * de[m]) / f.sum()),
                     float(-2.5 * np.log10(f.sum())), int(m.sum())))
    return Table(rows=rows, names=('source_id', 'ra', 'dec', 'g', 'n_members'))


def load(path):
    with fits.open(path) as h:
        d = h[1].data
        t, cube, q = np.asarray(d['TIME'], float), np.asarray(d['FLUX'], float), np.asarray(d['QUALITY'], int)
        sector = h[0].header.get('SECTOR')
        w = WCS(h[2].header)
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    return sector, t[ok], cube[ok], w


def psf_stack(xs, ys, shape, sig, dx=0.0, dy=0.0):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    return np.array([np.exp(-0.5 * ((xx - x - dx) ** 2 + (yy - y - dy) ** 2) / sig ** 2) for x, y in zip(xs, ys)])


def fit_out(img, xs, ys):
    """Out-of-eclipse model: common width, common offset, free amplitudes, flat background."""
    def res(p):
        sig, dx, dy, b = p[:4]
        return (np.tensordot(p[4:], psf_stack(xs, ys, img.shape, sig, dx, dy), 1) + b - img).ravel()
    b0 = np.median(img)
    a0 = [max(img[int(np.clip(round(y), 0, img.shape[0] - 1)), int(np.clip(round(x), 0, img.shape[1] - 1))] - b0, 1.0)
          for x, y in zip(xs, ys)]
    p0 = np.r_[0.8, 0.0, 0.0, b0, a0]
    lo = np.r_[0.4, -1.0, -1.0, -np.inf, np.zeros(len(xs))]
    hi = np.r_[1.6, 1.0, 1.0, np.inf, np.full(len(xs), np.inf)]
    r = least_squares(res, p0, bounds=(lo, hi))
    return r.x


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    tic, P, T0, dur_h, ra, dec = int(argv[0]), *map(float, argv[1:6])
    phase = float(argv[6]) if len(argv) > 6 else 0.0
    dur = dur_h / 24
    OUT.mkdir(parents=True, exist_ok=True)
    work = ROOT / 'data' / 'cache' / f'tesscut15_{tic}'
    work.mkdir(parents=True, exist_ok=True)
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    paths = sorted(work.glob('*.fits'))
    if not paths:
        dl = Tesscut.download_cutouts(coordinates=SkyCoord(ra, dec, unit='deg'), size=15, path=str(work))
        paths = [Path(p) for p in dl['Local Path']]
    src = gaia_sources(ra, dec)
    gt = src['g'][np.argmin((src['ra'] - ra) ** 2 + (src['dec'] - dec) ** 2)]
    src = src[src['g'] < gt + 4]
    src = merge_close(src, 8.0)
    summary = []
    for path in paths:
        try:
            sector, t, cube, w = load(path)
        except Exception as e:
            print('failed', path.name, repr(e)[:80]); continue
        xs, ys = w.all_world2pix(np.asarray(src['ra'], float), np.asarray(src['dec'], float), 0)
        inside = (xs > -2) & (xs < cube.shape[2] + 1) & (ys > -2) & (ys < cube.shape[1] + 1)
        s_in = src[inside]; xs, ys = xs[inside], ys[inside]
        n = np.round((t - T0 - phase * P) / P)
        dt = t - (T0 + (n + phase) * P)
        diffs, outs = [], []
        for k in np.unique(n):
            m = n == k
            ii = m & (np.abs(dt) < 0.25 * dur)
            oo = m & (np.abs(dt) > 0.75 * dur) & (np.abs(dt) < 1.5 * dur)
            if ii.sum() < 3 or (oo & (dt < 0)).sum() < 3 or (oo & (dt > 0)).sum() < 3:
                continue
            o = np.median(cube[oo], axis=0)
            diffs.append(o - np.median(cube[ii], axis=0)); outs.append(o)
        if not diffs:
            print(f'sector {sector}: no usable eclipse'); continue
        D = np.mean(diffs, axis=0); O = np.mean(outs, axis=0)
        p = fit_out(O, xs, ys)
        sig, dx, dy = p[:3]
        if sig > 1.55 or max(abs(dx), abs(dy)) > 0.95:
            print(f'sector {sector}: out-of-eclipse model failed (sigma {sig:.2f}, offset {dx:+.2f},{dy:+.2f}); skipped')
            continue
        A = psf_stack(xs, ys, O.shape, sig, dx, dy)
        F = p[4:]
        # single-source fits: which one star alone best reproduces D
        chi1, frac1, efrac1 = [], [], []
        for j in range(len(xs)):
            Xj = np.vstack([A[j].ravel(), np.ones(O.size)]).T
            c, *_ = np.linalg.lstsq(Xj, D.ravel(), rcond=None)
            chi1.append(np.sum((Xj @ c - D.ravel()) ** 2))
            frac1.append(c[0] / F[j] if F[j] > 0 else np.nan)
            ce = [np.linalg.lstsq(Xj, d.ravel(), rcond=None)[0][0] for d in diffs]
            efrac1.append(np.std(ce) / np.sqrt(len(ce)) / F[j] if len(ce) > 1 and F[j] > 0 else np.nan)
        chi1, frac1, efrac1 = np.array(chi1), np.array(frac1), np.array(efrac1)
        # flux-weighted centroid of the positive part of the difference image
        yy, xx = np.mgrid[:D.shape[0], :D.shape[1]]
        Dp = np.clip(D, 0, None)
        cx, cy = (Dp * xx).sum() / Dp.sum(), (Dp * yy).sum() / Dp.sum()
        scale = np.sqrt(np.abs(np.linalg.det(w.pixel_scale_matrix))) * 3600
        print(f'\nsector {sector}: {len(diffs)} eclipses, PSF sigma {sig:.2f} px, offset ({dx:+.2f}, {dy:+.2f}) px, '
              f'pixel {scale:.1f} arcsec')
        order = np.argsort(chi1)
        srt = np.sort(chi1)
        summary.append({'sector': int(sector), 'n_ecl': len(diffs), 'best': int(s_in['source_id'][order[0]]),
                        'best_g': float(s_in['g'][order[0]]),
                        'chi2_ratio_next': float(srt[1] / srt[0]) if len(srt) > 1 else None,
                        'depth_if_alone': {int(s_in['source_id'][j]): float(frac1[j]) for j in order[:4]},
                        'chi2_rel': {int(s_in['source_id'][j]): float(chi1[j] / chi1.min()) for j in order[:4]}})
        for j in order[:6]:
            sep = np.hypot(float(s_in['ra'][j] - ra) * np.cos(np.radians(dec)), float(s_in['dec'][j] - dec)) * 3600
            dcen = np.hypot(cx - xs[j] - dx, cy - ys[j] - dy) * scale
            print(f"  Gaia {s_in['source_id'][j]}{'+' + str(s_in['n_members'][j] - 1) if s_in['n_members'][j] > 1 else ''} G {s_in['g'][j]:5.2f} sep {sep:5.1f}\" : depth if alone {100*frac1[j]:7.3f} "
                  f"+- {100*efrac1[j]:.3f} %, single-source rel. chi2 {chi1[j]/chi1.min():6.2f}, diff-centroid offset {dcen:5.1f}\"")
        fig, ax = plt.subplots(1, 2, figsize=(9, 4.3))
        for a, img, ttl in ((ax[0], O, 'out of eclipse'), (ax[1], D, 'out minus in (light lost)')):
            a.imshow(img, origin='lower', cmap='Greys_r' if a is ax[0] else 'viridis')
            a.plot(xs + dx, ys + dy, 'r+', ms=8)
            a.plot(cx, cy, 'wx', ms=10, mew=2) if a is ax[1] else None
            a.set_title(f'TIC {tic} S{sector}: {ttl}', fontsize=9)
        j = order[0]
        ax[1].plot(xs[j] + dx, ys[j] + dy, 'o', mfc='none', mec='r', ms=14)
        fig.tight_layout(); fig.savefig(OUT / f'{tic}_s{sector}.png', dpi=100); plt.close(fig)


    # verdict: the source that is best in most sectors, and how decisively (median chi2 ratio to the runner-up)
    import json
    from collections import Counter
    v = {'tic': tic, 'sectors': summary}
    if summary:
        cnt = Counter(x['best'] for x in summary)
        best, nb = cnt.most_common(1)[0]
        tgt = int(src['source_id'][np.argmin((src['ra'] - ra) ** 2 + (src['dec'] - dec) ** 2)])
        rel_t = [x['chi2_rel'].get(tgt) for x in summary if x['best'] != tgt and x['chi2_rel'].get(tgt)]
        v.update({'best_source': best, 'n_best': nb, 'n_sectors': len(summary), 'target_source': tgt,
                  'on_target': best == tgt,
                  'median_chi2_ratio_next': float(np.median([x['chi2_ratio_next'] for x in summary if x['best'] == best
                                                             and x['chi2_ratio_next']] or [np.nan])),
                  'median_target_chi2_rel_when_not_best': float(np.median(rel_t)) if rel_t else None})
        print(f"\nVERDICT TIC {tic}: best source {best} in {nb}/{len(summary)} sectors "
              f"({'the target' if best == tgt else 'NOT the target'}); median chi2 ratio to next {v['median_chi2_ratio_next']:.2f}")
    (OUT / f'{tic}.json').write_text(json.dumps(v, indent=1))


if __name__ == '__main__':
    main()
