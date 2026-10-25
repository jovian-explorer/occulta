"""Uniform aperture photometry from TESS full-frame-image cutouts (TESScut), to compare eclipse depths across
sectors without the sector-to-sector differences of pipeline apertures and background estimates.
For every sector, an 11x11 cutout is downloaded; the background is the median of the pixels outside a 5x5 box
around the target; the flux is the sum over a fixed square aperture of side n_ap pixels centred on the target
pixel. The same aperture size is used in every sector, so that the dilution by neighbours is the same in all.
Usage: occulta binaries photometry <TIC> <ra> <dec> [n_ap=3] -> binaries/eb/ffi/<TIC>_s<sector>.npz"""
import sys
from pathlib import Path

from .. import paths

import numpy as np
from astropy.io import fits

ROOT = paths.HOME
OUT = paths.BIN / 'eb' / 'ffi'


def phot(path, n_ap):
    with fits.open(path) as h:
        d = h[1].data
        t, cube, q = np.asarray(d['TIME'], float), np.asarray(d['FLUX'], float), np.asarray(d['QUALITY'], int)
        sector = h[0].header.get('SECTOR')
    ny, nx = cube.shape[1:]
    cy, cx = ny // 2, nx // 2
    yy, xx = np.mgrid[:ny, :nx]
    bmask = (np.abs(yy - cy) > 2) | (np.abs(xx - cx) > 2)
    amask = (np.abs(yy - cy) <= n_ap // 2) & (np.abs(xx - cx) <= n_ap // 2)
    bkg = np.nanmedian(cube[:, bmask], axis=1)
    f = np.nansum(cube[:, amask], axis=1) - amask.sum() * bkg
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    return sector, t[ok], f[ok]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    tic, ra, dec = int(argv[0]), float(argv[1]), float(argv[2])
    n_ap = int(argv[3]) if len(argv) > 3 else 3
    OUT.mkdir(parents=True, exist_ok=True)
    work = ROOT / 'data' / 'cache' / f'tesscut_{tic}'
    work.mkdir(parents=True, exist_ok=True)
    c = SkyCoord(ra, dec, unit='deg')
    for p in Tesscut.download_cutouts(coordinates=c, size=11, path=str(work)):
        try:
            s, t, f = phot(p['Local Path'], n_ap)
        except Exception as e:
            print('failed', p['Local Path'], repr(e)[:80]); continue
        np.savez_compressed(OUT / f'{tic}_s{s}_ap{n_ap}.npz', t=t, f=f)
        print('sector', s, t.size, 'points')
        Path(p['Local Path']).unlink(missing_ok=True)


if __name__ == '__main__':
    main()
