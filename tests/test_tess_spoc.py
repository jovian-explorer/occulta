"""TESS-SPOC FFI light curves (HLSP tess-spoc): unit codes, target list, url, reader. Offline."""
import numpy as np
from astropy.io import fits
from occulta import sources, lightcurves as lc


def test_unit_code_round_trip_and_range():
    us = sources.units('TESS-SPOC')
    assert us[0] == 4104 and us[-1] == 4001
    assert not set(us) & set(sources.units('SPOC') + sources.units('QLP') + sources.units('K2') + sources.units('Kepler'))
    assert sources.unit_name(4014) == 'TESS S14 FFI SPOC'
    assert sources.unit_name(1014) == 'TESS S14 2-min' and sources.unit_name(3002) == 'Kepler file 2'
    assert sources.fmt_of_unit(4014) == 'spoc'


def test_url_matches_the_hlsp_layout():
    url = sources.url_for(4014, 7547522)[0]
    assert url.endswith('mast:HLSP/tess-spoc/s0014/target/0000/0000/0754/7522/'
                        'hlsp_tess-spoc_tess_phot_0000000007547522-s0014_tess_v1_lc.fits')


def test_target_list_parsing(monkeypatch):
    seen = []

    class R:
        status_code = 200
        text = '#TIC_ID,RA,DEC\n0000000007547522,272.19,47.84\n0000000007547527,272.31,47.83\n'
    monkeypatch.setattr(sources.requests, 'get', lambda u, timeout: seen.append(u) or R())
    assert sources.unit_targets(4014) == {7547522: None, 7547527: None}
    assert seen == ['https://archive.stsci.edu/hlsps/tess-spoc/target_lists/s0014.csv']
    R.status_code = 404
    assert sources.unit_targets(4014) is None


def test_reader_uses_pdcsap_and_quality_zero(tmp_path):
    t = np.linspace(1700, 1710, 500)
    f = np.full(500, 1000.0); f[100] = np.nan
    q = np.zeros(500, int); q[200] = 8
    cols = [fits.Column(name=n, format='D', array=a) for n, a in
            [('TIME', t), ('PDCSAP_FLUX', f), ('SAP_FLUX', f * 2), ('SAP_BKG', np.ones(500)),
             ('MOM_CENTR1', np.ones(500)), ('MOM_CENTR2', np.ones(500))]]
    cols.append(fits.Column(name='QUALITY', format='J', array=q))
    h0 = fits.PrimaryHDU(); h0.header['CAMERA'] = 2; h0.header['CCD'] = 3
    p = tmp_path / 'hlsp_tess-spoc_tess_phot_0000000000000007-s0014_tess_v1_lc.fits'
    fits.HDUList([h0, fits.BinTableHDU.from_columns(cols)]).writeto(p)
    tt, ff, _, _, _, meta = lc.read_lc(str(p), sources.fmt_of_unit(4014))
    assert tt.size == 498 and np.all(ff == 1000.0) and meta == {'cam': 2, 'ccd': 3}
