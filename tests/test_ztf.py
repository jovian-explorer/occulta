"""Offline tests of occulta.ztf: format, filtering, table parsing and sparse-cadence recovery."""
import os
import numpy as np
import pytest

from occulta import ztf

CSV = """oid,expid,hjd,mjd,mag,magerr,catflags,filtercode,ra,dec
1,10,2458200.5,58200.0,15.00,0.01,0,zr,10.0,20.0
1,11,2458201.5,58201.0,15.10,0.01,0,zr,10.0,20.0
1,12,2458202.5,58202.0,15.00,0.01,32768,zr,10.0,20.0
1,13,2458203.5,58203.0,15.00,0.30,0,zr,10.0,20.0
1,14,2458204.5,58204.0,15.00,0.01,0,zg,10.0,20.0
2,15,2458205.5,58205.0,15.00,0.01,0,zr,10.0,20.0
"""

IPAC = """\\fixlen = T
| oid   | hjd         | mag   | magerr | catflags | filtercode |
| long  | double      | float | float  | int      | char       |
  1       2458200.5     15.00   0.01     0          zr
  1       2458201.5     15.10   0.01     32768      zr
"""


def synth(seed=1, n_nights=1000, period=3.7, depth=0.03, dur=0.12, sig=0.004):
    """About 3 years, one point per night at a random time of the night, with a box transit injected."""
    rng = np.random.default_rng(seed)
    nights = np.sort(rng.choice(np.arange(int(n_nights * 1.1)), n_nights, replace=False))
    t = 2800.0 + nights + rng.uniform(0, 0.3, n_nights)
    flux = 1.0 + rng.normal(0, sig, t.size)
    ph = ((t - 2900.3 + 0.5 * period) % period) - 0.5 * period
    flux[np.abs(ph) < dur / 2] -= depth
    return {'t': t, 'flux': flux, 'err': np.full(t.size, sig), 'filt': np.full(t.size, 'zr')}


def test_parse_csv_and_ipac():
    c = ztf.parse_table(CSV)
    assert len(c['hjd']) == 6 and c['filtercode'][4] == 'zg' and c['catflags'][2] == 32768
    i = ztf.parse_table(IPAC)
    assert list(i['catflags']) == [0, 32768] and i['mag'][1] == pytest.approx(15.1)


def test_parse_errors():
    with pytest.raises(ztf.ZtfError):
        ztf.parse_table('')
    with pytest.raises(ztf.ZtfError):
        ztf.parse_table('a,b\n1,2\n')
    with pytest.raises(ztf.ZtfError):
        ztf.parse_table('ERROR: bad request')


def test_flag_filtering():
    g = ztf.good_mask(ztf.parse_table(CSV))
    assert list(g) == [True, True, False, False, True, True]   # catflags 32768 and magerr 0.30 are dropped


def test_mag_to_flux():
    f, e = ztf.mag_to_flux([15.0, 15.0, 16.0], [0.01, 0.01, 0.01])
    assert np.median(f) == pytest.approx(1.0)
    assert f[2] == pytest.approx(10 ** -0.4, rel=1e-9)         # one magnitude fainter than the median
    assert e[0] == pytest.approx(0.4 * np.log(10) * 0.01)


def test_build_lightcurve_separates_filters_and_btjd():
    rows = ['oid,hjd,mag,magerr,catflags,filtercode']
    for i in range(40):
        rows.append(f'1,{2458200.5 + i},15.0,0.01,0,zr')
        rows.append(f'3,{2458200.5 + i},17.0,0.01,0,zg')
    rows.append('2,2458300.5,15.0,0.01,0,zr')                   # second object, fewer points: not used
    lc = ztf.build_lightcurve(ztf.parse_table('\n'.join(rows)))
    assert lc['oid'] == {'zg': 3, 'zr': 1}
    assert lc['t'][0] == pytest.approx(1200.5)                  # HJD - 2457000
    for flt in ('zg', 'zr'):                                    # each filter has median 1 despite a 2 mag offset
        assert np.median(lc['flux'][lc['filt'] == flt]) == pytest.approx(1.0)
    with pytest.raises(ztf.ZtfError):
        ztf.build_lightcurve(ztf.parse_table(CSV))              # too few points


def test_query_url_validation():
    url, p = ztf.query_url(10.0, -20.0, 0.001)
    assert p['POS'].startswith('CIRCLE 10.000000 -20.000000') and p['FORMAT'] == 'csv'
    assert ztf.query_url(oid=686103400034440)[1]['ID'] == ['686103400034440']
    for bad in (dict(ra=1, dec=1, radius_deg=0.5), dict(ra=1, dec=95, radius_deg=0.01), dict()):
        with pytest.raises(ztf.ZtfError):
            ztf.query_url(**bad)


def test_recovers_injected_transit_in_sparse_data():
    r = ztf.search(synth())
    assert r is not None
    assert abs(r['period'] - 3.7) < 0.01
    assert r['depth'] == pytest.approx(0.03, rel=0.3) and r['snr'] > 10
    assert r['n_in_transit'] >= 5 and r['n_transits'] >= 3


def test_pure_noise_is_not_a_strong_detection():
    r = ztf.search(synth(depth=0.0))
    assert r is None or r['snr'] < 8


def test_too_sparse_returns_none():
    assert ztf.search({k: v[:10] for k, v in synth().items()}) is None


def test_alias_flag():
    assert ztf.alias_flags(1.003) == 'alias_1d' and ztf.alias_flags(3.7) == ''


# ---------------------------------------------------------------- wiring (unit codes, command line)
def test_ztf_unit_codes_are_free_and_resolve():
    from occulta import sources
    us = sources.units('ZTF')
    assert us == [5001, 5002, 5003] == [ztf.UNIT_OF[f] for f in ztf.FILTERS]
    others = sum((sources.units(m) for m in ('QLP', 'SPOC', 'TESS-SPOC', 'K2', 'Kepler')), [])
    assert not set(us) & set(others)
    assert [sources.unit_name(u) for u in us] == ['ZTF zg', 'ZTF zr', 'ZTF zi']
    assert {sources.fmt_of_unit(u) for u in us} == {'ztf'}
    assert sources.unit_targets(5002) is None
    assert sources.unit_name(4003) == 'TESS S3 FFI SPOC' and sources.fmt_of_unit(4003) == 'spoc'


def test_ztf_command_prints_the_search(monkeypatch, capsys):
    from occulta import cli
    lc = synth()
    lc = {**lc, 'oid': {'zr': 1}, 'n_raw': len(lc['t']), 'n_good': len(lc['t'])}
    monkeypatch.setattr(ztf, 'search_position', lambda ra, dec, radius, **kw: {'lc': lc, 'result': {'period': 3.7, 'snr': 12.0}})
    cli.main(['ztf', '10.0', '20.0', '--band', 'zr'])
    out = capsys.readouterr().out
    assert 'units [5002]' in out and '"period": 3.7' in out


def test_ztf_command_reports_a_failed_request(monkeypatch):
    from occulta import cli
    def boom(*a, **k): raise ztf.ZtfError('IRSA request failed')
    monkeypatch.setattr(ztf, 'search_position', boom)
    with pytest.raises(SystemExit, match='ZTF: IRSA request failed'):
        cli.main(['ztf', '10.0', '20.0'])


def test_survey_plan_accepts_ztf_without_a_target_list(tmp_path, monkeypatch):
    from occulta import survey
    monkeypatch.setitem(survey.CFG, 'plan', [{'mission': 'ZTF'}])
    assert [j[3] for j in survey.jobs()] == [5001, 5002, 5003]


@pytest.mark.skipif(os.environ.get('OCCULTA_ONLINE') != '1', reason='needs IRSA; set OCCULTA_ONLINE=1')
def test_live_ztf_light_curve_near_wasp_12():
    """The public IRSA table at WASP-12 (Dec +29.7, in the ZTF footprint) has the columns the parser needs, and the
    search runs on the cleaned light curve of the object with the most points inside 10 arcsec (the field is crowded,
    so no period is asserted; this tests the access path only, and the query takes about one minute)."""
    out = ztf.search_position(97.6366, 29.6722, 10.0 / 3600)
    lc = out['lc']
    assert lc['n_good'] > 100 and len(lc['t']) >= 30 and set(lc['oid']) <= set(ztf.FILTERS)
    assert 1000 < lc['t'].min() < lc['t'].max() < 4000          # BTJD of 2018 to 2026
    r = out['result']
    assert r is None or (r['period'] > 0.2 and r['n_transits'] >= 3)
