"""Offline checks of the microlensing search for free-floating planets (rogues.py): the point-lens magnification, the
recovery of injected events of known parameters, and the rejection of flares, moving sources, quality-flagged cadences,
pure noise and brightenings common to one CCD."""
import numpy as np
import pytest

from occulta import rogues as R

CAD = 10 / 1440


def _lc(seed, noise=1000e-6, span=27.0):
    rng = np.random.default_rng(seed)
    t = np.arange(2000.0, 2000.0 + span, CAD)
    t = t[(t < 2013.0) | (t > 2014.0)]                    # the orbit gap of a TESS sector
    return t, 1 + rng.normal(0, noise, t.size), rng


def test_magnification_point_and_finite_source():
    assert R.magnification(np.array([1.0]))[0] == pytest.approx(3 / np.sqrt(5), rel=1e-9)
    rho = 0.1
    assert R.magnification(np.array([0.0]), rho)[0] == pytest.approx(np.sqrt(1 + 4 / rho ** 2), rel=2e-3)
    # far from the source the finite source does not matter
    assert R.magnification(np.array([1.0]), 0.01)[0] == pytest.approx(3 / np.sqrt(5), rel=1e-3)


@pytest.mark.parametrize('te, u0, fs', [(0.15, 0.3, 0.02), (0.5, 0.3, 0.02), (0.5, 1.0, 0.03), (1.0, 0.5, 0.01)])
def test_injected_lens_is_recovered(te, u0, fs):
    t, f, _ = _lc(int(te * 100 + u0 * 10))
    t0 = 2007.3
    f = f + R.lens(t, t0, te, u0, fs)
    ev = [e for e in R.search_lc(t, f) if abs(e['t0'] - t0) < 0.1]
    assert ev, 'injected event not triggered'
    e = ev[0]
    assert e['flags'] == '', e['flags']
    assert e['te_d'] == pytest.approx(te, rel=0.5)        # tE, u0 and fs trade off in weak, blended events
    amp = fs * (R.magnification(np.array([u0]))[0] - 1)
    assert e['amp_ppm'] == pytest.approx(amp * 1e6, rel=0.25)


def test_flare_is_rejected():
    t, f, _ = _lc(3)
    f = f + R.flare(t, 2007.3, 0.02, 0.03)
    ev = [e for e in R.search_lc(t, f) if abs(e['t0'] - 2007.3) < 0.2]
    assert ev and all(('flare' in e['flags']) or ('asymmetric' in e['flags']) for e in ev)


def test_moving_source_is_rejected():
    """An asteroid crossing the aperture brightens the star and drags the centroid along."""
    t, f, rng = _lc(4)
    t0, te = 2007.3, 0.2
    f = f + R.lens(t, t0, te, 0.5, 0.02)
    bump = R.lens(t, t0, te, 0.5, 1.0) / R.lens(np.array([t0]), t0, te, 0.5, 1.0)[0]
    cx = 0.002 * rng.normal(size=t.size) + 0.05 * bump * np.tanh((t - t0) / te)
    cy = 0.002 * rng.normal(size=t.size)
    ev = [e for e in R.search_lc(t, f, cx=cx, cy=cy) if abs(e['t0'] - t0) < 0.1]
    assert ev and 'moving' in ev[0]['flags']


def test_quality_flagged_cadence_is_rejected():
    t, f, _ = _lc(5)
    t0 = 2007.3
    f = f + R.lens(t, t0, 0.3, 0.4, 0.02)
    ev = [e for e in R.search_lc(t, f, flagged=[t0 + 0.01]) if abs(e['t0'] - t0) < 0.1]
    assert ev and 'quality' in ev[0]['flags']


@pytest.mark.parametrize('seed', range(4))
def test_noise_gives_no_candidate(seed):
    t, f, _ = _lc(100 + seed)
    assert not [e for e in R.search_lc(t, f) if e['flags'] == '']


def test_common_mode_on_one_ccd_is_rejected(tmp_path):
    from occulta import survey
    con = survey.db(str(tmp_path / 'hunt.db'))
    row = dict(te_d=0.2, u0=0.5, fs=0.02, amp_ppm=20000, snr=20, dchi2_flat=500, dchi2_flare=50, asym_sigma=0,
               motion_sigma=0, offset_sigma=0, bkg_sigma=0, chi2r=1, flags='')
    for i in range(6):                                    # six stars on camera 1, CCD 2 within 0.05 d: systematics
        R.store(con, 100 + i, 10, 1, 2, [dict(row, t0=2007.30 + 0.01 * i)], tmp_path)
    R.store(con, 200, 10, 3, 4, [dict(row, t0=2007.30)], tmp_path)     # a lone event elsewhere
    cands = R.rank(con, tmp_path, common_n=5)
    assert [c['tic'] for c in cands] == [200]
    assert con.execute("SELECT COUNT(*) FROM brightening WHERE flags LIKE '%common%'").fetchone()[0] == 6
    assert (tmp_path / 'rogues' / 'report.md').exists() and (tmp_path / 'rogues' / 'rogues.csv').exists()
