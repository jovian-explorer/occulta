"""Offline checks of the moon search (moons.py) on synthetic transits: a moon's TTV and TDV in quadrature, a TTV alone,
TTV and TDV in phase, an orbital-sampling dip, pure noise with its mass and radius limits, and the files and table
that one search writes."""
import json

import numpy as np
import pytest

from occulta import moons as M
from occulta.planets import trapezoid_ss

CAD = 10 / 1440


def _transits(seed, P=5.0, N=30, ttv_min=0.0, tdv_min=0.0, f=0.27, dphi=90.0, noise=500e-6, moon=None):
    """Transit windows of +-1 d around N transits (depth 5000 ppm, 3 h) with sinusoidal TTV and TDV of the given
    amplitudes (min), frequency (cycles per transit) and phase offset, or with a moon dip (projected orbit radius in days, depth)
    at a random orbital phase in every transit."""
    rng = np.random.default_rng(seed)
    t0, t14, dep = 1000.0, 3 / 24, 5000e-6
    T, Y = [], []
    for n in range(N):
        tc = t0 + n * P + ttv_min / 1440 * np.cos(2 * np.pi * f * n)
        d = t14 + tdv_min / 1440 * np.cos(2 * np.pi * f * n - np.radians(dphi))
        t = np.arange(tc - 1, tc + 1, CAD)
        y = trapezoid_ss(t - tc, 0, dep, d, 0.1 * d, CAD) + rng.normal(0, noise, t.size)
        if moon:
            ts, dm = moon
            x = ts * np.cos(rng.uniform(0, 2 * np.pi))
            y = y - np.where(np.abs(t - tc - x) < d / 2, dm, 0)
        T.append(t); Y.append(y)
    s = {'tic': 1, 'period_d': P, 't0_bjd': t0 + 2457000, 't14_h': 3.0, 'depth_ppm': 5000, 'ingress_frac': 0.2,
         'rp_earth': 8.0, 'star': {'mass': 1.0, 'rad': 1.0}, 'status': 'ready'}
    return np.concatenate(T), np.concatenate(Y), s


def test_planet_mass_relation():
    assert M.planet_mass(1.0) == pytest.approx(1.0, rel=0.05)
    assert M.planet_mass(3.9) == pytest.approx(13.0, rel=0.35)              # Neptune-size: about 10-20 Earth masses
    assert M.planet_mass(15.0) == 317.8                                     # above the Neptunian branch


@pytest.mark.parametrize('seed', range(3))
def test_moon_timing_in_quadrature_is_found(seed):
    o = M.analyse(*_transits(seed, ttv_min=4, tdv_min=6))
    assert o['status'] == 'moon_signal' and 'moon_timing' in o['signals']
    assert abs(abs(o['phase']['dphi_deg']) - 90) < 30
    assert o['phase']['freq_per_epoch'] == pytest.approx(0.27, abs=0.01)


def test_ttv_alone_is_a_timing_variation_not_a_moon():
    o = M.analyse(*_transits(1, ttv_min=4))
    assert o['status'] == 'timing_variation' and 'moon_timing' not in o['signals']


def test_ttv_and_tdv_in_phase_are_not_a_moon():
    o = M.analyse(*_transits(2, ttv_min=4, tdv_min=6, dphi=0.0))
    assert 'moon_timing' not in o['signals']


@pytest.mark.parametrize('seed', range(3))
def test_noise_gives_no_moon_and_finite_limits(seed):
    o = M.analyse(*_transits(10 + seed))
    assert o['status'] == 'no_moon'
    lim = o['mass_limit']['ms_earth_any_orbit']
    assert np.isfinite(lim) and 0 < lim < 1e3
    assert o['ttv']['limit95_s'] > 0 and o['tdv']['limit95_s'] > 0


def test_orbital_sampling_dip_is_found_and_absent_in_noise():
    o = M.analyse(*_transits(3, P=40.0, N=12, moon=(0.08, 1500e-6)))
    assert 'ose_dip' in o['signals'] and o['status'] == 'moon_signal'
    o = M.analyse(*_transits(4, P=40.0, N=12))
    assert 'ose_dip' not in o['signals']
    assert 0 < o['ose']['rs_limit_earth_any_orbit'] < 5


def test_close_in_planet_has_no_ose_stack_and_says_so():
    o = M.analyse(*_transits(5))
    assert 'ose' not in o and any('OSE' in n for n in o['notes'])


def test_too_few_transits():
    o = M.analyse(*_transits(6, N=3))
    assert o['status'] == 'too_few_transits'


def test_run_one_writes_report_and_table(tmp_path):
    from occulta import survey
    t, y, s = _transits(7)
    d = tmp_path / 'TIC1'
    d.mkdir()
    (d / 'summary.json').write_text(json.dumps(s))
    np.savez_compressed(d / 'transit_windows.npz', t=t, f=y, period=s['period_d'], t0=s['t0_bjd'])
    con = survey.db(str(tmp_path / 'hunt.db'))
    assert [p.name for p in M.folders(tmp_path)] == ['TIC1']
    assert M.run_one(con, d) == 'no_moon'
    for f in ('moons.md', 'moons.json', 'moons.png'):
        assert (d / f).exists()
    assert con.execute('SELECT status, n_transits FROM moon WHERE tic=1').fetchone() == ('no_moon', 30)
