"""Transit-shape fit of classify: recovery of the stellar density from synthetic transits, and the decisions."""
import numpy as np
import pytest

from occulta import classify as c

batman = pytest.importorskip('batman')


def synth(P=3.0, ar=10.0, k=0.1, b=0.3, noise=5e-4, cad=120.0, seed=1, n=4):
    t = np.arange(-1, (n - 1) * P + 1, cad / 86400)                  # one continuous grid, windows of 1 d around each transit
    t = t[np.abs((t + P / 2) % P - P / 2) < 1]
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc = 0.0, P, k, ar, np.degrees(np.arccos(b / ar))
    p.ecc, p.w, p.u, p.limb_dark = 0.0, 90.0, [0.4, 0.25], 'quadratic'
    f = batman.TransitModel(p, t, supersample_factor=max(1, round(cad / 240)), exp_time=cad / 86400).light_curve(p) - 1
    return t, f + np.random.default_rng(seed).normal(0, noise, len(t))


@pytest.mark.parametrize('P,ar,b,cad,n', [(3.0, 10.0, 0.3, 120.0, 4), (3.1234, 15.0, 0.5, 1800.0, 25)])
def test_fit_recovers_density(P, ar, b, cad, n):
    k = 0.1
    t, f = synth(P, ar, k, b, cad=cad, n=n)
    dur = float(c.t14_h(k, ar, b, P))
    r = c.transit_fit(t, f, P, 0.0, dur, k * k * 1e6, u=(0.4, 0.25), n_boot=30)
    assert r['ok']
    lo, hi = np.percentile(r['rho'], [0.5, 99.5])
    assert lo < c.rho_from_ar(ar, P) < hi
    assert abs(np.median(r['k']) - k) < 0.01


def test_fit_refuses_uncovered_transit():
    t, f = synth()
    keep = np.abs((t + 1.5) % 3.0 - 1.5) > 0.02           # remove the centre of every transit
    r = c.transit_fit(t[keep], f[keep], 3.0, 0.0, 2.4, 1e4, n_boot=5)
    assert not r['ok']


def test_kepler_density_round_trip():
    assert np.isclose(c.rho_from_ar(c.ar_from_rho(1.41, 10.0), 10.0), 1.41)


def test_tables():
    u1, u2 = c.ld_coeffs(5750, 4.5, 'TESS')
    assert 0.2 < u1 < 0.5 and 0.1 < u2 < 0.4
    assert c.ld_coeffs(5750, 4.5, 'Kepler')[0] > u1           # stronger limb darkening in the bluer Kepler band
    assert 1.0 < c.ms_rho(5770) < 1.8                          # the Sun on the dwarf sequence
    assert 5500 < c.ms_teff(0.82) < 6000
    assert np.isclose(c.ms_radius(c.ms_rho(5770)), 1.0, atol=0.1)
    assert np.isnan(c.ms_radius(1e4))


def fake(k, b, rho, n=200, seed=0):
    g = np.random.default_rng(seed)
    return {'k': np.full(n, k) * (1 + 0.01 * g.normal(size=n)), 'b': np.full(n, b) + 0.005 * g.normal(size=n),
            'rho': np.full(n, rho) * (1 + 0.05 * g.normal(size=n))}


def test_class_from_rp():
    assert c.class_from_rp(np.full(100, 1.0)) == 'planet'
    assert c.class_from_rp(np.full(100, 2.5)) == 'brown dwarf or low-mass star'
    assert c.class_from_rp(np.full(100, 5.0)) == 'eclipsing binary'
    assert c.class_from_rp(np.r_[np.full(50, 1.9), np.full(50, 2.1)]) is None


def test_long_grazing_is_eb():
    assert c.shape_long(fake(0.3, 1.05, 1.0), 1.4, 1.0)[0] == 'eclipsing binary'


def test_long_low_density_is_eb():
    t, note = c.shape_long(fake(0.1, 0.2, 0.05), 1.4, 1.0)
    assert t == 'eclipsing binary' and 'e < 0.6' in note


def test_long_consistent_flat_fit_is_retested():
    t, _ = c.shape_long(fake(0.1, 0.2, 1.3), 1.4, 1.0)
    assert t == 'planet'                                       # Rp = 0.1 R_sun = 0.97 R_J


def test_long_eccentric_range_stays_undetermined():
    assert c.shape_long(fake(0.1, 0.2, 0.4), 1.4, 1.0)[0] is None   # within e < 0.6 of the host, but not equal


def test_norad_dwarf():
    rho = c.ms_rho(3300)
    t, note, rs = c.shape_norad(fake(0.1, 0.3, rho), 3300, 'catalogue')
    assert t == 'planet' and rs < 0.4


def test_norad_inconsistent_with_dwarf():
    assert c.shape_norad(fake(0.1, 0.3, 0.01), 5800, 'catalogue')[0] is None
    assert c.shape_norad(fake(0.1, 0.3, 1.0), None, '')[0] is None
