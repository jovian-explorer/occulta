"""Tests of the host-radius, spectroscopic-mass, eccentricity, any-phase secondary and phase-curve checks of classify."""
import numpy as np

from occulta import classify as c


def test_bc_g_solar_and_continuous():
    assert abs(float(c.bc_g(5772.0)) - 0.06) < 1e-9
    assert abs(float(c.bc_g(3999.9)) - float(c.bc_g(4000.0))) < 0.02     # the two Andrae et al. pieces meet at 4000 K


def test_gaia_radius_sun_at_10pc():
    g = {'plx': 100.0 + c.PLX_ZP_MAS, 'plx_e': 0.05, 'gmag': 4.67, 'ag': None}
    r, note = c.gaia_radius(g, 5772.0, 50.0)
    assert 0.95 < np.median(r) < 1.05 and 'lower limit' in note


def test_gaia_radius_needs_parallax_snr():
    r, note = c.gaia_radius({'plx': 0.4, 'plx_e': 0.2, 'gmag': 12.0}, 5000.0, 150.0)
    assert r is None and 'below' in note


def test_giant_is_evolved():
    # a 4800 K giant of 10 R_sun at 1 kpc: L = 10^2 (4800/5772)^4 = 47.8 L_sun
    lum = 100 * (4800 / 5772) ** 4
    mbol = c.MBOL_SUN - 2.5 * np.log10(lum)
    gm = mbol - float(c.bc_g(4800.0)) + 10.0
    r, _ = c.gaia_radius({'plx': 1.0 + c.PLX_ZP_MAS, 'plx_e': 0.02, 'gmag': gm, 'ag': 0.0}, 4800.0, 100.0)
    assert 9 < np.median(r) < 11 and c.host_state(r, 4800.0) == 'evolved'
    assert c.host_state(np.full(100, 0.8), 5500.0) == 'dwarf'


def test_orbit_inside_giant():
    assert 'orbit inside' in c.orbit_inside(np.full(100, 119.0), 70.0)      # a = 103 R_sun for 3 M_sun
    assert c.orbit_inside(np.full(100, 1.0), 10.0) == ''


def test_spec_mass():
    assert abs(c.spec_mass(4.438, 1.0) - 1.0) < 1e-9
    assert abs(c.spec_mass(2.438, 10.0) - 1.0) < 1e-9


def test_ecc_from_secondary_and_density_factor():
    ecw, esw, e = c.ecc_from_secondary(0.5, 1.0)
    assert ecw == 0 and esw == 0 and e == 0
    ecw, esw, e = c.ecc_from_secondary(0.6, 1.5)
    assert abs(ecw - np.pi / 20) < 1e-12 and abs(esw - 0.2) < 1e-12
    assert abs(c.rho_ecc_factor(-0.6, 0.6) - 0.125) < 1e-12     # apoastron transit at e = 0.6, as RHO_EB_RATIO


def test_occultation_max_hot_jupiter():
    # k 0.1, a/R 5 around a 6000 K star: reflection 400 ppm plus a ~2000 K dayside in the red band
    d = c.occultation_max(0.1, 5.0, 6000.0)
    assert 4e-4 < d < 2e-3


def _lc(P=7.3, dur_h=4.0, dep=0.01, sec_phase=None, sec_dep=0.0, sec_w=1.0, noise=3e-4, seed=0):
    rng = np.random.default_rng(seed)
    T = np.arange(0, 120, 10 / 1440)
    Y = rng.normal(0, noise, len(T))
    ph = (T / P) % 1
    d = np.minimum(ph, 1 - ph) * P * 24
    Y[d < dur_h / 2] -= dep
    if sec_phase is not None:
        Y[np.abs(ph - sec_phase) * P * 24 < sec_w * dur_h / 2] -= sec_dep
    return T, Y


def test_secondary_scan_finds_eccentric_eclipse():
    T, Y = _lc(sec_phase=0.31, sec_dep=0.003, sec_w=1.5)
    s = c.secondary_scan(T, Y, 7.3, 0.0, 4.0, 1e4)
    assert s['detected'] and s['stellar'] and abs(s['phase'] - 0.31) < 0.01
    assert s['ecosw'] < -0.2 and s['esinw'] > 0.1
    ok, note = c.lc_test(s, {})
    assert ok and 'secondary eclipse at phase' in note


def test_secondary_scan_quiet_planet():
    T, Y = _lc(dep=0.005)
    s = c.secondary_scan(T, Y, 7.3, 0.0, 4.0, 5e3)
    assert s['ok'] and not s['detected']
    assert c.lc_test(s, {})[0] is False


def test_secondary_scan_ignores_single_event():
    T, Y = _lc(dep=0.005)
    Y[np.abs(T - 40.0) < 0.1] -= 0.02          # one lone 2 per cent dip (systematic or a neighbour's single eclipse)
    s = c.secondary_scan(T, Y, 7.3, 0.0, 4.0, 5e3)
    assert not s['detected']


def test_phase_curve_ellipsoidal_star_and_long_period_skip():
    rng = np.random.default_rng(3)
    P, T = 0.8, np.arange(0, 27, 2 / 1440)
    ph = (T / P) % 1
    Y = -2e-3 * np.cos(4 * np.pi * ph) + rng.normal(0, 5e-4, len(T))
    pc = c.phase_curve(T, Y, P, 0.0, 1.5, 2e4, 1.0)
    assert pc['ok'] and pc['detected'] and pc['stellar'] and abs(pc['ell_ppm'] - 2000) < 100
    assert not c.phase_curve(T, Y, 3.0, 0.0, 3.0, 1e4)['ok']
