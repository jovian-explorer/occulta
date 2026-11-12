"""Transit duration against the host density (occulta.classify.duration_test)."""
from occulta import classify as c


def test_density_of_the_sun_and_unknown_hosts():
    assert abs(c.stellar_density(1.0, 1.0) - 1.41) < 0.01
    assert c.stellar_density(None, 1.0) is None and c.stellar_density(1.0, 0) is None


def test_earth_transit_of_the_sun_lasts_13_hours():
    assert abs(c.max_duration_h(365.25, c.stellar_density(1, 1)) - 13.0) < 0.2


def test_companion_radius_lengthens_the_maximum():
    rho = c.stellar_density(1, 1)
    assert c.max_duration_h(3.0, rho, 10000) > 1.09 * c.max_duration_h(3.0, rho, 0)


def test_contact_orbit_is_capped_at_half_the_period():
    assert abs(c.max_duration_h(1.0, 1e-4) - 12.0) < 1e-9


def test_hot_jupiter_duration_passes_and_a_long_one_fails():
    rho = c.stellar_density(1, 1)
    ok, note = c.duration_test(2.8, 3.0, rho, 10000)
    assert not ok and 'longest central transit' in note
    bad, _ = c.duration_test(8.0, 3.0, rho, 10000)
    assert bad
    assert not c.duration_test(8.0, 3.0, rho, 10000, factor=3.0)[0]


def test_no_density_is_not_flagged():
    assert c.duration_test(20.0, 0.5, None) == (False, 'duration test: no density')
