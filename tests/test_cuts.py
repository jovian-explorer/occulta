"""Offline checks of the remaining rejection cuts of the planet vetting (stellar size, V shape, step, depth change,
stellar density, pulsator), the Gaia checks, the cross-sector candidate rules, the common-mode flag, the pixel model,
and the light-travel-time fit. Archives are replaced by mocked answers and synthetic data of known content."""
import json
import sys
import types

import numpy as np
import pytest

from test_offline_science import _transits, _vet_synthetic

STARTS = (1000.0, 1030.0)
P, T0, DEP, DUR = 7.13, 1003.2, 2500e-6, 3 / 24
SUN = (10.0, 20.0, 10.0, 5700.0, 1.0, 1.0, 100.0, 'None', 0.0)       # ra, dec, Tmag, Teff, R, M, d, Gaia id, contamination
BLS = (P, T0, 3.0, DEP * 1e6)                                        # a BLS ephemeris stands in for the survey candidate


def _vet(tmp_path, monkeypatch, **kw):
    """One vetting in a fresh directory of tmp_path, so that a test can vet several variants."""
    d = tmp_path / f'v{len(list(tmp_path.glob("v*")))}'
    d.mkdir()
    return _vet_synthetic(d, monkeypatch, STARTS, **{'P': P, 't0': T0, 'depth': DEP, 'dur': DUR, **kw})


def _summary(tmp_path):
    return json.loads(next(tmp_path.glob('v*/findings/TIC1/summary.json')).read_text())


def _local(t):
    """Time from the nearest transit, and the epoch number."""
    n = np.round((t - T0) / P)
    return t - (T0 + n * P), n


# ---------------------------------------------------------------- stellar radius, depth and density
def test_planet_radius_follows_the_depth_and_the_tic_radius(tmp_path, monkeypatch):
    status, _, notes = _vet(tmp_path, monkeypatch, star=SUN)
    s = _summary(tmp_path)
    assert status == 'ready', notes
    assert abs(s['rp_earth'] - np.sqrt(s['depth_ppm'] / 1e6) * 109.1) < 0.5          # Rp = sqrt(depth) R*, 109.1 Re per Rsun


def test_vetting_rejects_a_stellar_size_radius_only_above_24_earth_radii(tmp_path, monkeypatch):
    big = SUN[:4] + (1.8, 1.0) + SUN[6:]          # 2 % depth on a 1.8 Rsun star: 27.8 Re
    status, _, notes = _vet(tmp_path, monkeypatch, depth=0.02, star=big)
    assert status == 'rejected' and 'stellar-size' in notes, notes
    ok = SUN[:4] + (1.4, 1.0) + SUN[6:]           # the same depth on a 1.4 Rsun star: 21.6 Re, below the limit
    status, _, notes = _vet(tmp_path, monkeypatch, depth=0.02, star=ok)
    assert 'stellar-size' not in notes, notes


def test_vetting_rejects_a_deep_dip_on_a_star_without_radius(tmp_path, monkeypatch):
    star = SUN[:4] + (None, None) + SUN[6:]
    status, _, notes = _vet(tmp_path, monkeypatch, depth=0.08, star=star)
    assert status == 'rejected' and 'no stellar radius' in notes, notes


@pytest.mark.parametrize('rad, expect', [(0.4, 'too long for the star at this period'), (0.2, 'far too long for the star')])
def test_vetting_rejects_a_transit_too_long_for_the_stellar_density(tmp_path, monkeypatch, rad, expect):
    """A 3 h transit at 7.13 d needs 2.2 g/cc (Seager & Mallen-Ornelas 2003); a dwarf of 0.4 or 0.2 Rsun and 1 Msun
    has 22 or 176 g/cc, so that the transit lasts far too long (ratio 0.10 and 0.013)."""
    star = SUN[:4] + (rad, 1.0) + SUN[6:]
    status, _, notes = _vet(tmp_path, monkeypatch, star=star)
    assert status == 'rejected' and expect in notes, notes


def test_transit_density_matches_the_central_transit_formula(tmp_path, monkeypatch):
    _vet(tmp_path, monkeypatch, star=SUN)
    s = _summary(tmp_path)
    k = np.sqrt(s['depth_ppm'] / 1e6); ar = P / (np.pi * s['t14_h'] / 24) * (1 + k)
    rho = 3 * np.pi * ar ** 3 / (6.674e-11 * (P * 86400) ** 2) / 1000
    assert abs(s['rho_transit_b0'] / rho - 1) < 0.01 and 1.5 < s['rho_transit_b0'] < 3.5     # about 2.2 g/cc


# ---------------------------------------------------------------- shapes and variability (through a BLS ephemeris)
def test_vetting_rejects_a_v_shaped_transit(tmp_path, monkeypatch):
    def tri(t):                                   # a triangle of 4 h and 3000 ppm instead of a flat-bottomed transit
        x, _ = _local(t)
        return np.clip(1 - np.abs(x) / (2 / 24), 0, None) * -3000e-6
    status, _, notes = _vet(tmp_path, monkeypatch, depth=0.0, extra=tri, cand=(P, T0, 4.0, 3000.0))
    assert status == 'rejected' and 'V-shaped' in notes, notes


def test_vetting_rejects_a_step_instead_of_a_transit(tmp_path, monkeypatch):
    def step(t):                                  # the flux stays one depth lower after every transit, until half a period later
        x, _ = _local(t)
        return -DEP * (x > 0)
    status, _, notes = _vet(tmp_path, monkeypatch, extra=step, cand=BLS)
    assert status == 'rejected' and 'step or ramp' in notes, notes


def test_vetting_rejects_depths_that_change_between_sectors(tmp_path, monkeypatch):
    def deeper_first(t):                          # 4000 ppm in the first sector, 1300 ppm in the second
        x, _ = _local(t)
        return np.where((t < 1026) & (np.abs(x) < DUR / 2), -2700e-6, 0.0)
    status, _, notes = _vet(tmp_path, monkeypatch, depth=1300e-6, extra=deeper_first, cand=(P, T0, 3.0, 2500.0))
    assert status == 'rejected' and 'depth changes between epochs' in notes, notes


def test_vetting_rejects_a_pulsator_but_not_a_weak_oscillation(tmp_path, monkeypatch):
    osc = lambda amp: (lambda t: amp * np.sin(2 * np.pi * t / (1 / 24)))      # one-hour period, above 1.5 / T14
    status, _, notes = _vet(tmp_path, monkeypatch, extra=osc(0.7 * DEP), cand=BLS)
    assert 'pulsator' in notes, notes
    status, _, notes = _vet(tmp_path, monkeypatch, extra=osc(0.2 * DEP), cand=BLS)
    assert 'pulsator' not in notes, notes
