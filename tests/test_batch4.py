"""Offline checks of the radius dilution correction, the Gaia match, the generalised-Gaussian eclipse fit, e cos w and
e sin w, the Roche ratio of the triage, the box durations and the false-alarm level, on small synthetic inputs."""
import json

import numpy as np
import pytest
from astropy.table import Table

from occulta import lightcurves as lc
from occulta.binaries import eclipses, gaia, report, triage
from test_cuts import SUN, _summary, _vet


# ---- planet radius with the TIC contamination ratio
@pytest.mark.parametrize('contr', [0.0, 0.5, 3.0, None, float('nan'), -0.2])
def test_planet_radius_is_corrected_for_the_tic_contamination_ratio(tmp_path, monkeypatch, contr):
    _vet(tmp_path, monkeypatch, star=SUN[:8] + (contr,))
    s = _summary(tmp_path)
    c = contr if contr is not None and np.isfinite(contr) and contr > 0 else 0.0
    base = np.sqrt(s['depth_ppm'] / 1e6) * 109.1                       # sqrt(depth) R*, 109.1 Re per Rsun
    assert abs(s['rp_earth_uncorrected'] - base) < 0.5
    assert s['rp_earth'] == pytest.approx(s['rp_earth_uncorrected'] * np.sqrt(1 + c))   # depth_true = depth_obs (1 + c)
    assert (s['rp_earth'] > s['rp_earth_uncorrected']) == (c > 0)
    report_md = next(tmp_path.rglob('report.md')).read_text(encoding='utf-8')
    assert 'without it' in report_md


# ---- Gaia match radius
def test_gaia_match_asks_for_the_brightest_source_within_3_arcsec(tmp_path, monkeypatch):
    (tmp_path / 'binaries' / 'eb').mkdir(parents=True)
    (tmp_path / 'findings' / 'TIC7').mkdir(parents=True)
    (tmp_path / 'binaries' / 'eb' / 'review.csv').write_text('tic,grade\n7,A\n8,C\n')
    (tmp_path / 'findings' / 'TIC7' / 'summary.json').write_text(json.dumps({'star': {'ra': 10.5, 'dec': -20.25}}))
    monkeypatch.setattr(gaia, 'B', tmp_path / 'binaries'); monkeypatch.setattr(gaia, 'ROOT', tmp_path)
    queries = []

    class Job:
        def get_results(self):
            return Table({'sep': [0.4, 2.1], 'source_id': [11, 22], 'gmag': [10.0, 14.0], 'bp_rp': [0.8, 1.2],
                          'parallax': [3.0, 1.0], 'parallax_error': [0.02, 0.1], 'ruwe': [1.1, 1.0],
                          'radial_velocity': [5.0, 6.0], 'rv_amplitude_robust': [1.0, 2.0], 'rv_nb_transits': [9, 8],
                          'non_single_star': [0, 0]})

    monkeypatch.setattr(gaia.Gaia, 'launch_job', lambda q, *a, **k: queries.append(q) or Job())
    gaia.main()
    q = queries[0]
    assert 'CIRCLE(ra,dec,3./3600)' in q and 'POINT(10.5,-20.25)' in q and q.endswith('order by phot_g_mean_mag')
    assert len(queries) == 1                                            # only the stars graded A or B are queried
    out = json.loads((tmp_path / 'binaries' / 'eb' / 'gaia_nss.json').read_text())
    assert out['7']['source_id'] == 11 and out['7']['gmag'] == 10.0     # the first row of the brightness order


# ---- generalised-Gaussian eclipse fit
def test_generalised_gaussian_fit_recovers_time_depth_and_width():
    rng = np.random.default_rng(5)
    t = np.arange(10, 10.8, 2 / 1440)
    tc, D, s, g, w = 10.403, 0.02, 0.05, 2.5, 0.12
    y = eclipses.model(t, 0.0, 0.002, tc, D, s, g) + rng.normal(0, 4e-4, t.size)
    r = eclipses.fit_one(t, y, tc + 0.01, w)                       # prediction 15 min off
    assert r is not None
    tc_f, tc_e, D_f, D_e, s_f, rms = r
    assert abs(tc_f - tc) < 4 * tc_e and tc_e < 5e-4
    assert abs(D_f - D) < 4 * D_e and abs(s_f - s) < 0.01
    assert 0.5 < rms * 1 / 4e-4 < 1.5                               # the residual is the noise level
    assert eclipses.fit_one(t[t < tc], y[t < tc], tc, w) is None    # an eclipse cut at its centre is not fitted


# ---- e cos w, e sin w
def test_orbit_solution_recovers_eccentricity_and_argument_of_periastron():
    assert report.phase_sep(0.0, 0.0) == pytest.approx(0.5)
    assert report.phase_sep(0.01, 0.0) == pytest.approx(0.5 + 2 * 0.01 / np.pi, rel=1e-3)   # first order: 2 e cos w / pi
    for e, w in [(0.2, 40), (0.05, 200), (0.4, 300), (0.1, 95)]:
        ecw, esw = e * np.cos(np.radians(w)), e * np.sin(np.radians(w))
        o = report.orbit(report.phase_sep(ecw, esw), (1 + esw) / (1 - esw))
        assert o['e_cos_w'] == pytest.approx(ecw, abs=1e-6) and o['e_sin_w'] == pytest.approx(esw, abs=1e-6)
        assert o['e'] == pytest.approx(e, abs=1e-6) and o['omega_deg'] == pytest.approx(w, abs=1e-4)
    assert report.orbit(1.2, 1.0) is None                        # no orbit with e < 0.99 gives this phase


# ---- Roche ratio in the triage
def test_roche_ratio_of_a_twin_mass_orbit():
    G, MSUN, RSUN = 6.674e-11, 1.989e30, 6.957e8
    P, M, R = 2.0, 1.0, 3.0
    a = (G * 2 * M * MSUN * (P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / RSUN      # Kepler's third law in Rsun
    assert triage.roche_ratio(R, M, P) == pytest.approx(R / (0.38 * a), rel=2e-3)
    assert 0.379 == pytest.approx(0.49 / (0.6 + np.log(2)), abs=5e-4)                   # Eggleton (1983) at q = 1
    assert triage.roche_ratio(R, float('nan'), P) == triage.roche_ratio(R, 1.0, P)      # no catalogue mass: 1 Msun
    assert triage.roche_ratio(2 * R, M, P) == pytest.approx(2 * triage.roche_ratio(R, M, P))


# ---- box durations
def test_box_search_uses_nine_widths_from_1_to_13_hours_and_matches_each_dip_width():
    assert lc.WIDTHS_H == [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.5, 9.0, 13.0]
    rng = np.random.default_rng(6)
    t = np.arange(1000, 1060, 10 / 1440)
    for width in (1.5, 3.0, 6.5, 13.0):
        y = 1 + rng.normal(0, 3e-4, t.size)
        y[np.abs(t - 1030) < width / 48] -= 5e-3
        tg, dt, (yg,) = lc.to_grid(t, y)
        top = max(lc.single_events(tg, dt, yg, lc.mad_std(yg)), key=lambda e: e['snr'])
        assert width / 1.6 <= top['dur_h'] <= width * 1.6, (width, top['dur_h'])
        assert 4.0e-3 < top['depth'] < 5.5e-3


# ---- false-alarm level
def test_false_alarm_level_is_the_strongest_inverted_dip_away_from_the_edges():
    rng = np.random.default_rng(7)
    t = np.arange(1000, 1026, 10 / 1440)
    tg, dt, _ = lc.to_grid(t, t)

    def fa(y, interior=None):
        yg = lc.to_grid(t, y)[2][0]
        return lc.single_events(tg, dt, yg, 400e-6, invert=True,
                                interior=np.ones(tg.size, bool) if interior is None else interior)

    noise = 1 + rng.normal(0, 400e-6, t.size)
    dip = np.where(np.abs(t - 1013) < 1.5 / 24, 3e-3, 0)
    assert 2 < fa(noise) < lc.SNR_EVENT                              # noise alone stays below the event threshold
    assert fa(noise - dip) > fa(noise)                               # a real dip raises it: boxes beside it see it in their flank
    assert fa(noise + dip) > 3 * lc.SNR_EVENT                        # a brightening is what the inverted search sees
    edge = np.ones(tg.size, bool); edge[np.abs(tg - 1013) < 0.5] = False
    assert fa(noise + dip, edge) < lc.SNR_EVENT                      # the same bump at a segment edge is ignored
