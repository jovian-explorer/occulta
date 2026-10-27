"""Offline checks of the thresholds of the detrending, the single-event search, and the planet vetting, on small synthetic
light curves with a known signal."""
import numpy as np

from occulta import lightcurves as lc
from test_offline_science import _transits, _vet_synthetic

STARTS = (1000.0, 1030.0)
PLANET = dict(P=7.13, t0=1003.2, depth=2500e-6, dur=3 / 24)


def _flux(rng, t, dips=(), noise=400e-6):
    """Relative flux with white noise and box dips (centre, depth, duration in days)."""
    y = 1 + rng.normal(0, noise, t.size)
    for c, d, w in dips:
        y[np.abs(t - c) < w / 2] -= d
    return y


def _clean(ev, tg, dt, yg, t):
    """Events that carry no flag: a box at the very start or end of a series is flagged 'edge' by characterise."""
    return [e for e in ev if lc.characterise(e, tg, dt, yg, yg * 0, yg * 0, yg * 0, [(t[0], t[-1])])['flags'] == '']


def test_detrend_removes_a_slow_trend_and_keeps_the_dip():
    rng = np.random.default_rng(1)
    t = np.arange(1000, 1026, 10 / 1440)
    f = _flux(rng, t, [(1013.0, 3000e-6, 3 / 24)]) * (1 + 0.01 * np.sin(2 * np.pi * (t - 1000) / 12))
    y, segs = lc.detrend(t, f * 1000)
    assert len(segs) == 1 and np.nanstd(y) < 600e-6                 # the 1 % trend is gone
    assert 2200e-6 < -np.nanmean(y[np.abs(t - 1013.0) < 1 / 24]) < 3300e-6    # a 3-h dip survives a 2-d window


def test_detrend_works_per_segment_across_a_gap_and_clips_only_upward_outliers():
    rng = np.random.default_rng(2)
    t = np.r_[np.arange(1000, 1010, 10 / 1440), np.arange(1011, 1021, 10 / 1440)]     # a 1-d gap splits the series
    f = _flux(rng, t) * np.where(t < 1010.5, 1.0, 1.02)                                # a step across the gap
    f[100] *= 1.05; f[200] *= 0.95
    y, segs = lc.detrend(t, f)
    assert len(segs) == 2 and lc.mad_std(y) < 600e-6      # the 2 % step across the gap is gone
    assert np.isnan(y[100]) and np.isfinite(y[200]) and y[200] < -0.03       # flares are clipped, dips are kept


def test_single_event_search_finds_a_dip_above_the_threshold():
    rng = np.random.default_rng(3)
    t = np.arange(1000, 1026, 10 / 1440)
    y, _ = lc.detrend(t, _flux(rng, t, [(1013.0, 2500e-6, 3 / 24)]))
    tg, dt, (yg,) = lc.to_grid(t, y)
    ev = lc.single_events(tg, dt, yg, lc.mad_std(yg))
    ev = _clean(ev, tg, dt, yg, t)
    assert len(ev) == 1 and ev[0]['snr'] >= lc.SNR_EVENT
    c = lc.characterise(ev[0], tg, dt, yg, yg * 0, yg * 0, yg * 0, [(t[0], t[-1])])
    assert abs(c['t0'] - 1013.0) < 0.05 and 2000 < c['depth_ppm'] < 3000 and c['flags'] == ''


def test_single_event_search_ignores_noise_and_dips_below_the_threshold():
    rng = np.random.default_rng(4)
    t = np.arange(1000, 1026, 10 / 1440)
    for dips in ([], [(1013.0, 150e-6, 3 / 24)]):                   # pure noise; a 150 ppm dip is below SNR 7
        y, _ = lc.detrend(t, _flux(rng, t, dips))
        tg, dt, (yg,) = lc.to_grid(t, y)
        assert _clean(lc.single_events(tg, dt, yg, lc.mad_std(yg)), tg, dt, yg, t) == []


def test_single_event_search_reports_each_of_several_dips_once():
    rng = np.random.default_rng(5)
    t = np.arange(1000, 1026, 10 / 1440)
    y, _ = lc.detrend(t, _flux(rng, t, [(1006.0, 3000e-6, 3 / 24), (1013.0, 3000e-6, 3 / 24), (1020.0, 3000e-6, 3 / 24)]))
    tg, dt, (yg,) = lc.to_grid(t, y)
    ev = lc.single_events(tg, dt, yg, lc.mad_std(yg))
    assert sorted(round(tg[e['j']] + e['k'] * dt / 2) for e in ev) == [1006, 1013, 1020]


def test_characterise_flags_a_v_shaped_dip_and_a_dip_at_a_gap_edge():
    rng = np.random.default_rng(6)
    t = np.arange(1000, 1026, 10 / 1440)
    y = _flux(rng, t, noise=300e-6) - 1 - 1e-2 * np.clip(1 - np.abs(t - 1013.0) / (0.5 / 24), 0, None) ** 2     # a spike
    tg, dt, (yg,) = lc.to_grid(t, y)
    ev = lc.single_events(tg, dt, yg, lc.mad_std(yg))
    assert ev and 'vshape' in lc.characterise(ev[0], tg, dt, yg, yg * 0, yg * 0, yg * 0, [(t[0], t[-1])])['flags']
    assert 'edge' in lc.characterise(ev[0], tg, dt, yg, yg * 0, yg * 0, yg * 0, [(t[0], 1013.1), (1013.3, t[-1])])['flags']


def test_characterise_flags_a_background_excursion():
    rng = np.random.default_rng(7)
    t = np.arange(1000, 1026, 10 / 1440)
    y, _ = lc.detrend(t, _flux(rng, t, [(1013.0, 3000e-6, 3 / 24)]))
    tg, dt, (yg,) = lc.to_grid(t, y)
    bkg = rng.normal(0, 1, tg.size) + 30 * (np.abs(tg - 1013.0) < 1.5 / 24)
    ev = lc.single_events(tg, dt, yg, lc.mad_std(yg))
    assert 'bkg' in lc.characterise(ev[0], tg, dt, yg, bkg, yg * 0, yg * 0, [(t[0], t[-1])])['flags']


def _vet(tmp_path, monkeypatch, **kw):
    return _vet_synthetic(tmp_path, monkeypatch, STARTS, **{**PLANET, **kw})


def test_vetting_accepts_a_clean_planet(tmp_path, monkeypatch):
    status, period, notes = _vet(tmp_path, monkeypatch)
    assert status == 'ready' and abs(period - 7.13) < 0.01, (status, notes)


def test_vetting_rejects_a_stellar_depth(tmp_path, monkeypatch):
    status, _, notes = _vet(tmp_path, monkeypatch, depth=0.3)
    assert status == 'rejected' and 'too deep' in notes, notes


def test_vetting_rejects_a_secondary_eclipse(tmp_path, monkeypatch):
    sec = lambda t: _transits(t, 7.13, 1003.2 + 7.13 / 2, 500e-6, 3 / 24)
    status, _, notes = _vet(tmp_path, monkeypatch, extra=sec)
    assert status == 'rejected' and 'secondary eclipse' in notes, notes


def test_vetting_rejects_alternating_depths(tmp_path, monkeypatch):
    def odd(t):                  # every second transit twice as deep
        n = np.round((t - 1003.2) / 7.13).astype(int)
        return np.where(n % 2 == 0, _transits(t, 7.13, 1003.2, 2500e-6, 3 / 24), 0.0)
    status, _, notes = _vet(tmp_path, monkeypatch, extra=odd)
    assert status == 'rejected' and 'odd/even' in notes, notes


def test_vetting_rejects_a_signal_with_fewer_than_two_detectable_dips(tmp_path, monkeypatch):
    status, _, notes = _vet(tmp_path, monkeypatch, depth=300e-6)      # below the single-event threshold
    assert status == 'rejected' and 'fewer than two consistent dips' in notes, notes


def test_vetting_halves_a_bls_period_that_is_twice_the_true_one(tmp_path, monkeypatch):
    """LHS 3844 b: BLS locks on 2 x P; the in-between epochs show the same depth, so P/2 is adopted."""
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, STARTS, 1.2, 1000.4, 2500e-6, 2 / 24,
                                           cand=(2.4, 1000.4, 2.0, 2500))
    assert abs(period - 1.2) < 0.005, (status, period, notes)


def test_vetting_does_not_halve_alternating_depths(tmp_path, monkeypatch):
    """A signal of period P whose in-between epochs have half the depth (an eclipsing binary) keeps P."""
    half = lambda t: _transits(t, 2.4, 1001.6, 1250e-6, 2 / 24)
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, STARTS, 2.4, 1000.4, 2500e-6, 2 / 24, extra=half,
                                           cand=(2.4, 1000.4, 2.0, 2500))
    assert abs(period - 2.4) < 0.005, (status, period, notes)


def test_vetting_keeps_the_true_bls_period(tmp_path, monkeypatch):
    status, period, notes = _vet(tmp_path, monkeypatch, cand=(7.13, 1003.2, 3.0, 2500))
    assert status == 'ready' and abs(period - 7.13) < 0.01, (status, period, notes)


LONG = (1000.0, 1700.0, 2400.0, 3000.0)      # four sectors over a 2,000 d baseline


def test_period_refit_recovers_a_bls_period_that_is_1e4_off(tmp_path, monkeypatch):
    """A BLS period 1e-4 d too long moves the last transit by about 0.17 d over 2,000 d; the refit brings it back."""
    from occulta import planets
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, LONG, 1.2, 1000.4, 1500e-6, 2 / 24,
                                           noise=1200e-6, cand=(1.2001, 1000.4, 2.0, 1500))
    assert status == 'ready' and abs(period - 1.2) < 2e-5, (status, period, notes)


def test_period_refit_with_the_doubled_period_off_by_2e4(tmp_path, monkeypatch):
    """LHS 3844 b case: BLS gives 2 x P with an error of 2e-4 d; P/2 is still adopted and exact."""
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, LONG, 1.2, 1000.4, 1500e-6, 2 / 24,
                                           noise=1200e-6, cand=(2.4002, 1000.4, 2.0, 1500))
    assert status == 'ready' and abs(period - 1.2) < 2e-5, (status, period, notes)


def test_refine_period_is_a_fixed_small_scan():
    from occulta import planets
    assert planets.REFINE_TRIALS <= 100
