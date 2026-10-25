"""Offline checks of the science code on synthetic light curves: no download is made, the archive calls are replaced
by light curves with a known signal."""
import subprocess
import sys

import numpy as np


def _transits(t, P, t0, depth, dur_d):
    ph = ((t - t0) / P + 0.5) % 1 - 0.5
    return np.where(np.abs(ph * P) < dur_d / 2, -depth, 0.0)


def _vet_synthetic(tmp_path, monkeypatch, starts, P, t0, depth, dur, cadence_min=10, noise=400e-6):
    """Vet star 1 by hand on synthetic sectors of 26 d that begin at `starts`; returns (status, period, notes)."""
    from occulta import planets, survey
    rng = np.random.default_rng(1)
    lcs = {}
    for k, start in enumerate(starts):
        t = np.arange(start, start + 26, cadence_min / 1440)
        lcs[str(k)] = (t, _transits(t, P, t0, depth, dur) + rng.normal(0, noise, t.size))
    monkeypatch.setattr(planets, 'fetch_all', lambda sid, work: [('QLP', k + 1, str(k)) for k in range(len(starts))])
    monkeypatch.setattr(planets, 'series', lambda author, path: lcs[path])
    monkeypatch.setattr(planets, 'get_star', lambda con, sid: None)
    con = survey.db(tmp_path / 'hunt.db')
    con.execute("INSERT INTO candidate(tic, kind, status) VALUES (1, 'manual', 'new')")
    planets.vet(con, 1, tmp_path / 'findings')
    return con.execute('SELECT status, period, notes FROM candidate WHERE tic=1').fetchone()


def test_vet_by_hand_finds_the_period(tmp_path, monkeypatch):
    """`occulta planets vet TIC` on a star that the survey has not promoted (no candidate depth) must still pair its
    dips and recover the period."""
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, (1000.0, 1030.0), 7.13, 1003.2, 2500e-6, 3 / 24)
    assert period is not None and abs(period - 7.13) < 0.01, (status, period, notes)


def test_vet_of_dips_spread_over_years_finishes_in_seconds(tmp_path, monkeypatch):
    """Dips in 30 sectors over 2000 d once made the pair-by-pair period search run for hours (TIC 348997611); it must
    now finish quickly and still find the period."""
    import time
    start = time.monotonic()
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, 1000.0 + 68.0 * np.arange(30), 9.71, 1001.3,
                                           5000e-6, 4 / 24, cadence_min=30, noise=600e-6)
    assert time.monotonic() - start < 600     # hours before the fix; under a minute on an idle PC
    assert period is not None and abs(period - 9.71) < 0.01, (status, period, notes)


def test_vet_that_runs_too_long_is_marked_vet_error(tmp_path, monkeypatch):
    from occulta import planets
    monkeypatch.setattr(planets, 'VET_LIMIT_S', -1)
    status, period, notes = _vet_synthetic(tmp_path, monkeypatch, (1000.0, 1030.0), 7.13, 1003.2, 2500e-6, 3 / 24)
    assert status == 'vet_error' and 'period search stopped' in notes


def test_binary_centroid_shift_uses_the_pixel_of_each_mission(tmp_path, monkeypatch):
    """The in-eclipse centroid shift is converted to arcsec with the TESS pixel (21") for TESS sectors and the Kepler
    pixel (3.98") for Kepler quarters; the period and the secondary phase of an eclipsing binary are recovered."""
    from occulta.binaries import analysis
    rng = np.random.default_rng(2)
    P, t0, w = 3.4, 1001.0, 4 / 24
    parts = []
    for u, t in ((5, np.arange(1000, 1027, 10 / 1440)), (3005, np.arange(1100, 1190, 30 / 1440))):
        ecl = _transits(t, P, t0, 0.2, w) + _transits(t, P, t0 + P / 2, 0.1, w)
        f = 1000 * (1 + ecl + rng.normal(0, 1e-3, t.size))
        cx = 50 + 0.01 * (ecl < -0.15) + rng.normal(0, 1e-4, t.size)     # 0.01 pixel shift in primary eclipse
        parts.append((t, f, np.full(t.size, u), cx, np.full(t.size, 50.0) + rng.normal(0, 1e-4, t.size)))
    data = tuple(np.concatenate([p[i] for p in parts]) for i in range(5))
    monkeypatch.setattr(analysis, 'load', lambda tic: data)
    monkeypatch.setattr(analysis, 'OUT', tmp_path)
    r = analysis.analyse(1, P, t0, 4.0, {})
    assert abs(r['period'] - P) < 1e-3 and not r['doubled']
    assert abs(r['secondary_phase_scan'] - 0.5) < 0.01
    shift = {c['unit']: c['shift_arcsec'] for c in r['centroid']}
    assert abs(shift[5] - 0.21) < 0.03 and abs(shift[3005] - 0.04) < 0.01, shift


def test_lighttime_scan_finds_an_injected_outer_orbit(tmp_path):
    """`occulta binaries lighttime-scan` reads the fitted eclipse times and recovers a sinusoidal light-travel-time
    term."""
    (tmp_path / 'binaries' / 'eb').mkdir(parents=True)
    P, T0, Pout, amp = 2.5, 1000.0, 400.0, 3 / 1440
    E = np.arange(0, 600, 3)
    t = T0 + P * E + amp * np.sin(2 * np.pi * (T0 + P * E) / Pout)
    rows = ['unit,type,cycle,t_btjd,t_err_d,depth,depth_err,sigma_d,rms']
    rows += [f'1001,P,{e},{float(x)!r},{30 / 86400!r},0.1,0.001,0.05,0.001' for e, x in zip(E, t)]
    (tmp_path / 'binaries' / 'eb' / '7_eclipses.csv').write_text('\n'.join(rows) + '\n')
    out = subprocess.run([sys.executable, '-m', 'occulta.cli', '--home', str(tmp_path), 'binaries', 'lighttime-scan',
                          '7', repr(P), repr(T0)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-500:]
    best = float(out.stdout.splitlines()[1].split()[1])
    assert abs(best / Pout - 1) < 0.02, out.stdout
