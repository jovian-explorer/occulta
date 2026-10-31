"""Offline checks of the single-transit and dipper search (singles.py): the period from the duration, the coverage
limit, the shape test of a transit against a dusty tail, and the classification of the dips in a survey database."""
import numpy as np
import pytest

from occulta import singles as S
from occulta.planets import trapezoid_ss


def test_period_from_duration_of_the_earth():
    # the Earth crosses the Sun in 13 h; the small-angle form gives about a year
    assert S.period_from_duration(13.0, 0.0092, 1.0, 1.0) == pytest.approx(365.25, rel=0.05)
    assert S.period_from_duration(13.0, 0.0092, 1.0, 1.0, 0.8) > 3 * S.period_from_duration(13.0, 0.0092, 1.0, 1.0)


def test_minimum_period_from_coverage():
    # a dip at 1010 in a window 1000-1027: a repeat after it must fall beyond 1027
    assert S.min_period(1010.0, 4.0, [(1000.0, 1027.0)]) == pytest.approx(17.0 - 4 / 48, rel=2e-3)
    # a second sector 50 d later closes the periods that would put a repeat inside it
    assert S.min_period(1010.0, 4.0, [(1000.0, 1027.0), (1060.0, 1087.0)]) > 30


@pytest.mark.parametrize('seed', range(3))
def test_shape_separates_transit_and_tail(seed):
    rng = np.random.default_rng(seed)
    t = np.arange(1000.0, 1027.0, 2 / 1440)
    tr = 1 + trapezoid_ss(t, 1010.0, 3000e-6, 5 / 24, 0.5 / 24, 2 / 1440) + rng.normal(0, 1000e-6, t.size)
    tl = 1 + S.tail(t, 1010.0, 0.05, 3000e-6, 0.15) + rng.normal(0, 1000e-6, t.size)
    assert S.search_lc(t, tr, [{'t0': 1010.0, 'dur_h': 5.0}])[0][3] == 'transit'
    assert S.search_lc(t, tl, [{'t0': 1010.1, 'dur_h': 5.0}])[0][3] == 'dipper'


def test_classify_dips_in_a_survey_database(tmp_path):
    from occulta import survey
    con = survey.db(str(tmp_path / 'hunt.db'))
    con.execute('INSERT INTO tic (tic, rad, mass) VALUES (1, 1.0, 1.0)')
    ev = [(1, 10, 1010.0, 5.0, 3000, 12.0, ''),            # a lone clean dip: single
          (2, 10, 1005.0, 3.0, 2000, 11.0, ''),            # two dips of very different depth: aperiodic
          (2, 10, 1020.0, 9.0, 9000, 15.0, ''),
          (3, 10, 1008.0, 4.0, 3000, 14.0, ''),            # already a periodic candidate
          (4, 10, 1008.0, 4.0, 3000, 14.0, 'vshape'),      # flagged
          (5, 10, 1008.0, 4.0, 3000, 8.0, ''),             # too weak for a lone dip
          (6, 10, 1008.0, 4.0, 3000, 12.0, '')]            # the inverted search found as strong a brightening
    for tic, s, t0, d, dep, snr, fl in ev:
        con.execute('INSERT INTO event VALUES (?,?,1,1,?,?,?,?,10,1.0,?)', (tic, s, t0, d, dep, snr, fl))
    for tic in range(1, 7):
        con.execute('INSERT INTO star_sector VALUES (?,10,1,1,3000,1000.0,1027.0,1000,300,?)', (tic, 13.0 if tic == 6 else 5.0))
    con.execute("INSERT INTO candidate (tic, kind) VALUES (3, 'planet')")
    rows = S.classify(con, inv_ratio=1.2)
    kinds = {(r['tic'], r['kind']) for r in rows}
    assert kinds == {(1, 'single'), (2, 'aperiodic')}
    one = next(r for r in rows if r['tic'] == 1)
    assert one['p_min_d'] == pytest.approx(17.0 - 5 / 48, rel=2e-3)
    assert one['p_dur_lo'] <= one['p_dur_d'] <= one['p_dur_hi']
    S.store(con, rows, tmp_path)
    assert (tmp_path / 'singles' / 'report.md').exists() and (tmp_path / 'singles' / 'singles.csv').exists()
    assert con.execute('SELECT COUNT(*) FROM single').fetchone()[0] == 3
