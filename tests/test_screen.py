from occulta import recovered, screen


def row(**k):
    r = dict(tic=1000, kind='multi', period=5.0, t0=100.0, dur_h=3.0, snr=40.0, n_events=3, sectors='10,11,12',
             objtype='planet', notes='class: planet (Rp 1.0 R_J from R_star 1.0 R_sun)')
    r.update(k)
    return r


def events(n=3, p=5.0, t0=100.0, depth=5000.0, flags='', sector=10):
    return [(sector + i, t0 + i * p, depth, flags) for i in range(n)]


def codes(r, ev, star=None, sib=0):
    return {c for c, _ in screen.reject(r, ev, star or {}, sib)}


def test_clean_candidate_is_kept():
    assert codes(row(), events()) == set()


def test_few_transits_in_one_sector():
    assert 'FEW' in codes(row(n_events=2, snr=20.0, sectors='10,10'), events(2, sector=10))
    assert 'FEW' not in codes(row(n_events=2, snr=60.0, sectors='10,10'), events(2))
    assert 'FEW' not in codes(row(n_events=2, snr=20.0, sectors='10,11'), events(2))


def test_provisional_sectors_only():
    assert 'PROV' in codes(row(sectors='100,101,101'), events(3, sector=100))
    assert 'PROV' not in codes(row(sectors='100,101,40'), events(3, sector=100))


def test_edge_and_orbit_alias():
    assert 'GAP' in codes(row(), events(flags='edge'))
    r = row(period=13.7)
    assert 'GAP' in codes(r, events(3, p=13.7, flags='edge')[:1] + events(3, p=13.7)[1:])
    assert 'GAP' not in codes(r, events(3, p=13.7))


def test_centroid_background_vshape():
    assert 'CENTROID' in codes(row(), events(flags='centroid_x'))
    assert 'BKG' in codes(row(), events(flags='bkg'))
    assert 'VSHAPE' in codes(row(), events(flags='vshape'))
    assert 'VSHAPE' not in codes(row(objtype='eclipsing binary'), events(flags='vshape'))


def test_depth_changes_between_sectors():
    ev = [(10, 100.0, 1000.0, ''), (11, 105.0, 6000.0, ''), (12, 110.0, 1100.0, '')]
    assert 'DEPTH' in codes(row(), ev)


def test_cross_star_and_blend():
    assert 'CROSS' in codes(row(), events(), sib=2)
    assert 'CROSS' not in codes(row(), events(), sib=1)
    assert 'BLEND' in codes(row(), events(), {'contratio': 1.4})
    assert 'BLEND' not in codes(row(), events(), {'contratio': 0.2})


def test_stored_statistics_repeated():
    assert 'ODDEVEN' in codes(row(notes='odd/even 4.2 sigma'), events())
    assert 'ODDEVEN' not in codes(row(objtype='eclipsing binary', notes='odd/even 4.2 sigma'), events())
    assert 'SECOND' in codes(row(notes='no secondary at any phase (best 12.0 sigma, threshold 10.0)'), events())


def test_oversized_planet_is_relabelled():
    big = row(notes='class: planet (Rp 1.8 R_J from R_star 1.0 R_sun)')
    assert screen.relabel(big, {'contratio': 0.5})[0] == 'unclassified'
    assert screen.relabel(big, {'contratio': 0.0}) is None
    assert screen.relabel(row(objtype='eclipsing binary'), {}) is None


def test_siblings_share_period_and_epoch():
    rs = [row(tic=1, period=5.0, t0=100.0), row(tic=2, period=5.005, t0=100.05), row(tic=3, period=5.001, t0=110.04),
          row(tic=4, period=7.0, t0=100.0)]
    assert screen.siblings_of(rs) == {1: 2, 2: 2, 3: 2, 4: 0}


def test_period_relation_and_phase():
    assert recovered.period_relation(5.0, 5.04) == 'P'
    assert recovered.period_relation(10.0, 5.0) == '2P'
    assert recovered.period_relation(5.0, 10.0) == 'P/2'
    assert not recovered.period_relation(5.0, 7.0)
