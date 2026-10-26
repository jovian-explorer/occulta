"""Offline checks of the Gaia checks, the cross-sector candidate rules, the common-mode flag, the pixel model, and the
light-travel-time fit. Archives are replaced by mocked answers and synthetic data of known content."""
import json
import re
import subprocess
import sys
import types

import numpy as np
import pytest

from test_cuts import SUN, _vet


# ---------------------------------------------------------------- Gaia checks
def _fake_gaia(monkeypatch, rows, eb_ids=()):
    """Replace the Gaia archive: `rows` are (source_id, ra, dec, G, ruwe, non_single_star, variable flag, sep arcsec)."""
    from astropy.table import Table
    tab = Table(rows=rows, names=('source_id', 'ra', 'dec', 'phot_g_mean_mag', 'ruwe', 'non_single_star',
                                  'phot_variable_flag', 'sep'))

    class Job:
        def __init__(self, r): self.r = r
        def get_results(self): return self.r

    class Archive:
        @staticmethod
        def launch_job(q):
            if 'vari_eclipsing_binary' in q:
                sid = int(q.split('source_id=')[1])
                return Job(Table({'source_id': np.array([sid] if sid in eb_ids else [], dtype=np.int64)}))
            return Job(tab.copy())

    mod = types.ModuleType('astroquery.gaia'); mod.Gaia = Archive
    monkeypatch.setitem(sys.modules, 'astroquery.gaia', mod)


def test_gaia_checks_report_ruwe_binarity_and_the_blends_that_could_cause_the_dip(monkeypatch):
    """A neighbour 2 mag fainter holds 13.7 % of the pair's light, so that it can cause a 5 % dip (it would need an
    eclipse 37 % deep) but not a 20 % one (it would need 146 %); the target itself is never listed."""
    from occulta import planets
    _fake_gaia(monkeypatch, [(1, 10.0, 20.0, 10.0, 2.0, 1, 'NOT_AVAILABLE', 0.5),
                             (2, 10.0, 20.0, 12.0, 1.0, 0, 'NOT_AVAILABLE', 8.0)], eb_ids=(1,))
    r = planets.gaia_checks(10.0, 20.0, '1', 0.05)
    assert r['ruwe'] == 2.0 and r['non_single_star'] == 1 and r['gaia_eb'] is True
    assert [b['sep_arcsec'] for b in r['possible_blends']] == [8.0]
    assert abs(r['possible_blends'][0]['eclipse_needed'] - 0.365) < 0.01
    assert planets.gaia_checks(10.0, 20.0, '1', 0.20)['possible_blends'] == []
    _fake_gaia(monkeypatch, [(1, 10.0, 20.0, 10.0, 1.0, 0, 'NOT_AVAILABLE', 0.5)])
    assert planets.gaia_checks(10.0, 20.0, '1', 0.05)['gaia_eb'] is False


def test_gaia_checks_use_the_tic_gaia_identifier_not_the_nearest_source(monkeypatch):
    """A high-proper-motion star can lie farther from the catalogue position than a faint neighbour; the identifier
    that the TIC gives decides which source is the target."""
    from occulta import planets
    _fake_gaia(monkeypatch, [(5, 10.0, 20.0, 16.0, 1.0, 0, 'NOT_AVAILABLE', 0.5),
                             (6, 10.0, 20.0, 10.0, 2.5, 0, 'NOT_AVAILABLE', 5.0)], eb_ids=(6,))
    r = planets.gaia_checks(10.0, 20.0, '6', 0.01)
    assert r['ruwe'] == 2.5 and r['gaia_eb'] is True and all(b['sep_arcsec'] != 5.0 for b in r['possible_blends'])
    assert planets.gaia_checks(10.0, 20.0, 'None', 0.01)['ruwe'] == 1.0          # no identifier: the nearest source


def test_vetting_rejects_a_gaia_eclipsing_binary_and_a_high_ruwe_but_not_ruwe_1_4(tmp_path, monkeypatch):
    status, _, notes = _vet(tmp_path, monkeypatch, star=SUN, gaia={'gaia_eb': True})
    assert status == 'rejected' and 'Gaia DR3 eclipsing binary' in notes, notes
    status, _, notes = _vet(tmp_path, monkeypatch, star=SUN, gaia={'ruwe': 1.41})
    assert status == 'rejected' and 'RUWE 1.41' in notes, notes
    status, _, notes = _vet(tmp_path, monkeypatch, star=SUN, gaia={'ruwe': 1.4})
    assert status == 'ready', notes


def test_vetting_rejects_a_signal_shared_with_another_candidate_within_the_tolerance(tmp_path, monkeypatch):
    """Another candidate within 0.3 % in period and 0.3 d in epoch points to a systematic of the spacecraft."""
    status, _, notes = _vet(tmp_path, monkeypatch, other=(7.13 * 1.002, 1003.4))
    assert status == 'rejected' and 'sector-synchronised' in notes, notes
    status, _, notes = _vet(tmp_path, monkeypatch, other=(7.13 * 1.004, 1003.2))          # period 0.4 % away
    assert 'sector-synchronised' not in notes, notes
    status, _, notes = _vet(tmp_path, monkeypatch, other=(7.13, 1003.6))                  # epoch 0.4 d away
    assert 'sector-synchronised' not in notes, notes


# ---------------------------------------------------------------- candidates across sectors and units
def _promoted(tmp_path, evs, inv=3.0, known=None):
    """Promote star 1 with events (sector, t0, duration h, depth ppm, S/N); returns its candidate rows."""
    from occulta import survey
    con = survey.db(tmp_path / f'h{len(list(tmp_path.glob("h*")))}.db')
    for sector in {e[0] for e in evs}:
        con.execute('INSERT INTO star_sector(tic, sector, inv_snr) VALUES (1,?,?)', (sector, inv))
    for sector, t0, dur, depth, snr in evs:
        con.execute("INSERT INTO event VALUES (1,?,1,1,?,?,?,?,20,0.5,'')", (sector, t0, dur, depth, snr))
    if known:
        con.execute('INSERT INTO known(tic, src, period, epoch, dur_h, depth_ppm) VALUES (1,?,?,?,3,2500)', known)
    con.commit()
    survey.promote(con)
    return con.execute('SELECT kind, n_events FROM candidate').fetchall()


@pytest.mark.parametrize('second, inv, kept', [
    ((2, 1040.0, 3.0, 3000, 12.0), 3.0, True),
    ((2, 1040.0, 3.0, 3000, 8.9), 3.0, False),          # S/N below 9
    ((2, 1040.0, 3.0, 3000, 12.0), 11.0, False),        # S/N below 1.2 x the inverted-curve level (13.2)
    ((2, 1040.0, 3.0, 4500, 12.0), 3.0, True),          # depth ratio 1.5 (limit 1.6)
    ((2, 1040.0, 3.0, 6000, 12.0), 3.0, False),         # depth ratio 2
    ((2, 1040.0, 5.4, 3000, 12.0), 3.0, True),          # duration ratio 1.8 (limit 2)
    ((2, 1040.0, 7.5, 3000, 12.0), 3.0, False),         # duration ratio 2.5
])
def test_stitching_needs_strong_similar_dips_above_the_false_alarm_level(tmp_path, second, inv, kept):
    got = _promoted(tmp_path, [(1, 1003.2, 3.0, 3000, 13.5), second], inv=inv)
    assert got == ([('multi', 2)] if kept else [])


def test_stitching_skips_dips_of_a_known_ephemeris_and_catalogued_binaries(tmp_path):
    evs = [(1, 1003.2, 3.0, 3000, 13.5), (2, 1040.0, 3.0, 3000, 12.0)]
    assert _promoted(tmp_path, evs, known=('TOI', 18.4, 1003.3)) == []          # both dips are transits of a TOI
    assert _promoted(tmp_path, evs, known=('EB catalogue Prsa', None, None)) == []
    assert _promoted(tmp_path, evs, known=('TOI', 5.0, 900.5)) == [('multi', 2)]  # a TOI whose ephemeris does not fit


def test_bls_candidates_need_snr_sde_three_transits_and_a_consistent_odd_even(tmp_path):
    """The limits are S/N >= 10, SDE >= 9, at least 3 transits, odd/even below 3 sigma."""
    from occulta import survey
    con = survey.db(tmp_path / 'bls.db')
    rows = [(1, 12.0, 11.0, 5, 0.5),      # accepted
            (2, 10.0, 9.0, 3, 2.9),       # exactly at the limits: accepted
            (3, 9.9, 11.0, 5, 0.5),       # S/N
            (4, 12.0, 8.9, 5, 0.5),       # SDE
            (5, 12.0, 11.0, 2, 0.5),      # transits
            (6, 12.0, 11.0, 5, 3.0),      # odd/even
            (7, 12.0, 11.0, 5, 0.5)]      # as 1, but a known planet transits there
    for tic, snr, sde, ntr, oe in rows:
        con.execute('INSERT INTO bls VALUES (?,1,2.5,1001.0,2.0,1500,?,?,?,?,0)', (tic, snr, sde, ntr, oe))
    con.execute("INSERT INTO known(tic, src, period, epoch, dur_h, depth_ppm) VALUES (7,'TOI',2.5,1001.0,2,1500)")
    con.commit()
    survey.promote(con)
    assert [r[0] for r in con.execute("SELECT tic FROM candidate WHERE kind='bls' ORDER BY tic")] == [1, 2]


# ---------------------------------------------------------------- common-mode systematics
def test_common_mode_flag_needs_five_events_within_a_fifth_of_a_day_on_one_ccd(tmp_path):
    from occulta import survey
    con = survey.db(tmp_path / 'cm.db')

    def ins(cam, ccd, t0, tic, fl=''):
        con.execute("INSERT INTO event VALUES (?,5,?,?,?,3,2500,12,20,0.5,?)", (tic, cam, ccd, t0, fl))
    for k in range(5): ins(1, 1, 1000 + 0.02 * k, 100 + k, 'edge' if k == 0 else '')     # five together: all flagged
    for k in range(4): ins(1, 2, 1000 + 0.02 * k, 200 + k)                                # four: none
    for k in range(5): ins(1, 3, 1000 + 0.15 * k, 300 + k)                                # five, spread over 0.6 d: none
    for k in range(5): ins(2, 1, 1000 + 0.02 * k, 400 + k)                                # another camera: its own count
    con.commit()
    survey.mark_systematics(con, 5)
    flagged = {r[0]: r[1] for r in con.execute("SELECT tic, flags FROM event WHERE flags LIKE '%common%'")}
    assert set(flagged) == set(range(100, 105)) | set(range(400, 405))
    assert flagged[100] == 'edge,common' and flagged[101] == 'common'


# ---------------------------------------------------------------- pixel model
def _tesscut(path, src_pos, eclipsed, P=2.0, T0=1001.0, dur=3 / 24, sig=1.0, seed=3):
    """A 15 x 15 TESScut-like file: two stars of amplitude 300 (target, at pixel 7,7) and 800 (at 10,7), a background of
    50, and a 25 % eclipse of the star `eclipsed` (0 or 1) every P days."""
    from astropy.io import fits
    rng = np.random.default_rng(seed)
    t = np.arange(1000.0, 1012.0, 10 / 1440)
    yy, xx = np.mgrid[:15, :15]
    psf = lambda x, y: np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / sig ** 2)
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    inecl = np.abs(ph * P) < dur / 2
    cube = np.empty((t.size, 15, 15))
    for i in range(t.size):
        a = [300.0, 800.0]
        if inecl[i]:
            a[eclipsed] *= 0.75
        cube[i] = 50 + sum(ai * psf(x, y) for ai, (x, y) in zip(a, src_pos)) + rng.normal(0, 1.0, (15, 15))
    wcs = fits.Header({'NAXIS': 2, 'NAXIS1': 15, 'NAXIS2': 15, 'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CRPIX1': 8.0,
                       'CRPIX2': 8.0, 'CRVAL1': 100.0, 'CRVAL2': 0.0, 'CDELT1': -21 / 3600, 'CDELT2': 21 / 3600})
    cols = [fits.Column(name='TIME', format='D', array=t), fits.Column(name='FLUX', format='225E', dim='(15,15)', array=cube),
            fits.Column(name='QUALITY', format='J', array=np.zeros(t.size, int))]
    h0 = fits.PrimaryHDU(); h0.header['SECTOR'] = 9
    fits.HDUList([h0, fits.BinTableHDU.from_columns(cols), fits.ImageHDU(np.zeros((15, 15), np.float32), header=wcs)]
                 ).writeto(path)


def test_sources_closer_than_8_arcsec_are_merged_into_one_blend():
    from astropy.table import Table
    from occulta.binaries import pixels
    src = Table(rows=[(1, 100.0, 0.0, 10.0), (2, 100.0 + 5 / 3600, 0.0, 10.0), (3, 100.0 + 30 / 3600, 0.0, 12.0)],
                names=('source_id', 'ra', 'dec', 'g'))
    m = pixels.merge_close(src, 8.0)
    assert list(m['n_members']) == [2, 1] and list(m['source_id']) == [1, 3]
    assert abs(m['g'][0] - (10.0 - 2.5 * np.log10(2))) < 1e-9 and abs(m['ra'][0] - (100.0 + 2.5 / 3600)) < 1e-9


@pytest.mark.parametrize('eclipsed, on_target', [(0, True), (1, False)])
def test_pixel_model_finds_which_star_eclipses(tmp_path, monkeypatch, eclipsed, on_target):
    """The difference image of a 25 % eclipse is fitted by the point-spread function of one Gaia source: the faint
    target lies 3 pixels (63 arcsec) from a brighter neighbour, whose light must not be mistaken for the eclipse."""
    from astropy.table import Table
    from occulta.binaries import pixels
    work = tmp_path / 'data' / 'cache' / 'tesscut15_9'; work.mkdir(parents=True)
    _tesscut(work / 'cut.fits', [(7, 7), (10, 7)], eclipsed)
    dra = 3 * 21 / 3600
    gaia = Table(rows=[(1, 100.0, 0.0, 13.0, 1.0, 1.0), (2, 100.0 - dra, 0.0, 12.0, 1.0, 1.0)],
                 names=('source_id', 'ra', 'dec', 'g', 'parallax', 'ruwe'))
    monkeypatch.setattr(pixels, 'ROOT', tmp_path); monkeypatch.setattr(pixels, 'OUT', tmp_path / 'out')
    monkeypatch.setattr(pixels, 'gaia_sources', lambda ra, dec, **k: gaia)
    pixels.main(['9', '2.0', '1001.0', '3.0', '100.0', '0.0'])
    v = json.loads((tmp_path / 'out' / '9.json').read_text())
    assert v['target_source'] == 1 and v['best_source'] == (1 if on_target else 2) and v['on_target'] is on_target, v
    assert v['n_sectors'] == 1 and v['median_chi2_ratio_next'] > 3, v
    d = v['sectors'][0]['depth_if_alone']
    assert abs(d[str(v['best_source'])] - 0.25) < 0.05, d            # the best source loses 25 % of its light


# ---------------------------------------------------------------- light-travel time of a Gaia orbit
def test_lighttime_fit_recovers_a_gaia_shaped_term_and_omega_from_thiele_innes(tmp_path):
    (tmp_path / 'binaries' / 'eb').mkdir(parents=True)
    P, T0, Pout, e, w_true, K = 2.5, 1000.0, 300.0, 0.4, np.radians(40.0), 3 / 1440
    Tp = 1100.0
    # Thiele-Innes constants of an orbit with omega 40, Omega 70, i 60 degrees (Halbwachs et al. 2023)
    om, Om, inc = w_true, np.radians(70.0), np.radians(60.0)
    A = np.cos(om) * np.cos(Om) - np.sin(om) * np.sin(Om) * np.cos(inc)
    B = np.cos(om) * np.sin(Om) + np.sin(om) * np.cos(Om) * np.cos(inc)
    F = -np.sin(om) * np.cos(Om) - np.cos(om) * np.sin(Om) * np.cos(inc)
    G = -np.sin(om) * np.sin(Om) + np.cos(om) * np.cos(Om) * np.cos(inc)
    E = np.arange(0, 700, 2)
    t0 = T0 + P * E
    M = 2 * np.pi * (t0 - Tp) / Pout
    Ea = M.copy()
    for _ in range(50):
        Ea -= (Ea - e * np.sin(Ea) - M) / (1 - e * np.cos(Ea))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Ea / 2), np.sqrt(1 - e) * np.cos(Ea / 2))
    t = t0 + K * ((1 - e ** 2) / (1 + e * np.cos(nu)) * np.sin(nu + w_true) + e * np.sin(w_true))     # Irwin (1952)
    rows = ['unit,type,cycle,t_btjd,t_err_d,depth,depth_err,sigma_d,rms']
    rows += [f'1001,P,{k},{float(x)!r},{30 / 86400!r},0.1,0.001,0.05,0.001' for k, x in zip(E, t)]
    (tmp_path / 'binaries' / 'eb' / '8_eclipses.csv').write_text('\n'.join(rows) + '\n')
    args = ['8', repr(P), repr(T0), repr(Pout), repr(e), repr(float(Tp - 389.0)), *(repr(float(x)) for x in (A, B, F, G))]
    out = subprocess.run([sys.executable, '-m', 'occulta.cli', '--home', str(tmp_path), 'binaries', 'lighttime', *args],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-500:]
    w = float(re.search(r'omega (-?[\d.]+) deg', out.stdout).group(1))
    assert abs(((w - 40.0) + 90) % 180 - 90) < 0.1, out.stdout                    # omega is known up to 180 degrees
    chi_lin = float(re.search(r'linear: chi2 ([\d.]+)', out.stdout).group(1))
    chi, k = map(float, re.search(r'Gaia-shape LTTE: chi2 ([\d.]+), K ([+-][\d.]+) min', out.stdout).groups())
    assert chi < 1e-3 * chi_lin and abs(abs(k) - 3.0) < 0.05, out.stdout
