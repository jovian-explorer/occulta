"""Offline checks: the package imports, the command line answers, and a working directory can be created."""
import importlib
import pkgutil
import subprocess
import sys

import numpy as np
import pytest

import occulta


def test_every_module_imports():
    for m in pkgutil.walk_packages(occulta.__path__, 'occulta.'):
        importlib.import_module(m.name)


def test_version_and_help():
    out = subprocess.run([sys.executable, '-m', 'occulta.cli', '--version'], capture_output=True, text=True)
    assert out.stdout.strip() == occulta.__version__
    out = subprocess.run([sys.executable, '-m', 'occulta.cli', 'binaries', '--help'], capture_output=True, text=True)
    assert 'pixels' in out.stdout and 'lighttime' in out.stdout


def test_init_creates_working_directory(tmp_path):
    subprocess.run([sys.executable, '-m', 'occulta.cli', '--home', str(tmp_path), 'init'], check=True)
    assert (tmp_path / 'config.json').exists()
    assert (tmp_path / 'binaries' / 'eb').is_dir()


def test_status_and_map_of_a_new_working_directory(tmp_path):
    _small_footprint(tmp_path)         # no 7-minute tess-point run when it is installed
    for sub in ('status', 'map'):
        out = subprocess.run([sys.executable, '-W', 'ignore', '-m', 'occulta.cli', '--home', str(tmp_path), 'survey', sub],
                             capture_output=True, text=True)
        assert out.returncode == 0, out.stderr[-500:]
    assert 'stars 0' in (tmp_path / 'logs' / 'survey.log').read_text()
    assert (tmp_path / 'findings' / 'progress' / 'progress_map.png').exists()


def test_survey_stops_cleanly_when_the_known_lists_cannot_be_downloaded(tmp_path, monkeypatch):
    import json
    import requests
    from occulta import survey
    for name in ('ROOT', 'RAW', 'CACHE', 'FIND', 'LOGS'):
        monkeypatch.setattr(survey, name, tmp_path / name.lower())
    monkeypatch.setattr(survey, 'DB', tmp_path / 'hunt.db')
    monkeypatch.setattr(survey, 'LOCK', tmp_path / 'run.lock')
    def fail(con):
        raise requests.exceptions.SSLError('archive unreachable')
    monkeypatch.setattr(survey, 'refresh_known', fail)
    survey.run(1)
    assert 'could not be downloaded' in (tmp_path / 'logs' / 'survey.log').read_text()


def _small_footprint(home):
    """A footprint of 300 positions with a known sector count, in place of the 7-minute tess-point run."""
    rng = np.random.default_rng(3)
    (home / 'data').mkdir(parents=True, exist_ok=True)
    np.savez(home / 'data' / 'footprint.npz', ra=rng.uniform(0, 360, 300), dec=np.degrees(np.arcsin(rng.uniform(-1, 1, 300))),
             nsec=rng.integers(0, 12, 300))


def test_map_is_shaded_by_sectors_from_a_cached_footprint(tmp_path):
    _small_footprint(tmp_path)
    out = subprocess.run([sys.executable, '-W', 'ignore', '-m', 'occulta.cli', '--home', str(tmp_path), 'survey', 'map'],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-500:]
    assert 'TESS covers' in out.stdout and 'QLP 107.6 M' in out.stdout
    assert (tmp_path / 'findings' / 'progress' / 'progress_map.png').stat().st_size > 10000


def test_map_without_tess_point_is_drawn_without_shading(tmp_path, monkeypatch, capsys):
    from occulta import progress
    monkeypatch.setitem(sys.modules, 'tess_stars2px', None)       # import fails as when tess-point is not installed
    monkeypatch.setattr(progress, 'ROOT', tmp_path)
    assert progress.footprint() is None
    progress.main()
    assert 'tess-point is not installed' in capsys.readouterr().out
    assert (tmp_path / 'findings' / 'progress' / 'progress_map.png').exists()
    assert not (tmp_path / 'data' / 'footprint.npz').exists()


def test_footprint_from_tess_point_is_cached_in_the_home(tmp_path, monkeypatch):
    t2p = pytest.importorskip('tess_stars2px')
    from occulta import progress
    monkeypatch.setattr(progress, 'ROOT', tmp_path)
    ra, dec, nsec = progress.footprint(200)
    assert ra.size == dec.size == nsec.size == 200 and nsec.max() > 0 and nsec.max() <= progress.NS
    assert (tmp_path / 'data' / 'footprint.npz').exists()
    monkeypatch.setattr(t2p, 'tess_stars2px_function_entry', None)    # a second call must read the cache
    assert (progress.footprint(200)[2] == nsec).all()


def test_refresh_totals_counts_at_mast_keeps_the_cache_and_falls_back_on_failure(tmp_path, monkeypatch, capsys):
    """`--refresh-totals` replaces the dated numbers of the products that MAST counts, writes them to the home, and
    keeps the dated number of a product whose count fails; without the option the cache is read and nothing is queried."""
    import json
    import requests
    from occulta import progress
    monkeypatch.setattr(progress, 'ROOT', tmp_path)
    assert progress.mast_totals() == progress.MAST_TOTALS                       # no cache: the dated numbers

    def post(url, data, timeout):
        flt = json.loads(data['request'])['params']['filters']
        if any(f['values'] == ['K2SFF'] for f in flt):
            raise requests.exceptions.ReadTimeout('slow')
        resp = type('R', (), {'raise_for_status': lambda s: None, 'json': lambda s: {'data': [{'Column1': 2_345_678}]}})()
        return resp
    monkeypatch.setattr(requests, 'post', post)
    tot = progress.mast_totals(refresh=True)
    assert tot['SPOC 2-min'][0] == 2.346 and tot['Kepler'][0] == 2.346 and tot['SPOC 2-min'][1] > '2026-10-25'
    assert tot['K2SFF'] == progress.MAST_TOTALS['K2SFF'] and tot['QLP'] == progress.MAST_TOTALS['QLP']
    assert 'K2SFF: the MAST count failed' in capsys.readouterr().out
    monkeypatch.setattr(requests, 'post', None)                                 # a plain call reads the file
    assert progress.mast_totals()['Kepler'][0] == 2.346


@pytest.mark.parametrize('tool', sorted(__import__('occulta.cli', fromlist=['x']).BINARY_TOOLS))
def test_every_binaries_tool_prints_its_usage_for_help(tool, tmp_path):
    """`occulta binaries <tool> --help` printed a ValueError for tools that read their arguments by position."""
    out = subprocess.run([sys.executable, '-W', 'ignore', '-m', 'occulta.cli', '--home', str(tmp_path), 'binaries', tool, '--help'],
                         capture_output=True, text=True)
    assert out.returncode == 0 and 'Usage' in out.stdout and 'Traceback' not in out.stderr, (out.stdout, out.stderr[-300:])


def _cli(*args, cwd=None, env=None):
    import os
    e = {k: v for k, v in os.environ.items() if k != 'OCCULTA_HOME'}
    e.update(env or {})
    return subprocess.run([sys.executable, '-W', 'ignore', '-m', 'occulta.cli', *args], capture_output=True, text=True,
                          cwd=cwd, env=e)


def test_no_working_directory_gives_a_hint_not_a_new_folder(tmp_path):
    out = _cli('survey', 'status', cwd=tmp_path)
    assert out.returncode != 0 and 'occulta init' in out.stderr and '--home' in out.stderr
    assert list(tmp_path.iterdir()) == []          # nothing was created behind the user's back
    assert _cli('init', cwd=tmp_path).returncode == 0
    assert _cli('survey', 'status', cwd=tmp_path).returncode == 0     # after init the current directory is a home


def test_a_name_instead_of_a_tic_number_is_explained(tmp_path):
    for cmd in (('reproduce', 'WASP-18'), ('planets', 'vet', 'TOI-185')):
        out = _cli('--home', str(tmp_path), *cmd)
        assert out.returncode != 0 and 'is not a TIC number' in out.stderr and 'Traceback' not in out.stderr


def test_unknown_tic_and_no_network_are_explained(tmp_path, monkeypatch):
    import requests
    from occulta import planets, reproduce
    monkeypatch.setenv('OCCULTA_HOME', str(tmp_path))
    monkeypatch.setattr(reproduce.paths, 'init', lambda: tmp_path)
    monkeypatch.setattr(reproduce, 'ROOT', tmp_path)
    monkeypatch.setattr(planets, 'fetch_all', lambda tic, d: [])
    with pytest.raises(SystemExit) as e:
        reproduce.main(['99999999999'])
    assert 'No TESS, K2 or Kepler light curves of TIC 99999999999' in str(e.value)

    def offline(tic, d):
        raise requests.exceptions.ConnectionError('no route')
    monkeypatch.setattr(planets, 'fetch_all', offline)
    with pytest.raises(SystemExit) as e:
        reproduce.main(['410153553'])
    assert 'Could not reach the MAST archive' in str(e.value) and 'internet connection' in str(e.value)
