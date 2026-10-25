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
