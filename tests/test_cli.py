"""Offline checks: the package imports, the command line answers, and a working directory can be created."""
import importlib
import pkgutil
import subprocess
import sys

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
