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
