"""Where occulta keeps its work: the database, the cache of light curves, the findings, and the analyses of
eclipsing binaries all live under one working directory, chosen with the environment variable OCCULTA_HOME or the
--home option of the command line, and the current directory otherwise."""
import json
import os
from importlib import resources
from pathlib import Path

HOME = Path(os.environ.get('OCCULTA_HOME') or Path.cwd()).resolve()
BIN = HOME / 'binaries'


def config():
    """The search configuration: <home>/config.json when present, the packaged default otherwise."""
    p = HOME / 'config.json'
    if p.exists():
        return json.loads(p.read_text())
    return json.loads(resources.files('occulta').joinpath('default_config.json').read_text())


def init():
    """Create the working directories and a copy of the default configuration that the user can edit."""
    for d in (HOME / 'data' / 'raw', HOME / 'data' / 'cache', HOME / 'findings', HOME / 'logs', BIN / 'eb'):
        d.mkdir(parents=True, exist_ok=True)
    p = HOME / 'config.json'
    if not p.exists():
        p.write_text(resources.files('occulta').joinpath('default_config.json').read_text())
    return HOME
