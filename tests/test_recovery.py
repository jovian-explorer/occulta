"""Recovery tests: the pipeline must find known objects from public data with no prior information. They download
about 200 MB from MAST and take about 15 minutes, so they run only when OCCULTA_ONLINE=1 is set."""
import json
import os
import subprocess
import sys

import pytest

CASES = [  # (TIC, kind, expected period, tolerance, description)
    (402710280, 'planet', 36.3401, 0.01, 'TIC 402710280.01, a published planet candidate'),
    (181149914, 'eb', 10.75981, 0.001, 'CPD-37 4813, an eccentric eclipsing binary found by this search'),
]


@pytest.mark.skipif(os.environ.get('OCCULTA_ONLINE') != '1', reason='needs MAST; set OCCULTA_ONLINE=1')
@pytest.mark.parametrize('tic,kind,period,tol,what', CASES)
def test_recovers_known_object(tmp_path, tic, kind, period, tol, what):
    out = tmp_path / str(tic)
    subprocess.run([sys.executable, '-W', 'ignore', '-m', 'occulta.cli', '--home', str(tmp_path), 'reproduce', str(tic),
                    '--out', str(out), '--mode', 'both'], check=True)
    m = json.loads((out / 'manifest.json').read_text())
    got = (m['planet_vetting'] or {}).get('period_d') if kind == 'planet' else (m['eclipsing_binary'] or {}).get('period')
    assert got is not None and abs(got - period) < tol, f'{what}: expected {period} d, got {got}'
