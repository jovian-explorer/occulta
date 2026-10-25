"""Re-analyse stars from the cached light curves at their CURRENT adopted period and epoch (eb/<TIC>.json), so that
period corrections made earlier (divided or doubled periods) are kept. Usage: occulta binaries rerun [TIC ...]"""
import json, sys
from pathlib import Path

from .. import paths
from . import analysis as eb_vet
ROOT = paths.HOME


def redo(tic):
    s = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())
    old = json.loads((eb_vet.OUT / f'{tic}.json').read_text())
    w = min(old.get('width_h') or s['t14_h'], old['period'] * 24 / 4)
    r = eb_vet.analyse(tic, old['period'], old['t0_btjd'], w, s.get('gaia', {}))
    return old, r


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    for tic in [int(x) for x in argv]:
        old, r = redo(tic)
        print(tic, 'P %.7f -> %.7f' % (old['period'], r['period']), 'T0 %.5f -> %.5f (%.1f min)' % (
            old['t0_btjd'], r['t0_btjd'], ((r['t0_btjd'] - old['t0_btjd'] + r['period'] / 2) % r['period'] - r['period'] / 2) * 1440),
            'sec phase', (r.get('secondary') or {}).get('phase'), flush=True)


if __name__ == '__main__':
    main()
