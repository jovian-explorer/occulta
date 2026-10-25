"""Re-analyse stars from the cached light curves at their CURRENT adopted period and epoch (eb/<TIC>.json), so that
period corrections made earlier (divided or doubled periods) are kept. Without a TIC, every analysed star is redone (after
a change to the analysis); --divide N analyses at 1/N of the adopted period, for a survey period that is N times the true one.
Usage: occulta binaries rerun [--divide N] [TIC ...]"""
import argparse, json, sys

from .. import paths
from . import analysis as eb_vet
ROOT = paths.HOME


def redo(tic, divide=1):
    s = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())
    old = json.loads((eb_vet.OUT / f'{tic}.json').read_text())
    P = old['period'] / divide
    w = min(old.get('width_h') or s['t14_h'], P * 24 / 4)
    r = eb_vet.analyse(tic, P, old['t0_btjd'], w, s.get('gaia', {}))
    return old, r


def main(argv=None):
    ap = argparse.ArgumentParser(prog='occulta binaries rerun', description=__doc__.split('Usage:')[0])
    ap.add_argument('tic', nargs='*', type=int, help='default: every analysed star')
    ap.add_argument('--divide', type=int, default=1, help='analyse at the adopted period divided by this integer')
    a = ap.parse_args(sys.argv[1:] if argv is None else list(argv))
    for tic in a.tic or sorted(int(p.stem) for p in eb_vet.OUT.glob('*.json') if p.stem.isdigit()):
        try:
            old, r = redo(tic, a.divide)
        except Exception as e:      # one failing star must not stop a pass over all of them
            print(tic, 'failed:', repr(e)[:120], flush=True); continue
        if 'period' not in r:
            print(tic, r.get('status')); continue
        print(tic, 'P %.7f -> %.7f' % (old['period'], r['period']), 'T0 %.5f -> %.5f (%.1f min)' % (
            old['t0_btjd'], r['t0_btjd'], ((r['t0_btjd'] - old['t0_btjd'] + r['period'] / 2) % r['period'] - r['period'] / 2) * 1440),
            'sec phase', (r.get('secondary') or {}).get('phase'), flush=True)


if __name__ == '__main__':
    main()
