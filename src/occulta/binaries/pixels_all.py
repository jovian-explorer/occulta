"""Assign the source of the eclipses at the pixel level (pixels.py) for many eclipsing binaries at once, each at its
adopted ephemeris (binaries/eb/<TIC>.json), and collect the verdicts in binaries/eb/diffimg/verdicts.csv.
The stars are those given on the command line, or, when none is given, every star of class A or B in
binaries/eb/review.csv. Usage: occulta binaries pixels-all [TIC ...]"""
import json
import sys

import pandas as pd

from .. import paths
from . import pixels


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    B, ROOT = paths.BIN, paths.HOME
    if argv:
        tics, cls = [int(t) for t in argv], {}
    else:
        rv = pd.read_csv(B / 'eb' / 'review.csv')
        sel = rv[rv.iloc[:, 1].isin(['A', 'B'])]
        tics = [int(t) for t in sel.iloc[:, 0]]
        cls = dict(zip(tics, sel.iloc[:, 1]))
    rows = []
    for tic in tics:
        e = json.loads((B / 'eb' / f'{tic}.json').read_text())
        st = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())['star']
        w = min(e.get('width_h') or 3.0, e['period'] * 24 / 5)
        out = B / 'eb' / 'diffimg' / f'{tic}.json'
        if not out.exists():
            pixels.main([str(tic), repr(e['period']), repr(e['t0_btjd']), repr(w), repr(st['ra']), repr(st['dec'])])
        v = json.loads(out.read_text()) if out.exists() else {}
        rows.append({'tic': tic, 'class': cls.get(tic), 'on_target': v.get('on_target'), 'n_best': v.get('n_best'),
                     'n_sectors': v.get('n_sectors'), 'best_source': v.get('best_source'),
                     'target_source': v.get('target_source'), 'median_chi2_ratio_next': v.get('median_chi2_ratio_next'),
                     'target_chi2_rel_when_not_best': v.get('median_target_chi2_rel_when_not_best')})
        print(rows[-1], flush=True)
        pd.DataFrame(rows).to_csv(B / 'eb' / 'diffimg' / 'verdicts.csv', index=False)


if __name__ == '__main__':
    main()
