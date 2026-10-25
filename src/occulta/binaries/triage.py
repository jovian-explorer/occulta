"""One-line triage of analysed eclipsing-binary candidates: period, depths, secondary phase, widths, host radius,
Roche-lobe ratio of the host for a twin-mass orbit, Gaia RUWE, centroid shifts. Usage: occulta binaries triage TIC ..."""
import json, sys
from pathlib import Path

from .. import paths
import numpy as np

B = paths.BIN
ROOT = paths.HOME
def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    for tic in [int(x) for x in argv]:
        p = B / 'eb' / f'{tic}.json'
        if not p.exists():
            print(tic, 'not analysed'); continue
        e = json.loads(p.read_text())
        if e.get('status') != 'ok':
            print(tic, e.get('status')); continue
        st = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())['star']
        P = e['period']; R = st.get('rad'); M = st.get('mass') or 1.0
        roche = None
        if R and R == R:
            a = ((2 * (M if M == M else 1.0)) * (P / 365.25) ** 2) ** (1 / 3) * 215.03
            roche = R / (0.38 * a)
        sec = e.get('secondary') or {}
        cen = max((c['sigma'] for c in e.get('centroid', []) if c['unit'] < 99 or c['unit'] > 1000), default=0)
        print(f"{tic} T{st['tmag']:.2f} Teff {st.get('teff')} R {R if R == R else '-'} | P {P:.5f} dbl {e['doubled']} "
              f"D1 {100*e['primary_depth']:.2f}% D2 {100*e['secondary_depth']:.2f}% ph2 {sec.get('phase', e['secondary_phase_scan']):.3f} "
              f"w {e['width_h']:.1f}h | Roche {roche if roche is None else round(roche, 2)} RUWE {(e.get('gaia') or {}).get('ruwe')} "
              f"cen_max_sig {cen} n_ecl {e['eclipses_seen']}/{e['eclipses_covered']}")



if __name__ == '__main__':
    main()
