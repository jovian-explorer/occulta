"""List every dip deeper than k times the scatter, per sector, for one star (from the cached light curve)."""
import sys
from pathlib import Path
import numpy as np
from scipy.ndimage import median_filter
from . import analysis as eb_vet
def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    tic = int(argv[0]); k = float(argv[1]) if len(argv) > 1 else 5
    T, F, U, *_ = eb_vet.load(tic)
    for u in np.unique(U):
        m = U == u; t = T[m]; f = F[m] / median_filter(F[m], size=301, mode='nearest') - 1
        s = eb_vet.lcsearch.mad_std(f)
        idx = np.where(f < -k * s)[0]
        groups = [g for g in np.split(idx, np.where(np.diff(idx) > 5)[0] + 1) if g.size > 3]
        print(u, f'{t[0]:.1f}-{t[-1]:.1f}', 'rms %.4f' % s, [(round(float(t[g][np.argmin(f[g])]), 3), round(float(f[g].min()), 4)) for g in groups])



if __name__ == '__main__':
    main()
