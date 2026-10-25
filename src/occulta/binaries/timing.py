"""O-C diagram of the per-eclipse times from eclipses (binaries/eb/<TIC>_eclipses.csv), with a linear
ephemeris fitted to all times outside QLP Sectors 99-104 (timestamps not yet redelivered in corrected form; provisional), which are drawn grey.
Usage: occulta binaries timing <TIC> <P> <T0> <out.png>   -> prints the ephemeris and the mean O-C per sector"""
import sys
import numpy as np
from .. import paths
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time

def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    tic, P, t0, out = argv[0], float(argv[1]), float(argv[2]), argv[3]
    d = np.atleast_1d(np.genfromtxt(paths.BIN / 'eb' / f'{tic}_eclipses.csv', delimiter=',', names=True, dtype=None, encoding=None))
    p = d[d['type'] == 'P']
    E = np.round((p['t_btjd'] - t0) / P)
    bad = (p['unit'] >= 99) & (p['unit'] <= 104)
    A = np.vstack([np.ones((~bad).sum()), E[~bad]]).T
    c, *_ = np.linalg.lstsq(A, p['t_btjd'][~bad], rcond=None)
    oc = (p['t_btjd'] - (c[0] + c[1] * E)) * 1440
    err = p['t_err_d'] * 1440
    spoc = p['unit'] > 1000
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.errorbar(p['t_btjd'][~bad & ~spoc], oc[~bad & ~spoc], err[~bad & ~spoc], fmt='o', ms=5, color='0.1', label='QLP')
    ax.errorbar(p['t_btjd'][spoc], oc[spoc], err[spoc], fmt='s', ms=5, mfc='white', color='0.1', label='SPOC 2-min')
    ax.errorbar(p['t_btjd'][bad], oc[bad], err[bad], fmt='^', ms=5, color='0.55', label='QLP S99-104 (provisional timestamps)')
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('BTJD (BJD - 2457000)'); ax.set_ylabel('O - C (min)'); ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(out, dpi=100)
    print('T0 %.5f P %.7f' % (c[0], c[1]))
    for u in sorted(set(p['unit'].tolist())):
        m = p['unit'] == u
        print(u, Time(np.median(p['t_btjd'][m]) + 2457000, format='jd').iso[:7], 'n', m.sum(), 'O-C %.1f' % oc[m].mean(),
              'se %.1f' % (oc[m].std() / np.sqrt(m.sum())), 'depth %.2f%%' % (100 * np.median(p['depth'][m])))



if __name__ == '__main__':
    main()
