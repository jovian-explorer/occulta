"""Fit the eclipse times of an eclipsing binary with a light-travel-time term whose shape is fixed by a Gaia DR3
astrometric orbit (period, eccentricity, epoch of periastron, argument of periastron from the Thiele-Innes elements),
leaving only the linear ephemeris and the light-travel-time semi-amplitude K free (its sign absorbs the 180-degree
ambiguity of omega between the photocentre and the eclipsing pair).
Usage: occulta binaries lighttime <TIC> <P_d> <T0_btjd> <Pout_d> <ecc> <t_peri_gaia_d> <A> <B> <F> <G>
t_peri_gaia_d is Gaia's t_periastron (days from J2016.0 = JD 2457389.0)."""
import sys
import numpy as np, pandas as pd

def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    tic, P0, T00, Pout, ecc, tp, A, B, F, G = int(argv[0]), *map(float, argv[1:10])
    d = pd.read_csv(paths.BIN / 'eb' / f'{tic}_eclipses.csv'); d = d[(d.type == 'P') & ~d.unit.between(99, 104) & (d.t_err_d > 0)]
    t = d.t_btjd.values; E = np.round((t - T00) / P0); err = np.maximum(d.t_err_d.values, 20 / 86400)
    Tp = tp + 2457389.0 - 2457000.0                                   # BTJD of periastron
    w = 0.5 * (np.arctan2(B - F, A + G) + np.arctan2(-B - F, A - G))  # omega (photocentre), Thiele-Innes inversion


    def shape(t, Pout, Tp, w):
        M = 2 * np.pi * (t - Tp) / Pout
        Ea = M.copy()
        for _ in range(50):
            Ea = Ea - (Ea - ecc * np.sin(Ea) - M) / (1 - ecc * np.cos(Ea))
        nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(Ea / 2), np.sqrt(1 - ecc) * np.cos(Ea / 2))
        return (1 - ecc ** 2) / (1 + ecc * np.cos(nu)) * np.sin(nu + w) + ecc * np.sin(w)


    def fit(cols):
        X = np.vstack(cols).T
        c, *_ = np.linalg.lstsq(X / err[:, None], t / err, rcond=None)
        r = t - X @ c
        return float(np.sum((r / err) ** 2)), c, r


    chi_lin, _, r_lin = fit([np.ones_like(t), E])
    chi, c, r = fit([np.ones_like(t), E, shape(t, Pout, Tp, w)])
    print(f'N {t.size}; omega {np.degrees(w):.1f} deg; T_peri {Tp:.1f} BTJD')
    print(f'linear: chi2 {chi_lin:.1f}, rms {np.std(r_lin)*1440:.2f} min')
    print(f'Gaia-shape LTTE: chi2 {chi:.1f}, K {c[2]*1440:+.2f} min, rms {np.std(r)*1440:.2f} min')
    # let the outer period and periastron move within their Gaia errors (and a little beyond)
    best = (chi, Pout, Tp, c[2])
    for Pq in np.linspace(Pout - 3, Pout + 3, 61):
        for Tq in np.linspace(Tp - 25, Tp + 25, 101):
            x2, cq, _ = fit([np.ones_like(t), E, shape(t, Pq, Tq, w)])
            if x2 < best[0]:
                best = (x2, Pq, Tq, cq[2])
    print(f'best within Pout+-3 d, T_peri+-25 d: chi2 {best[0]:.1f} at Pout {best[1]:.2f} d, T_peri {best[2]:.1f}, K {best[3]*1440:+.2f} min')

    if len(argv) > 10:                                            # optional figure: O-C with the Gaia-shaped curve
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        k = c[2]
        oc = (t - c[0] - c[1] * E) * 1440
        tt = np.linspace(t.min() - 30, t.max() + 30, 4000)
        spoc = d.unit.values > 1000
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.plot(tt, k * shape(tt, Pout, Tp, w) * 1440, '-', color='0.5', lw=1, label='Gaia DR3 orbit (shape and phase fixed)')
        ax.errorbar(t[~spoc], oc[~spoc], err[~spoc] * 1440, fmt='o', ms=6, color='0.1', label='QLP')
        ax.errorbar(t[spoc], oc[spoc], err[spoc] * 1440, fmt='s', ms=6, mfc='none', color='0.1', label='SPOC 2-min')
        ax.set_xlabel('BTJD'); ax.set_ylabel('O - C (min)'); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(argv[10], dpi=110)
        # zoomed panels per TESS visit
        o = np.argsort(t); grp_id = np.empty(t.size, int); grp_id[o] = np.cumsum(np.r_[0, np.diff(t[o]) > 100])
        yrs = sorted(set(grp_id))
        fig, axs = plt.subplots(1, len(yrs), figsize=(3 * len(yrs), 3.2), sharey=True)
        for a, y in zip(np.atleast_1d(axs), yrs):
            m = grp_id == y
            x = np.linspace(t[m].min() - 15, t[m].max() + 15, 400)
            a.plot(x, k * shape(x, Pout, Tp, w) * 1440, '-', color='0.5')
            a.errorbar(t[m & ~spoc], oc[m & ~spoc], err[m & ~spoc] * 1440, fmt='o', ms=5, color='0.1')
            a.errorbar(t[m & spoc], oc[m & spoc], err[m & spoc] * 1440, fmt='s', ms=5, mfc='none', color='0.1')
            a.set_xlabel('BTJD')
        np.atleast_1d(axs)[0].set_ylabel('O - C (min)')
        fig.tight_layout(); fig.savefig(argv[10].replace('.png', '_zoom.png'), dpi=110)
        print('figures written')



if __name__ == '__main__':
    main()
