"""Light-travel-time scan of eclipse times: fits T0 + P E + a sinusoid of period Pout to the individually fitted primary
eclipses (eb/<TIC>_eclipses.csv, QLP Sectors 99-104 excluded) over a grid of Pout, and lists the best periods and any
fixed periods given (e.g. a Gaia DR3 orbital period). Usage: occulta binaries lighttime-scan <TIC> <P_d> <T0_btjd> [Pout_fixed ...]"""
import numpy as np, pandas as pd, sys

from .. import paths


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    tic=int(argv[0]); P0=float(argv[1]); T00=float(argv[2])
    d=pd.read_csv(paths.BIN / 'eb' / f'{tic}_eclipses.csv'); d=d[(d.type=='P')&~d.unit.between(99,104)]
    d=d[d.t_err_d>0]
    t=d.t_btjd.values; E=np.round((t-T00)/P0); e=np.maximum(d.t_err_d.values, 20/86400)
    def fit(Pout):
        cols=[np.ones_like(t),E]
        if Pout: cols+= [np.sin(2*np.pi*t/Pout), np.cos(2*np.pi*t/Pout)]
        A=np.vstack(cols).T/e[:,None]; c,*_=np.linalg.lstsq(A,t/e,rcond=None)
        r=(t-np.vstack(cols).T@c); return np.sum((r/e)**2), c, r
    chi0,_,r0=fit(None); print('N',len(t),'chi2 linear',round(chi0,1),'rms min',round(np.std(r0)*1440,2))
    grid=np.exp(np.linspace(np.log(20),np.log(3000),20000))
    chis=np.array([fit(p)[0] for p in grid])
    o=np.argsort(chis)
    best=[];
    for i in o:
        if all(abs(np.log(grid[i]/b))>0.02 for b in best): best.append(grid[i])
        if len(best)==8: break
    for b in best:
        c2,c,r=fit(b); print('Pout %.1f chi2 %.1f amp %.2f min rms %.2f'%(b,c2,np.hypot(c[2],c[3])*1440,np.std(r)*1440))
    for pg in [float(x) for x in argv[3:]]:
        c2,c,r=fit(pg); print('FIXED Pout %.2f chi2 %.1f amp %.2f min rms %.2f'%(pg,c2,np.hypot(c[2],c[3])*1440,np.std(r)*1440))



if __name__ == '__main__':
    main()
