"""Gaia DR3 cross-match of reported eclipsing binaries: brightest source within 3 arcsec, its RUWE, non_single_star
flag, RV statistics, and any non-single-star solution (orbit, acceleration, spectroscopic). Writes eb/gaia_nss.json.
Usage: occulta binaries gaia"""
import json
from pathlib import Path

from .. import paths

import pandas as pd
from astroquery.gaia import Gaia

B = paths.BIN
ROOT = paths.HOME
NSS = ['nss_two_body_orbit', 'nss_acceleration_astro', 'nss_non_linear_spectro', 'nss_vim_fl']
KEEP = ['nss_solution_type', 'period', 'period_error', 'eccentricity', 'eccentricity_error', 't_periastron',
        'a_thiele_innes', 'b_thiele_innes', 'f_thiele_innes', 'g_thiele_innes', 'semi_amplitude_primary',
        'center_of_mass_velocity', 'accel_ra', 'accel_dec', 'significance', 'goodness_of_fit', 'parallax']


def main(argv=None):
    rv = pd.read_csv(B / 'eb' / 'review.csv')
    tics = [int(t) for t in rv[rv.iloc[:, 1].isin(['A', 'B'])].iloc[:, 0]]
    out = {}
    for tic in tics:
        s = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())['star']
        q = (f"select DISTANCE(POINT({s['ra']},{s['dec']}),POINT(ra,dec))*3600 as sep, source_id, phot_g_mean_mag as gmag, "
             f"bp_rp, parallax, parallax_error, ruwe, radial_velocity, rv_amplitude_robust, rv_nb_transits, non_single_star "
             f"from gaiadr3.gaia_source where 1=CONTAINS(POINT({s['ra']},{s['dec']}),CIRCLE(ra,dec,3./3600)) order by phot_g_mean_mag")
        r = Gaia.launch_job(q).get_results()
        if not len(r):
            continue
        g = {k: (None if str(r[k][0]) == '--' else (int(r[k][0]) if k in ('source_id', 'non_single_star', 'rv_nb_transits') else float(r[k][0])))
             for k in r.colnames}
        g['nss'] = {}
        if g['non_single_star']:
            for tab in NSS:
                n = Gaia.launch_job(f"select * from gaiadr3.{tab} where source_id={g['source_id']}").get_results()
                if len(n):
                    g['nss'][tab] = {k: (None if str(n[k][0]) == '--' else (str(n[k][0]) if k == 'nss_solution_type' else float(n[k][0])))
                                     for k in KEEP if k in n.colnames}
        out[tic] = g
        nss = '; '.join(f"{v.get('nss_solution_type')} P={v.get('period')}" for v in g['nss'].values())
        print(tic, 'G %.2f RUWE %.2f NSS %s RVamp %s (%s) %s' % (g['gmag'], g['ruwe'] or 0, g['non_single_star'],
              g['rv_amplitude_robust'], g['rv_nb_transits'], nss), flush=True)
    (B / 'eb' / 'gaia_nss.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
