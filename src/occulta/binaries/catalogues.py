"""Is this star already catalogued as variable or as a planet host?  lookup(tic, ra, dec) -> dict.
Searched: VSX (B/vsx/vsx), GCVS (B/gcvs), Gaia DR3 eclipsing binaries (I/358/veb) and variability classes
(I/358/vclassre), the TESS EB catalogue of Prsa et al. 2022 (J/ApJS/258/16), the TESS OBAF eclipsing
binaries (J/A+A/691/A242), the TESS Ten Thousand Catalog (J/ApJS/279/50, vetted table3; unvetted table2 is reported only), SIMBAD, and the ExoFOP TOI/CTOI lists.
Targets within 10 arcsec count as the star itself; variables within 63 arcsec (3 TESS pixels) are listed as
neighbours that could contaminate the aperture."""
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import os
# the CDS server is sometimes unreachable for hours; the CfA mirror holds the same catalogues
Vizier.VIZIER_SERVER = os.environ.get('VIZIER_SERVER', 'vizier.cfa.harvard.edu')

VIZ = {'VSX': 'B/vsx/vsx', 'GCVS': 'B/gcvs/gcvs_cat', 'GaiaDR3_EB': 'I/358/veb', 'GaiaDR3_class': 'I/358/vclassre',
       'TESS_EB': 'J/ApJS/258/16/tess-ebs', 'TESS_OBAF_EB': 'J/A+A/691/A242/obaf-eb1',
       'TESS_T10K_EB': 'J/ApJS/279/50/table3', 'TESS_T10K_unvetted': 'J/ApJS/279/50/table2'}


def lookup(tic, ra, dec):
    c = SkyCoord(ra, dec, unit='deg')
    out = {'self': {}, 'neighbours_63arcsec': []}
    v = Vizier(columns=['*', '_r'], row_limit=20)
    for name, cat in VIZ.items():
        try:
            res = v.query_region(c, radius=63 * u.arcsec, catalog=cat)
        except Exception as e:
            out['self'][name] = f'query failed: {e!r}'[:80]
            continue
        rows = res[0] if len(res) else []
        for r in rows:
            sep = float(r['_r'])   # arcsec, the unit of the query radius
            label = {k: str(r[k]) for k in r.colnames if k in ('Name', 'Type', 'Period', 'GCVS', 'VarType', 'Class',
                                                                  'Source', 'TIC', 'Per')}
            if sep <= 10:
                out['self'][name] = label
            else:
                out['neighbours_63arcsec'].append({'catalogue': name, 'sep_arcsec': round(sep, 1), **label})
    try:
        from astroquery.simbad import Simbad
        Simbad.add_votable_fields('otype')
        s = Simbad.query_region(c, radius=10 * u.arcsec)
        if s is not None and len(s):
            cols = {k.lower(): k for k in s.colnames}
            out['self']['SIMBAD'] = {'id': str(s[0][cols['main_id']]), 'otype': str(s[0][cols['otype']])}
    except Exception as e:
        out['self']['SIMBAD'] = f'query failed: {e!r}'[:80]
    try:
        j = requests.get('https://exofop.ipac.caltech.edu/tess/target.php', params={'id': tic, 'json': ''}, timeout=60).json()
        if j.get('tois') or j.get('ctois'):
            out['self']['ExoFOP'] = {'tois': [t.get('toi') for t in j.get('tois', [])],
                                     'ctois': [t.get('ctoi') for t in j.get('ctois', [])]}
    except Exception as e:
        out['self']['ExoFOP'] = f'query failed: {e!r}'[:80]
    var = {'VSX', 'GCVS', 'GaiaDR3_EB', 'TESS_EB', 'TESS_OBAF_EB', 'TESS_T10K_EB', 'ExoFOP'}   # T10K table2 is an unvetted list: reported, not counted
    st = out['self']
    out['new'] = not any(k in var and isinstance(st[k], dict) for k in st) and \
        not ('SIMBAD' in st and isinstance(st['SIMBAD'], dict) and any(
            x in st['SIMBAD']['otype'] for x in ('EB', 'V*', 'Var', 'Pu', 'dS', 'RR', 'Ro', 'El', 'SB', 'gD', 'BY', 'RS')))
    return out
