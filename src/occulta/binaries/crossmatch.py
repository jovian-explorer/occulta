"""Byproduct search: eclipsing binaries and other stellar signals that the planet vetting rejected and that no
catalogue lists yet.

Selects every candidate rejected for a stellar reason (stellar-size radius, depth above 25%, V shape, odd/even
difference, or secondary eclipse), cross-matches it within 10 arcsec against VSX (B/vsx), Gaia DR3 eclipsing
binaries (I/358/veb), the TESS EB catalogues (J/ApJS/258/16; J/A+A/691/A242) and SIMBAD variable classes, and then removes stars
with a catalogued variable within 63 arcsec (three TESS pixels), or a VSX variable at least 2 mag brighter within 180
arcsec, which could contaminate the aperture.
Writes binaries/new_variables.csv (all, with the catalogue that knows each) and new_variables_isolated.csv.
Usage: occulta binaries crossmatch"""
import sqlite3, sys, time
from pathlib import Path

from .. import paths

import astropy.units as u
from astropy.table import Table
from astroquery.xmatch import XMatch

ROOT = paths.HOME
STELLAR = ('secondary', 'odd/even', 'V-shaped', 'stellar-size', 'too deep')
CATS = {'VSX': 'vizier:B/vsx/vsx', 'GaiaEB': 'vizier:I/358/veb', 'TESS_EB': 'vizier:J/ApJS/258/16/tess-ebs',
        'OBAF_EB': 'vizier:J/A+A/691/A242/obaf-eb1', 'T10K_EB': 'vizier:J/ApJS/279/50/table3', 'SIMBAD': 'simbad'}
VAR_TYPES = ('EB', 'Variable', 'V*', 'Al', 'bL', 'WU', 'Puls', 'dS', 'RR', 'Ro', 'gD', 'El', 'SB')


def xmatch(tab, cat, radius):
    for _ in range(3):
        try:
            return XMatch.query(cat1=tab['tic', 'ra', 'dec'], cat2=cat, max_distance=radius * u.arcsec,
                                colRA1='ra', colDec1='dec')
        except Exception as e:
            print(cat, 'retry', repr(e)[:80]); time.sleep(10)
    sys.exit(f'{cat} cross-match failed')


REVIEW_STELLAR = ('eclipsing binary', 'v-shaped', 'stellar companion', 'grazing', 'eb,', 'eb)', 'eb ')
REVIEW_NOT = ('instrumental', 'systematic', 'pulsator', 'no transit visible', 'no clear dip', 'not a transit', 'variability')


def main(argv=None):
    con = sqlite3.connect(ROOT / 'data' / 'hunt.db', timeout=60)
    rows = []
    for t, p, n, st in con.execute("SELECT tic, period, notes, status FROM candidate WHERE status IN ('rejected', 'rejected_review') "
                                   "AND period > 0 AND tic < 2e10"):
        # candidates rejected after review are kept when the review found a stellar eclipse, not a systematic
        if st == 'rejected_review' and n and (any(x in n.lower() for x in REVIEW_STELLAR)
                                              and not any(x in n.lower() for x in REVIEW_NOT)):
            n = n + ' | ' + STELLAR[0]
        if n and any(x in n for x in STELLAR):
            r = con.execute('SELECT ra, dec, tmag, teff, rad FROM tic WHERE tic=?', (t,)).fetchone()
            if r and r[0] is not None:
                rows.append((t, p, *r))
    tab = Table(rows=rows, names=('tic', 'period', 'ra', 'dec', 'tmag', 'teff', 'rad'))
    hit, stype = {}, {}
    for name, cat in CATS.items():
        res = xmatch(tab, cat, 10)
        hit[name] = set(int(x) for x in res['tic'])
        if name == 'SIMBAD':
            for x in res:
                stype.setdefault(int(x['tic']), str(x['main_type']))
    out = []
    for r in tab:
        t = int(r['tic']); st = stype.get(t, '')
        known = [n for n in ('VSX', 'GaiaEB', 'TESS_EB', 'OBAF_EB', 'T10K_EB') if t in hit[n]] + ([f'SIMBAD:{st}'] if any(v in st for v in VAR_TYPES) else [])
        out.append((t, r['period'], r['ra'], r['dec'], r['tmag'], r['teff'], r['rad'], ';'.join(known) or 'NEW', st))
    allv = Table(rows=out, names=('tic', 'period_d', 'ra', 'dec', 'tmag', 'teff', 'rad', 'known_in', 'simbad_type'))
    allv.write(paths.BIN / 'new_variables.csv', overwrite=True)
    new = allv[allv['known_in'] == 'NEW']
    near = set()
    for cat in ('vizier:B/vsx/vsx', 'vizier:I/358/veb'):
        near |= set(int(v) for v in xmatch(new, cat, 63)['tic'])
    # a much brighter variable contaminates the aperture from farther away (HD 71433, 102 arcsec from
    # TYC 4861-1879-1): VSX variables within 180 arcsec whose maximum is at least 2 mag brighter than the target
    wide = xmatch(new, 'vizier:B/vsx/vsx', 180)
    tm = {int(t): float(m) for t, m in zip(new['tic'], new['tmag'])}
    for x in wide:
        try:
            if float(x['max']) < tm[int(x['tic'])] - 2:
                near.add(int(x['tic']))
        except (KeyError, ValueError, TypeError):
            continue
    iso = new[[int(v) not in near for v in new['tic']]]
    iso.write(paths.BIN / 'new_variables_isolated.csv', overwrite=True)
    print(f'{len(allv)} stellar rejects; {len(new)} not catalogued; {len(iso)} with no catalogued variable within 63 arcsec')


if __name__ == '__main__':
    main()
