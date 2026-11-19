"""Which vetted candidates are already known?  Period-aware cross-match of the ready candidates with the planet and
eclipsing-binary catalogues.

A candidate is "recovered" (independently found by occulta, but not new) when a catalogue lists the same star, by TIC
number or within MATCH_ARCSEC of its position, with a period equal to the candidate period P, to 2P or to P/2 within
REL_TOL (a BLS or a QLP search often locks on a harmonic, and the catalogue period may be the half or the double).
The catalogues are read from one folder of plain files (fetch() writes them):

    planet:  NASA Exoplanet Archive confirmed planets (pscomppars), TOI table, K2 candidates and Kepler cumulative
             KOI table; ExoFOP TOI and CTOI lists; the internal known table of the survey database
    EB:      Villanova TESS EB catalogue (Prsa et al. 2022), Kepler EB catalogue (Kirk et al. 2016), TESS Ten Thousand
             Catalog (Kostov et al. 2025), TESS OBAF eclipsing binaries, Gaia DR3 eclipsing binaries (period = 1 /
             frequency), VSX eclipsing types and SIMBAD eclipsing-binary types (no period: star-level only)

Matches are returned as a table with one row per (candidate, catalogue entry); decide() reduces it to one verdict per
candidate: 'recovered' (period-compatible match), 'known_star' (the star is listed with another period; kept as
a note, never as a new detection when the catalogue says eclipsing binary).

    occulta recovered [--catalogues DIR] [--out FILE.csv]
"""
import csv
import io
import re
import sys
from pathlib import Path

import numpy as np

REL_TOL = 0.015            # fractional period agreement
MATCH_ARCSEC = 10.0        # within this separation the catalogue star is the candidate star (3 TESS pixels would blend)
RELATIONS = (('P', 1.0), ('2P', 2.0), ('P/2', 0.5))
EB_VSX = ('E', 'EA', 'EB', 'EW', 'EC', 'ED', 'ELL', 'GS', 'AR', 'D', 'DM', 'DS', 'DW', 'K', 'KE', 'KW', 'SD', 'DC', 'PN')
EB_SIMBAD = ('EB*', 'Al*', 'bL*', 'WU*', 'EP*', 'EB?', 'Al?', 'SB*', 'SB?')


def period_relation(p, q, tol=REL_TOL):
    """'P', '2P' or 'P/2' when the candidate period p equals q, 2q or q/2 within tol, otherwise None.
    '2P' reads: the candidate period is twice the catalogue period."""
    try:
        p, q = float(p), float(q)
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(p) and np.isfinite(q)) or p <= 0 or q <= 0:
        return None
    best = None
    for name, f in RELATIONS:
        err = abs(p / q - f) / f
        if err < tol and (best is None or err < best[0]):
            best = (err, name)
    return best[1] if best else None


def phase_consistent(t0_a, p_a, t0_b, p_b, rel, dur_d=0.2):
    """Whether two ephemerides (BTJD) put a transit at the same time: at the shorter common period the transit times
    must agree to within dur_d. None when an epoch is missing."""
    try:
        t0_a, t0_b, p_a, p_b = float(t0_a), float(t0_b), float(p_a), float(p_b)
    except (TypeError, ValueError):
        return None
    if not all(np.isfinite([t0_a, t0_b, p_a, p_b])):
        return None
    p = min(p_a, p_b)
    d = (t0_a - t0_b + 0.5 * p) % p - 0.5 * p
    return bool(abs(d) < max(dur_d, 0.02 * p))


def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def _read_asu(path):
    """VizieR ASU tab-separated file: comment lines, a header, a units line and a dashes line."""
    lines = [l for l in Path(path).read_text(encoding='utf8', errors='ignore').splitlines() if l.strip() and not l.startswith('#')]
    head = [h.strip() for h in lines[0].split('\t')]
    rows = [[x.strip() for x in l.split('\t')] for l in lines[3:] if not l.startswith('-')]
    return [dict(zip(head, r)) for r in rows if len(r) == len(head)]


def _ads(html):
    m = re.search(r'abs/([^/]+)/abstract', html or '')
    return m.group(1) if m else ''


def _short_ref(html):
    return re.sub(r'<[^>]+>', '', html or '').strip()


def load(folder):
    """List of catalogue entries: dicts with cat, kind ('planet', 'planet candidate', 'EB', 'BD', 'variable'), id (TIC
    or None), ra, dec, period, t0 (BTJD), cls (catalogue class text), ref (reference text) and bib (ADS bibcode)."""
    d = Path(folder)
    out = []

    def add(**k):
        k.setdefault('id', None); k.setdefault('period', None); k.setdefault('t0', None); k.setdefault('cls', '')
        k.setdefault('ref', ''); k.setdefault('bib', ''); k.setdefault('name', '')
        out.append(k)

    def csvrows(name):
        p = d / name
        return list(csv.DictReader(io.StringIO(p.read_text(encoding='utf8', errors='ignore')))) if p.exists() else []

    for r in csvrows('pscomppars.csv'):
        t = re.sub(r'\D', '', r.get('tic_id') or '')
        add(cat='NASA Exoplanet Archive (confirmed)', kind='planet', id=int(t) if t else None, ra=_num(r['ra']),
            dec=_num(r['dec']), period=_num(r['pl_orbper']), cls='confirmed planet', name=r['pl_name'],
            ref=_short_ref(r['disc_refname']), bib=_ads(r['disc_refname']))
    for r in csvrows('toi_nea.csv'):
        add(cat='NASA Exoplanet Archive (TOI)', kind='planet candidate', id=int(r['tid']), ra=_num(r['ra']),
            dec=_num(r['dec']), period=_num(r['pl_orbper']), t0=(_num(r['pl_tranmid']) or 0) - 2457000 or None,
            cls=f"TOI {r['tfopwg_disp']}", name=r['toidisplay'], ref='TESS Objects of Interest', bib='')
    for r in csvrows('exofop_toi.csv'):
        add(cat='ExoFOP TOI', kind='planet candidate', id=int(r['TIC ID']), ra=_num(r.get('RA')), dec=_num(r.get('Dec')),
            period=_num(r['Period (days)']), t0=(_num(r['Epoch (BJD)']) or 0) - 2457000 or None,
            cls=f"TOI {r.get('TFOPWG Disposition', '')}", name=f"TOI-{r['TOI']}", ref='ExoFOP TESS TOI list')
    for r in csvrows('exofop_ctoi.csv'):
        add(cat='ExoFOP CTOI', kind='planet candidate', id=int(r['TIC ID']), ra=_num(r['RA']), dec=_num(r['Dec']),
            period=_num(r['Period (days)']), t0=(_num(r['Transit Epoch (BJD)']) or 0) - 2457000 or None,
            cls=f"CTOI {r.get('TFOPWG Disposition', '')}", name=f"CTOI {r['CTOI']}",
            ref=f"ExoFOP CTOI ({r.get('User', '')})")
    for r in csvrows('arxiv_entries.csv'):                 # tables copied from arXiv sources (columns cat kind id period cls ref bib)
        add(cat=r['cat'], kind=r['kind'], id=int(r['id']), ra=None, dec=None, period=_num(r['period']), cls=r['cls'],
            name=f"TIC {r['id']}", ref=r['ref'], bib=r['bib'])
    for r in csvrows('k2pandc.csv'):
        add(cat='K2 planets and candidates', kind='planet candidate', ra=_num(r['ra']), dec=_num(r['dec']),
            period=_num(r['pl_orbper']), cls=r['disposition'].lower(), name=r['pl_name'], ref=r['disp_refname'])
    for r in csvrows('koi.csv'):
        add(cat='Kepler KOI', kind='planet candidate', ra=_num(r['ra']), dec=_num(r['dec']), period=_num(r['koi_period']),
            cls=r['koi_disposition'].lower(), name=r['kepoi_name'], ref='Kepler cumulative KOI table')
    if (d / 'prsa2022.tsv').exists():
        for r in _read_asu(d / 'prsa2022.tsv'):
            add(cat='TESS EB catalogue', kind='EB', id=int(r['TIC']), ra=_num(r['RAJ2000']), dec=_num(r['DEJ2000']),
                period=_num(r['Per']), cls='EB', name=f"TIC {int(r['TIC'])}", ref='Prša et al. 2022', bib='2022ApJS..258...16P')
    if (d / 'kirk2016.tsv').exists():
        for r in _read_asu(d / 'kirk2016.tsv'):
            add(cat='Kepler EB catalogue', kind='EB', ra=_num(r['_RA']), dec=_num(r['_DE']), period=_num(r['Per']),
                cls='EB', name=f"KIC {r['KIC'].strip()}", ref='Kirk et al. 2016', bib='2016AJ....151...68K')
    for f, tab, dis in (('t10k_vetted.tsv', 'vetted', 'EB'), ('t10k_unvetted.tsv', 'unvetted', 'EB candidate')):
        if (d / f).exists():
            for r in _read_asu(d / f):
                if _num(r.get('Per')) is None:
                    continue
                add(cat=f'TESS Ten Thousand Catalog ({tab})', kind='EB', id=int(r['TIC']), ra=_num(r['RAJ2000']),
                    dec=_num(r['DEJ2000']), period=_num(r['Per']), cls=dis, name=f"TIC {int(r['TIC'])}",
                    ref='Kostov et al. 2025', bib='2025ApJS..279...50K')
    if (d / 'obaf.tsv').exists():
        for r in _read_asu(d / 'obaf.tsv'):
            per = _num(r.get('Per') or r.get('Porb') or r.get('P'))
            tic = _num(str(r.get('TIC')).strip())
            add(cat='TESS OBAF eclipsing binaries', kind='EB', id=int(tic) if tic else None,
                ra=_num(r.get('RAJ2000')), dec=_num(r.get('DEJ2000')), period=per, cls='EB',
                name=f"TIC {int(tic)}" if tic else '', ref='IJspeert et al. 2024', bib='2024A&A...691A.242I')
    for r in csvrows('gaia_nb.csv'):
        fq = _num(r.get('ebfreq'))
        if fq and _num(r.get('sep_arcsec')) is not None and float(r['sep_arcsec']) <= 3.0:
            add(cat='Gaia DR3 eclipsing binaries', kind='EB', ra=None, dec=None, period=1 / fq, cls='EB',
                name=f"Gaia DR3 {r['source_id']}", ref='Mowlavi et al. 2023', bib='2023A&A...674A..16M',
                sid=int(r['sid']))
    p = d / 'vsx.tsv'
    if p.exists():
        for l in p.read_text(encoding='utf8', errors='ignore').splitlines():
            f = l.split('\t')
            if len(f) < 7 or f[1] == 'ERR':
                continue
            typ = f[2].strip()
            base = typ.split('|')[0].split('+')[0]
            kind = 'EB' if base in EB_VSX else 'variable'
            add(cat='VSX', kind=kind, ra=_num(f[4]), dec=_num(f[5]), period=_num(f[3]), cls=typ,
                name=f[1].strip(), ref='AAVSO VSX', sid=int(f[0]))
    p = d / 'simbad.csv'
    if p.exists():
        for r in csv.DictReader(io.StringIO(p.read_text(encoding='utf8'))):
            ot = r['otype']
            kind = 'EB' if ot in EB_SIMBAD else ('planet candidate' if ot in ('Pl', 'Pl?') else 'variable')
            m = re.search(r'(TIC|EPIC|KIC) (\d+)', r['id'])
            add(cat='SIMBAD', kind=kind, ra=_num(r['ra']), dec=_num(r['dec']), cls=ot, name=r['main_id'],
                ref='SIMBAD', sidname=r['id'])
    return out


def match(cands, entries, arcsec=MATCH_ARCSEC):
    """One row per (candidate, catalogue entry) on the same star. cands: dicts with tic, name, ra, dec, period, t0.
    Entries match by TIC number, by candidate sid (VSX and Gaia, queried per star), by name (SIMBAD), or by position."""
    byid, bysid, byname = {}, {}, {}
    pos = []
    for e in entries:
        if e.get('id'):
            byid.setdefault(e['id'], []).append(e)
        if e.get('sid'):
            bysid.setdefault(e['sid'], []).append(e)
        if e.get('sidname'):
            byname.setdefault(e['sidname'], []).append(e)
        if e.get('ra') is not None and e.get('dec') is not None and not e.get('sid') and not e.get('sidname'):
            pos.append(e)
    pra = np.radians([e['ra'] for e in pos])
    pde = np.radians([e['dec'] for e in pos])
    rows = []
    for c in cands:
        hits = {}
        for e in byid.get(c['tic'], []):
            hits[id(e)] = (e, 'TIC', 0.0)
        for e in bysid.get(c['tic'], []):
            hits[id(e)] = (e, 'coord', None)
        for e in byname.get(c['name'], []):
            hits[id(e)] = (e, 'name', 0.0)
        if c.get('ra') is not None and len(pos):
            ra, de = np.radians(c['ra']), np.radians(c['dec'])
            cs = np.sin(de) * np.sin(pde) + np.cos(de) * np.cos(pde) * np.cos(ra - pra)
            sep = np.degrees(np.arccos(np.clip(cs, -1, 1))) * 3600
            for i in np.where(sep <= arcsec)[0]:
                e = pos[i]
                if id(e) not in hits:
                    hits[id(e)] = (e, 'coord', float(sep[i]))
        for e, how, sep in hits.values():
            rel = period_relation(c['period'], e['period'])
            rows.append({'tic': c['tic'], 'cat': e['cat'], 'kind': e['kind'], 'how': how, 'sep': sep,
                         'cat_name': e['name'], 'cat_period': e['period'], 'cat_class': e['cls'], 'ref': e['ref'],
                         'bib': e['bib'], 'relation': rel,
                         'phase_ok': phase_consistent(c.get('t0'), c['period'], e.get('t0'), e['period'] or 0, rel)
                         if rel else None})
    return rows


PRIORITY = {'NASA Exoplanet Archive (confirmed)': 0, 'ExoFOP TOI': 1, 'NASA Exoplanet Archive (TOI)': 1, 'ExoFOP CTOI': 2}


def decide(rows):
    """{tic: verdict dict} from match() rows. 'recovered' when an entry has a period relation to the candidate, with
    confirmed planets, then TOIs, then CTOIs, then EB catalogues preferred as the reported entry; 'known_eb_star'
    when only a catalogue of eclipsing binaries lists the star (another period, or none); 'known_star' otherwise."""
    by = {}
    for r in rows:
        by.setdefault(r['tic'], []).append(r)
    out = {}
    for tic, rs in by.items():
        rec = [r for r in rs if r['relation']]
        # a period match to a variable that is not eclipsing (a pulsator) is not a recovery of a transit or eclipse
        rec = [r for r in rec if r['kind'] in ('planet', 'planet candidate', 'EB')]
        if rec:
            rec.sort(key=lambda r: (PRIORITY.get(r['cat'], 5), r['sep'] if r['sep'] is not None else 0))
            r = rec[0]
            out[tic] = dict(verdict='recovered', match=r, others=sorted({x['cat'] for x in rec}))
            continue
        eb = [r for r in rs if r['kind'] == 'EB']
        if eb:
            out[tic] = dict(verdict='known_eb_star', match=eb[0], others=sorted({x['cat'] for x in eb}))
        elif any(r['kind'] in ('planet', 'planet candidate') for r in rs):
            m = [r for r in rs if r['kind'] in ('planet', 'planet candidate')]
            out[tic] = dict(verdict='known_star', match=m[0], others=sorted({x['cat'] for x in m}))
    return out


def main(argv=None):
    import argparse
    import sqlite3
    from . import paths
    ap = argparse.ArgumentParser(prog='occulta recovered')
    ap.add_argument('--catalogues', default=str(paths.HOME / 'xmatch' / 'cat'))
    ap.add_argument('--out', default=str(paths.HOME / 'xmatch' / 'matches.csv'))
    a = ap.parse_args(argv)
    from . import sources
    con = sqlite3.connect(paths.HOME / 'data' / 'hunt.db')
    cands = []
    for tic, kind, p, t0, ot in con.execute("SELECT tic, kind, period, t0, objtype FROM candidate WHERE status='ready'"):
        m = sources.mission_of(tic)
        pos = con.execute('SELECT ra, dec FROM tic WHERE tic=?' if m == 'TESS' else 'SELECT ra, dec FROM star WHERE sid=?',
                          (tic,)).fetchone()
        cands.append({'tic': tic, 'name': sources.name_of(tic), 'ra': pos[0] if pos else None,
                      'dec': pos[1] if pos else None, 'period': p, 't0': t0, 'objtype': ot})
    entries = load(a.catalogues)
    for tic, src, per in con.execute('SELECT tic, src, period FROM known'):        # the survey's own copy of the lists
        entries.append({'cat': 'survey known table', 'kind': 'planet candidate', 'id': tic, 'ra': None, 'dec': None,
                        'period': per, 't0': None, 'cls': src, 'ref': 'ExoFOP/NASA lists at survey time',
                        'bib': '', 'name': src})
    rows = match(cands, entries)
    with open(a.out, 'w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    v = decide(rows)
    n = {}
    for x in v.values():
        n[x['verdict']] = n.get(x['verdict'], 0) + 1
    print(len(cands), 'candidates;', n)


if __name__ == '__main__':
    main()
