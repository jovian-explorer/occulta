"""Object type of every vetted candidate. The stellar parameters (radius, effective temperature, mass, TESS magnitude,
contamination ratio) are filled from TIC v8.2 at MAST, and from the star table for Kepler and K2 stars; the companion
radius follows as Rp = sqrt(depth) * R_star. A candidate is a planet below 2 R_J, a brown dwarf or low-mass star
(ambiguous) between 2 and 3 R_J, and an eclipsing binary above 3 R_J or when the vetting found a secondary eclipse or an
odd/even depth difference above 3 sigma. A host larger than 3 R_sun is flagged as a giant. Each host is also
cross-matched with Gaia DR3: a two-body orbit at the transit period (or half or twice it) makes the candidate an
eclipsing binary, with the companion mass from the spectroscopic mass function at sin i = 1; RUWE above 1.4, an orbit at
another period, or an astrometric acceleration is noted. The same holds for the Gaia DR3 photometric eclipsing-binary
catalogue (vari_eclipsing_binary), from Gaia's own sparse light curves. Radial velocities test the companion mass at
the transit period (sin i = 1): Gaia DR3 RVS variability (K ~ half the robust amplitude, for hosts that fail the
constant-RV test) and the APOGEE DR17 visit scatter (K ~ sqrt(2) vscatter); a companion above 0.08 M_sun makes the
candidate an eclipsing binary. For eclipses deeper than 3000 ppm, the ZTF light curve (2018 onward, dec > -30) is
searched in a narrow window around the transit period, and its odd/even depths and phase-0.5 depth are tested as in
the vetting; a ZTF depth more than 1.5 times the TESS/Kepler depth (contaminating light blended in the large pixels)
replaces it in the radius. The label goes into the new candidate column objtype and the reason into notes.

    occulta classify [--status ready] [--no-gaia] [--no-external] [TIC ...]
"""
import argparse
import json
import time

import numpy as np

from . import sources
from . import survey

RJ_RSUN = 71492.0 / 695700.0       # equatorial radius of Jupiter in solar radii
PLANET_MAX_RJ, BD_MAX_RJ = 2.0, 3.0
GIANT_RSUN = 3.0
OE_SIGMA = SEC_SIGMA = 3.0          # the cuts of planets.vet


def fill_star(con, sids):
    """TIC v8.2 parameters of Kepler and K2 stars, from the nearest TIC entry (within 3 arcsec) at the position that
    MAST resolves for the KIC or EPIC name; kept in the star table."""
    from astroquery.mast import Catalogs
    for sid in sids:
        try:
            tab = Catalogs.query_object(sources.name_of(sid), catalog='Tic', radius=3 / 3600)
        except Exception as e:                   # name not resolved
            print(f'{sources.name_of(sid)}: {e}')
            continue
        if not len(tab):
            continue
        r = tab[int(np.argmin(tab['dstArcSec']))]
        v = lambda k: None if np.ma.is_masked(r[k]) else float(r[k])
        con.execute('INSERT OR REPLACE INTO star(sid, mag, teff, rad, mass, ra, dec, dist) VALUES(?,?,?,?,?,?,?,?)',
                    (sid, v('Tmag'), v('Teff'), v('rad'), v('mass'), v('ra'), v('dec'), v('d')))
        con.commit()


def fill_tic(con, tics):
    """TIC v8.2 parameters of the given TESS stars that have no radius yet, in batches."""
    from astroquery.mast import Catalogs
    fill_star(con, [t for t in tics if sources.mission_of(t) != 'TESS'])
    tics = [t for t in tics if sources.mission_of(t) == 'TESS']
    for i in range(0, len(tics), 200):
        tab = Catalogs.query_criteria(catalog='Tic', ID=tics[i:i + 200])
        for r in tab:
            v = lambda k: None if np.ma.is_masked(r[k]) else float(r[k])
            con.execute('INSERT INTO tic(tic, ra, dec, tmag, teff, rad, mass, dist, gaia, contratio) VALUES(?,?,?,?,?,?,?,?,?,?) '
                        'ON CONFLICT(tic) DO UPDATE SET tmag=excluded.tmag, teff=excluded.teff, rad=excluded.rad, '
                        'mass=excluded.mass, contratio=excluded.contratio, ra=COALESCE(tic.ra, excluded.ra), '
                        'dec=COALESCE(tic.dec, excluded.dec), dist=COALESCE(tic.dist, excluded.dist), gaia=COALESCE(tic.gaia, excluded.gaia)',
                        (int(r['ID']), v('ra'), v('dec'), v('Tmag'), v('Teff'), v('rad'), v('mass'), v('d'),
                         None if np.ma.is_masked(r['GAIA']) else str(r['GAIA']), v('contratio')))
        con.commit()


def stellar_radius(con, tic):
    if sources.mission_of(tic) == 'TESS':
        r = con.execute('SELECT rad FROM tic WHERE tic=?', (tic,)).fetchone()
    else:
        r = con.execute('SELECT rad FROM star WHERE sid=?', (tic,)).fetchone()
    return r[0] if r and r[0] and r[0] > 0 else None


def vet_tests(tic):
    """Odd/even and secondary-eclipse significance from the vetting summary, when there is one."""
    p = survey.FIND / sources.name_of(tic).replace(' ', '') / 'summary.json'
    if not p.exists():
        return np.nan, np.nan
    s = json.loads(p.read_text())
    f = lambda x: float(x) if x is not None else np.nan
    return f(s.get('odd_even_sigma')), f(s.get('secondary_sigma'))


def label(depth_ppm, rstar, oe, sec):
    """(objtype, reason) of one candidate."""
    giant = rstar is not None and rstar > GIANT_RSUN
    if np.isfinite(sec) and sec > SEC_SIGMA:
        t, why = 'eclipsing binary', f'secondary eclipse at {sec:.1f} sigma'
    elif np.isfinite(oe) and oe > OE_SIGMA:
        t, why = 'eclipsing binary', f'odd/even depth differ at {oe:.1f} sigma'
    elif rstar is None or not depth_ppm or depth_ppm <= 0:
        return 'unclassified', 'no stellar radius in TIC v8.2' if rstar is None else 'no transit depth'
    else:
        rp = np.sqrt(depth_ppm / 1e6) * rstar / RJ_RSUN
        if rp < PLANET_MAX_RJ:
            t = 'planet'
        elif rp <= BD_MAX_RJ:
            t = 'brown dwarf or low-mass star'
        else:
            t = 'eclipsing binary'
        why = f'Rp {rp:.2f} R_J from R_star {rstar:.2f} R_sun'
    if giant:
        why += f'; giant host (R_star {rstar:.1f} R_sun)'
    return t, why


GAIA_RADIUS_DEG = 2 / 3600          # brightest Gaia DR3 source within 2 arcsec is the host
RUWE_MAX = 1.4                      # above this the single-star astrometric fit is poor (unresolved companion)
P_TOL = 0.02                        # Gaia orbit period agrees with the transit period (or 1/2, 2x) to 2 per cent


def positions(con, tics):
    """{tic: (ra, dec)} of the hosts, from the tic table (TESS) or the star table (Kepler, K2)."""
    pos = {}
    for t in tics:
        r = con.execute('SELECT ra, dec FROM tic WHERE tic=?' if sources.mission_of(t) == 'TESS'
                        else 'SELECT ra, dec FROM star WHERE sid=?', (t,)).fetchone()
        if r and r[0] is not None:
            pos[t] = r
    return pos


def gaia_dr3(con, tics):
    """{tic: dict} of Gaia DR3 RUWE, non_single_star, any two-body orbit or acceleration solution, and the radial-velocity
    variability statistics of each host, from one asynchronous upload cross-match at the ESA archive. Returns {} when
    the archive cannot be reached."""
    import tempfile
    from pathlib import Path
    from astropy.table import Table
    from astroquery.gaia import Gaia
    pos = positions(con, tics)
    if not pos:
        return {}
    f = Path(tempfile.mkdtemp()) / 'pos.xml'
    Table(rows=[(k, *v) for k, v in pos.items()], names=('sid', 'ra', 'dec'),
          dtype=('i8', 'f8', 'f8')).write(f, format='votable')
    q = f"""SELECT u.sid, g.phot_g_mean_mag AS gmag, g.ruwe, g.non_single_star, g.rv_amplitude_robust AS rvamp,
        g.rv_chisq_pvalue AS rvp, g.rv_renormalised_gof AS rvgof, g.rv_nb_transits AS rvn,
        o.nss_solution_type AS orb, o.period, o.eccentricity, o.semi_amplitude_primary AS k1,
        a.nss_solution_type AS acc, e.frequency AS ebfreq
        FROM tap_upload.pos AS u JOIN gaiadr3.gaia_source AS g
        ON 1=CONTAINS(POINT(u.ra, u.dec), CIRCLE(g.ra, g.dec, {GAIA_RADIUS_DEG}))
        LEFT JOIN gaiadr3.nss_two_body_orbit AS o ON o.source_id=g.source_id
        LEFT JOIN gaiadr3.nss_acceleration_astro AS a ON a.source_id=g.source_id
        LEFT JOIN gaiadr3.vari_eclipsing_binary AS e ON e.source_id=g.source_id"""
    try:
        r = Gaia.launch_job_async(q, upload_resource=str(f), upload_table_name='pos').get_results()
    except Exception as e:                       # archive down or overloaded: classify without Gaia
        print(f'Gaia DR3 cross-match skipped: {e}')
        return {}
    num = lambda v: None if np.ma.is_masked(v) else float(v)
    out = {}
    for row in sorted(r, key=lambda x: (int(x['sid']), float(x['gmag']) if not np.ma.is_masked(x['gmag']) else 99)):
        out.setdefault(int(row['sid']), {'ruwe': num(row['ruwe']), 'nss': int(row['non_single_star'] or 0),
                                         'orb': None if np.ma.is_masked(row['orb']) else str(row['orb']),
                                         'period': num(row['period']), 'ecc': num(row['eccentricity']),
                                         'k1': num(row['k1']),
                                         'rvamp': num(row['rvamp']), 'rvp': num(row['rvp']),
                                         'rvgof': num(row['rvgof']), 'rvn': num(row['rvn']),
                                         'ebper': 1 / num(row['ebfreq']) if num(row['ebfreq']) else None,
                                         'acc': None if np.ma.is_masked(row['acc']) else str(row['acc'])})
    return out


def companion_mass(p_day, k1_kms, ecc, m1_sun):
    """Companion mass (M_sun) from the spectroscopic mass function with sin i = 1, which holds for a transiting or
    eclipsing system: M2^3 / (M1 + M2)^2 = P K^3 (1 - e^2)^1.5 / (2 pi G)."""
    if not (p_day and k1_kms and m1_sun):
        return None
    fm = 1.0361e-7 * (1 - (ecc or 0) ** 2) ** 1.5 * k1_kms ** 3 * p_day     # mass function in M_sun
    m2 = fm
    for _ in range(100):
        m2 = (fm * (m1_sun + m2) ** 2) ** (1 / 3)
    return m2


def gaia_test(g, period, m1):
    """(is_binary, note) from the Gaia DR3 entry of one host: an orbit at the transit period (or half or twice it)
    means a stellar-mass companion that Gaia could detect, so the signal is an eclipsing binary."""
    if not g:
        return False, ''
    notes = []
    if g['orb'] and g['period'] and period:
        ratio = g['period'] / period
        hit = next((m for m in (1, 0.5, 2) if abs(ratio / m - 1) < P_TOL), None)
        if hit:
            m2 = companion_mass(g['period'], g['k1'], g['ecc'], m1)
            mtxt = f', M2 {m2:.2f} M_sun' if m2 else ''
            return True, f"Gaia DR3 {g['orb']} orbit P {g['period']:.3f} d ({hit}x transit period){mtxt}"
        notes.append(f"Gaia DR3 {g['orb']} orbit at P {g['period']:.2f} d (not the transit period)")
    if g.get('ebper') and period:
        ratio = g['ebper'] / period
        hit = next((m for m in (1, 0.5, 2) if abs(ratio / m - 1) < P_TOL), None)
        if hit:
            return True, f"Gaia DR3 eclipsing binary in its own photometry, P {g['ebper']:.3f} d ({hit}x transit period)"
        notes.append(f"Gaia DR3 photometric eclipsing binary at P {g['ebper']:.2f} d")
    if g['acc']:
        notes.append(f"Gaia DR3 astrometric {g['acc']}")
    if g['ruwe'] and g['ruwe'] > RUWE_MAX:
        notes.append(f"Gaia RUWE {g['ruwe']:.2f} (possible unresolved companion)")
    return False, '; '.join(notes)


MSTAR_MIN = 0.08                    # companion above the hydrogen-burning limit (M_sun) is a star
RV_P, RV_GOF, RV_NMIN = 0.01, 4.0, 10   # Gaia DR3 RV-variability criteria (Katz et al. 2023)
APOGEE_VSCAT = 1.0                  # APOGEE visit-to-visit RV scatter (km/s) above which the host is RV variable


def rv_test(k_kms, period, m1, src):
    """(is_binary, note) from an RV semi-amplitude estimate at the transit period: a companion above 0.08 M_sun is a
    star, so the transit is an eclipse. A planet induces K of order 0.1 km/s or less."""
    m2 = companion_mass(period, k_kms, 0, m1)
    if m2 is None:
        return False, ''
    note = f'{src} RV K ~ {k_kms:.1f} km/s, M2 ~ {m2:.2f} M_sun at the transit period'
    return m2 > MSTAR_MIN, note


def gaia_rv(g, period, m1):
    """Gaia DR3 RVS variability: K ~ half the robust peak-to-peak amplitude, only for a host that passes the
    constant-RV rejection (p < 0.01, renormalised GoF > 4, at least 10 transits)."""
    if not g or not g.get('rvamp') or g.get('rvp') is None or g.get('rvgof') is None:
        return False, ''
    if g['rvp'] < RV_P and g['rvgof'] > RV_GOF and (g['rvn'] or 0) >= RV_NMIN:
        return rv_test(g['rvamp'] / 2, period, m1, 'Gaia DR3')
    return False, ''


def apogee(pos):
    """{tic: (vscatter, nvisits)} from the SDSS DR17 APOGEE allStar table (nearest entry within 2 arcsec), one query
    per host through the SkyServer SQL service."""
    import requests
    from concurrent.futures import ThreadPoolExecutor
    url = 'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch'
    d = 2 / 3600

    def one(item):
        t, (ra, dec) = item
        c = max(np.cos(np.radians(dec)), 1e-3)
        q = (f'select top 1 vscatter, nvisits from apogeeStar where ra between {ra - d / c} and {ra + d / c} '
             f'and dec between {dec - d} and {dec + d} and nvisits > 1 order by nvisits desc')
        try:
            rows = requests.get(url, params={'cmd': q, 'format': 'json'}, timeout=60).json()[0]['Rows']
        except Exception:
            return t, None
        return t, (rows[0]['vscatter'], rows[0]['nvisits']) if rows else None

    with ThreadPoolExecutor(8) as ex:
        return {t: v for t, v in ex.map(one, pos.items()) if v}


def apogee_test(a, period, m1):
    """APOGEE multi-epoch RVs: for a sinusoidal orbit the visit RV standard deviation is K / sqrt(2)."""
    if not a or a[0] is None or a[0] < APOGEE_VSCAT:
        return False, ''
    return rv_test(np.sqrt(2) * a[0], period, m1, f'APOGEE ({a[1]} visits)')


LAMOST_FLOOR = 5.0                  # LAMOST low-resolution RV systematic floor (km/s), added in quadrature to e_HRV
LAMOST_P = 0.01                     # chi-square probability below which the LAMOST epochs are not constant


def lamost(pos):
    """{tic: [(hrv, e_hrv), ...]} from LAMOST DR5 low-resolution A/F/G/K spectra (V/164/stellar5) within 3 arcsec,
    read from the CfA VizieR mirror (the CDS server is often unreachable)."""
    import requests
    from concurrent.futures import ThreadPoolExecutor
    url = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'

    def one(item):
        t, (ra, dec) = item
        try:
            txt = requests.get(url, params={'-source': 'V/164/stellar5', '-c': f'{ra} {dec}', '-c.rs': 3,
                                            '-out': 'HRV,e_HRV', '-out.max': 200}, timeout=60).text
        except Exception:
            return t, None
        ep = []
        for line in txt.splitlines():
            p = line.split('\t')
            if len(p) == 2 and not line.startswith('#'):
                try:
                    ep.append((float(p[0]), float(p[1])))
                except ValueError:
                    pass
        return t, ep if len(ep) > 1 else None

    with ThreadPoolExecutor(8) as ex:
        return {t: v for t, v in ex.map(one, pos.items()) if v}


def lamost_test(ep, period, m1):
    """LAMOST multi-epoch RVs: variable when the epochs reject a constant RV (chi-square p < 0.01 with the 5 km/s
    floor); K is then half the peak-to-peak spread."""
    if not ep:
        return False, ''
    from scipy.stats import chi2
    v = np.array([e[0] for e in ep])
    s = np.hypot([e[1] for e in ep], LAMOST_FLOOR)
    w = 1 / s ** 2
    x2 = np.sum(w * (v - np.sum(w * v) / np.sum(w)) ** 2)
    if chi2.sf(x2, len(v) - 1) >= LAMOST_P:
        return False, ''
    return rv_test(np.ptp(v) / 2, period, m1, f'LAMOST ({len(v)} epochs)')


ZTF_MIN_PPM = 3000                  # ZTF per-point precision (~1 per cent) only tests eclipses deeper than this
ZTF_BRIGHT = 13.0                   # ZTF saturates for brighter hosts
ZTF_DILUTION = 1.5                  # ZTF depth this much deeper than the TESS/Kepler depth replaces it for the radius
BJD_OFFSET = {'TESS': 2457000.0, 'Kepler': 2454833.0, 'K2': 2454833.0}


def ztf_lc(ra, dec):
    """ZTF DR light curve (all bands, clean epochs only) from the IRSA light-curve service: (hjd, rel. flux, err)."""
    import io
    import pandas as pd
    import requests
    r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves',
                     params={'POS': f'CIRCLE {ra} {dec} 0.0005', 'FORMAT': 'csv'}, timeout=300)
    df = pd.read_csv(io.StringIO(r.text))
    df = df[(df.catflags == 0) & (df.mag > ZTF_BRIGHT)]
    t, f, e = [], [], []
    for _, g in df.groupby(['oid']):           # each field/CCD/band object normalised by its own median
        if len(g) < 20:
            continue
        fl = 10 ** (-0.4 * (g.mag.values - np.median(g.mag.values)))
        t.append(g.hjd.values); f.append(fl); e.append(0.921 * g.magerr.values * fl)
    if not t:
        return None
    return np.concatenate(t), np.concatenate(f), np.concatenate(e)


def ztf_test(lc, mission, period, t0, dur_h, depth_ppm):
    """(is_binary, note) from ZTF. A box search over a narrow period window (+-0.05 per cent, enough to absorb the
    ephemeris drift over the ZTF baseline) recovers the eclipse; then the odd and even depths and the depth at phase
    0.5 are compared as in planets.vet."""
    from astropy.timeseries import BoxLeastSquares
    if lc is None or not period or not dur_h:
        return False, '', 0
    t, f, e = lc
    if len(t) < 100:
        return False, '', 0
    dur = dur_h / 24
    span = np.ptp(t)
    step = 0.25 * dur * period / span
    periods = np.arange(period * (1 - 5e-4), period * (1 + 5e-4), step)
    if len(periods) > 200000:
        return False, '', 0
    bls = BoxLeastSquares(t, f, e)
    res = bls.power(periods, dur, objective='snr')
    i = int(np.argmax(res.depth_snr))
    if res.depth_snr[i] < 7:
        return False, f'ZTF: eclipse not recovered ({len(t)} epochs)', 0
    st = bls.compute_stats(res.period[i], dur, res.transit_time[i])
    (do, eo), (de, ee) = st['depth_odd'], st['depth_even']
    oe = abs(do - de) / np.hypot(eo, ee) if eo > 0 and ee > 0 else 0
    sd, se = st['depth_phased']
    sec = sd / se if se > 0 else 0
    zd = res.depth[i] * 1e6
    note = f'ZTF depth {zd:.0f} ppm (TESS/Kepler {depth_ppm:.0f}), odd/even {oe:.1f} sigma, secondary {sec:.1f} sigma'
    if oe > OE_SIGMA:
        return True, f'ZTF odd/even depths differ at {oe:.1f} sigma; ' + note, zd
    if sec > SEC_SIGMA:
        return True, f'ZTF secondary eclipse at {sec:.1f} sigma; ' + note, zd
    return False, note, zd


def ztf_all(pos, cands):
    """{tic: lc} for the deep candidates north of dec -30, fetched in parallel."""
    from concurrent.futures import ThreadPoolExecutor
    todo = [t for t, dep in cands if dep and dep >= ZTF_MIN_PPM and t in pos and pos[t][1] > -30]

    def one(t):
        try:
            return t, ztf_lc(*pos[t])
        except Exception:
            return t, None

    with ThreadPoolExecutor(8) as ex:
        return {t: lc for t, lc in ex.map(one, todo) if lc is not None}


def stellar_mass(con, tic):
    r = con.execute('SELECT mass FROM tic WHERE tic=?' if sources.mission_of(tic) == 'TESS'
                    else 'SELECT mass FROM star WHERE sid=?', (tic,)).fetchone()
    return r[0] if r and r[0] and r[0] > 0 else None


def run(status='ready', tics=None, gaia=True, external=True):
    con = survey.db()
    if 'objtype' not in [r[1] for r in con.execute('PRAGMA table_info(candidate)')]:
        con.execute('ALTER TABLE candidate ADD COLUMN objtype TEXT')
    rows = con.execute('SELECT tic, depth_ppm, notes, period, t0, dur_h FROM candidate WHERE status=?',
                       (status,)).fetchall()
    if tics:
        rows = [r for r in rows if r[0] in set(tics)]
    fill_tic(con, [r[0] for r in rows if stellar_radius(con, r[0]) is None])
    pos = positions(con, [r[0] for r in rows])
    gd = gaia_dr3(con, [r[0] for r in rows]) if gaia else {}
    ap = apogee(pos) if external else {}
    lm = lamost(pos) if external else {}
    zl = ztf_all(pos, [(r[0], r[1]) for r in rows]) if external else {}
    print(f'Gaia {len(gd)}, APOGEE {len(ap)}, LAMOST {len(lm)}, ZTF {len(zl)} hosts matched')
    counts = {}
    for tic, dep, notes, period, t0, dur in rows:
        zeb, znote, zd = ztf_test(zl.get(tic), sources.mission_of(tic), period, t0, dur, dep or 0)
        if dep and zd > ZTF_DILUTION * dep:      # ZTF's 1 arcsec seeing resolves contaminants that TESS pixels blend
            t, why = label(zd, stellar_radius(con, tic), *vet_tests(tic))
            why += f' (ZTF undiluted depth; TESS/Kepler depth {dep:.0f} ppm)'
        else:
            t, why = label(dep, stellar_radius(con, tic), *vet_tests(tic))
        m1 = stellar_mass(con, tic)
        tests = [gaia_test(gd.get(tic), period, m1), gaia_rv(gd.get(tic), period, m1),
                 apogee_test(ap.get(tic), period, m1), lamost_test(lm.get(tic), period, m1), (zeb, znote)]
        for eb, note in tests:
            if eb:
                t, why = 'eclipsing binary', note + '; ' + why
            elif note:
                why += '; ' + note
        kept = [n for n in (notes or '').split('; ') if n and not n.startswith('class: ')]
        con.execute('UPDATE candidate SET objtype=?, notes=?, updated=? WHERE tic=?',
                    (t, '; '.join(kept + [f'class: {t} ({why})']), time.strftime('%Y-%m-%d %H:%M'), tic))
        counts.setdefault(t, []).append(tic)
    con.commit()
    for t, v in sorted(counts.items()):
        print(f'{t:30s} {len(v)}' + (': ' + ' '.join(map(str, v)) if t == 'unclassified' else ''))
    return counts


def main(argv):
    ap = argparse.ArgumentParser(prog='occulta classify', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--status', default='ready', help='candidates of this status (default: ready)')
    ap.add_argument('--no-gaia', action='store_true', help='skip the Gaia DR3 cross-match')
    ap.add_argument('--no-external', action='store_true', help='skip the APOGEE and ZTF checks')
    ap.add_argument('tic', type=int, nargs='*')
    a = ap.parse_args(argv)
    run(a.status, a.tic, gaia=not a.no_gaia, external=not a.no_external)
