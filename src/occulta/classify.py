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
replaces it in the radius. The observed duration is compared with the longest central (b = 0) transit that the host
density (from the TIC radius and mass) allows at the period, including the companion radius; a planet or brown-dwarf
candidate whose duration exceeds DUR_FACTOR (1.5) times that maximum becomes unclassified, since its size rests on a host
that cannot be the eclipsed star. Without a density the test is skipped and noted. Three further checks then try to
name the eclipsed object. A host with no TIC radius or mass takes the Gaia DR3 FLAME (else GSP-Phot) values. A dip
that lasts more than DUTY_CONTACT (0.25) of the period needs (R1 + R2) / a > 0.7 whatever the host, which only a
contact or near-contact pair gives, so it is an eclipsing binary. A candidate that fails the duration test and has a
Gaia neighbour closer than 8 arcsec (which the difference-image test of the vetting cannot separate) that could
produce the depth with an eclipse of at most half its light is a blended eclipsing binary on that neighbour. What is still
unclassified then goes through a transit-shape fit (shape_test): the folded transit windows of the vetting are fitted
with a quadratic limb-darkened transit (batman; Claret 2017 and Claret & Bloemen 2011 coefficients) for Rp/R_star,
a/R_star, b and the mid-time, with bootstrap errors, which gives the density of the eclipsed star, rho_LC. For a host
without a radius, rho_LC must agree with a dwarf of the host's Teff (TIC, Gaia GSP-Phot, or Gaia Bp-Rp on the Pecaut &
Mamajek 2013 sequence), and the dwarf radius at rho_LC then gives Rp. For a signal too long for its host, a grazing
V-shaped fit, or a rho_LC below what any orbit with e < 0.6 allows around the host, makes it an eclipsing binary, and a
flat-bottomed fit consistent with the catalogue density gives Rp = k R_star. Every decision needs 95 per cent of the
bootstrap fits on its side. What is left unclassified has neither a stellar radius nor an identified eclipsed star, nor
a transit shape that decides between them. ZTF depths and odd/even tests need at least three in-eclipse epochs
whose median shows the eclipse. The label goes into the new candidate column objtype and the reason into notes.

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


def label(depth_ppm, rstar, oe, sec, src='TIC v8.2'):
    """(objtype, reason) of one candidate; src names the catalogue of the stellar radius."""
    giant = rstar is not None and rstar > GIANT_RSUN
    if np.isfinite(sec) and sec > SEC_SIGMA:
        t, why = 'eclipsing binary', f'secondary eclipse at {sec:.1f} sigma'
    elif np.isfinite(oe) and oe > OE_SIGMA:
        t, why = 'eclipsing binary', f'odd/even depth differ at {oe:.1f} sigma'
    elif rstar is None or not depth_ppm or depth_ppm <= 0:
        return 'unclassified', 'no stellar radius in TIC v8.2 or Gaia DR3' if rstar is None else 'no transit depth'
    else:
        rp = np.sqrt(depth_ppm / 1e6) * rstar / RJ_RSUN
        if rp < PLANET_MAX_RJ:
            t = 'planet'
        elif rp <= BD_MAX_RJ:
            t = 'brown dwarf or low-mass star'
        else:
            t = 'eclipsing binary'
        why = f'Rp {rp:.2f} R_J from R_star {rstar:.2f} R_sun' + ('' if src == 'TIC v8.2' else f' ({src})')
    if giant:
        why += f'; giant host (R_star {rstar:.1f} R_sun)'
    return t, why


DUR_FACTOR = 1.5                    # observed duration above this multiple of the longest central transit fails
G_CGS = 6.674e-8                    # gravitational constant (cm^3 g^-1 s^-2)


def stellar_density(rstar, mstar):
    """Mean stellar density (g/cm^3) from R_star and M_star in solar units, or None when either is unknown."""
    if not rstar or not mstar or rstar <= 0 or mstar <= 0:
        return None
    return 1.40980 * mstar / rstar ** 3      # mean solar density 1.4098 g/cm^3


def max_duration_h(period_d, rho_cgs, depth_ppm=0.0):
    """Longest first-to-fourth-contact duration (h) of a central (b = 0) transit on a circular orbit (Seager &
    Mallen-Ornelas 2003; Winn 2010, eq. 14): T = (P / pi) arcsin[(1 + k) / (a / R_star)], with k = sqrt(depth) and
    a / R_star = (G rho_star P^2 / 3 pi)^(1/3) from Kepler's third law (companion mass neglected)."""
    if not period_d or not rho_cgs or period_d <= 0 or rho_cgs <= 0:
        return None
    p_s = period_d * 86400.0
    a_r = (G_CGS * rho_cgs * p_s ** 2 / (3 * np.pi)) ** (1 / 3)
    k = np.sqrt(depth_ppm / 1e6) if depth_ppm and depth_ppm > 0 else 0.0
    x = min(1.0, (1 + k) / a_r)              # contact orbit: the transit cannot last longer than half the period
    return period_d / np.pi * np.arcsin(x) * 24


def duration_test(dur_h, period_d, rho_cgs, depth_ppm=0.0, factor=DUR_FACTOR):
    """(too_long, note): the observed duration against the longest central transit that the host density allows. A
    duration above factor x that maximum is unlikely to come from a companion of this star (only e > ~0.4 near
    apoastron stretches it by 1.5), so the host parameters, or the star that is eclipsed, are probably not the
    catalogued ones. No density: no test."""
    if rho_cgs is None:
        return False, 'duration test: no density'
    tmax = max_duration_h(period_d, rho_cgs, depth_ppm)
    if tmax is None or not dur_h or dur_h <= 0:
        return False, 'duration test: no period or duration'
    r = dur_h / tmax
    note = f'duration {dur_h:.2f} h is {r:.2f} x the longest central transit ({tmax:.2f} h, rho_star {rho_cgs:.2f} g/cc)'
    return r > factor, note


DUTY_CONTACT = 0.25                 # dip longer than this fraction of the period: (R1 + R2) / a > sin(pi / 4)
BLEND_SEP_ARCSEC = 8.0              # closer neighbours are not separated by the difference-image test of the vetting
BLEND_MAX_ECLIPSE = 0.5             # a neighbour must dim by at most half its light to give the observed depth


def contact_test(dur_h, period_d, frac=DUTY_CONTACT):
    """(is_binary, note): a dip that lasts more than frac of the period. For a central eclipse on a circular orbit
    T / P = arcsin((R1 + R2) / a) / pi, so T / P > 0.25 needs (R1 + R2) / a > 0.71 for any host: the two bodies
    (nearly) touch, as in a W UMa or near-contact binary, which a planet cannot do outside its Roche limit."""
    if not dur_h or not period_d or period_d <= 0:
        return False, ''
    q = dur_h / 24 / period_d
    if q <= frac:
        return False, ''
    return True, (f'contact or near-contact binary: the dip lasts {q:.2f} of the period, (R1 + R2) / a > '
                  f'{np.sin(np.pi * min(q, 0.5)):.2f}')


def vet_blends(tic):
    """Resolved Gaia DR3 neighbours that could produce the depth, from the vetting summary ([] when unknown)."""
    p = survey.FIND / sources.name_of(tic).replace(' ', '') / 'summary.json'
    if not p.exists():
        return []
    g = json.loads(p.read_text()).get('gaia') or {}
    return g.get('possible_blends') or []


def blend_test(blends, sep_max=BLEND_SEP_ARCSEC, need_max=BLEND_MAX_ECLIPSE):
    """(is_blend, note): the nearest neighbour closer than sep_max whose eclipse of at most need_max of its own light
    would give the observed depth."""
    ok = [b for b in blends if b.get('sep_arcsec') is not None and b['sep_arcsec'] <= sep_max
          and b.get('eclipse_needed') is not None and b['eclipse_needed'] <= need_max]
    if not ok:
        return False, ''
    b = min(ok, key=lambda x: x['sep_arcsec'])
    return True, (f"blended eclipsing binary: Gaia DR3 neighbour at {b['sep_arcsec']} arcsec (dG {b['dG']:+.2f}) "
                  f"gives the depth with a {100 * b['eclipse_needed']:.1f} per cent eclipse")


CLASS_PARTS = ('Rp ', 'giant host', 'Gaia ', 'TESS/Kepler depth', 'ZTF', 'APOGEE', 'LAMOST', 'duration ',
               'too long', 'no stellar radius', 'no transit depth', 'secondary eclipse at', 'odd/even depth',
               'was ', 'contact or near-contact', 'blended eclipsing')


def strip_class(notes):
    """The vetting notes without the previous class note. The class reason itself contains '; ', so its pieces
    (and those left by older runs) are dropped by their opening words."""
    pre = (notes or '').split('class: ')[0]
    return [n for n in pre.split('; ') if n and not n.startswith(CLASS_PARTS)]


GAIA_RADIUS_DEG = 2 / 3600         # brightest Gaia DR3 source within 2 arcsec is the host
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
        a.nss_solution_type AS acc, e.frequency AS ebfreq,
        p.radius_flame AS rflame, p.radius_gspphot AS rgsp, p.mass_flame AS mflame,
        g.bp_rp, p.teff_gspphot AS tgsp, g.parallax AS plx, g.parallax_error AS plxe, p.ag_gspphot AS ag
        FROM tap_upload.pos AS u JOIN gaiadr3.gaia_source AS g
        ON 1=CONTAINS(POINT(u.ra, u.dec), CIRCLE(g.ra, g.dec, {GAIA_RADIUS_DEG}))
        LEFT JOIN gaiadr3.nss_two_body_orbit AS o ON o.source_id=g.source_id
        LEFT JOIN gaiadr3.nss_acceleration_astro AS a ON a.source_id=g.source_id
        LEFT JOIN gaiadr3.vari_eclipsing_binary AS e ON e.source_id=g.source_id
        LEFT JOIN gaiadr3.astrophysical_parameters AS p ON p.source_id=g.source_id"""
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
                                         'acc': None if np.ma.is_masked(row['acc']) else str(row['acc']),
                                         'rad': num(row['rflame']) or num(row['rgsp']), 'mass': num(row['mflame']),
                                         'bp_rp': num(row['bp_rp']), 'teff': num(row['tgsp']),
                                         'gmag': num(row['gmag']), 'plx': num(row['plx']), 'plx_e': num(row['plxe']),
                                         'ag': num(row['ag'])})
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
ZTF_MIN_IN = 3                      # in-eclipse ZTF epochs needed before a ZTF depth or odd/even test is trusted
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
    # the box mean can rest on one bad epoch; the median of the in-eclipse epochs must show the same eclipse
    ph = (t - res.transit_time[i] + 0.5 * res.period[i]) % res.period[i] - 0.5 * res.period[i]
    inn = np.abs(ph) < dur / 2
    zmed = 1 - np.median(f[inn]) / np.median(f[~inn]) if inn.sum() else 0.0
    if inn.sum() < ZTF_MIN_IN or zmed < 0.5 * res.depth[i]:
        return False, (f'ZTF: eclipse not recovered ({len(t)} epochs, {int(inn.sum())} in eclipse, median depth '
                       f'{zmed * 1e6:.0f} ppm)'), 0
    st = bls.compute_stats(res.period[i], dur, res.transit_time[i])
    (do, eo), (de, ee) = st['depth_odd'], st['depth_even']
    oe = abs(do - de) / np.hypot(eo, ee) if eo > 0 and ee > 0 else 0
    sd, se = st['depth_phased']
    sec = sd / se if se > 0 else 0
    zd = zmed * 1e6
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
        for _ in range(2):                       # a dropped IRSA request must not change the class between runs
            try:
                return t, ztf_lc(*pos[t])
            except Exception:
                time.sleep(5)
        return t, None

    with ThreadPoolExecutor(8) as ex:
        return {t: lc for t, lc in ex.map(one, todo) if lc is not None}


# ---- transit-shape fit: the class of the signals that the radius or the duration leaves undetermined ----------------
# The folded transit windows of the vetting are fitted with a quadratic limb-darkened transit (batman; Mandel & Agol
# 2002, Kreidberg 2015) on a circular orbit with the period fixed. Free: k = Rp/R*, a/R*, b, the mid-time and the
# flux zero point. The density of the eclipsed star follows from Kepler's third law, rho_LC = 3 pi (a/R*)^3 / (G P^2)
# (Seager & Mallen-Ornelas 2003). Errors come from N_BOOT bootstrap resamplings of the points, and every decision
# needs CONF of the resampled fits on its side, so that an ambiguous fit leaves the signal undetermined.
FIT_SNR_MIN = 10.0                  # folded vetting SNR below this: no shape fit
FIT_BINS, FIT_COVER_MIN = 10, 0.8   # at least 80 per cent of ten bins across the transit must hold data
FIT_OOT_MIN = 10                    # points needed before and after the transit
N_BOOT = 100
CONF = 0.95                         # fraction of the bootstrap fits that a decision needs
GRAZE_B = 0.9                       # grazing: b > 0.9 and b + k > 1 (V shape, no second and third contact)
ECC_MAX = 0.6                       # largest eccentricity allowed for a long-period companion
RHO_CAT_FAC = 1.3                   # catalogue density uncertainty (TIC radius ~7 per cent, mass ~15 per cent)
# a transit at apoastron (omega = -90 deg) lasts longer by (1 + e)/(1 - e))^(1/2), so the circular fit gives
# rho_LC / rho_star = ((1 - e) / (1 + e))^(3/2) (Kipping 2010; Dawson & Johnson 2012): 0.125 at e = 0.6, and 0.096
# once the catalogue density may be 1.3 times too high
RHO_EB_RATIO = ((1 - ECC_MAX) / (1 + ECC_MAX)) ** 1.5 / RHO_CAT_FAC
MS_RHO_FAC = 1.5                    # spread of the main-sequence density at fixed Teff (age, metallicity)
CADENCES_S = (20.0, 58.85, 120.0, 200.0, 600.0, 1765.5, 1800.0)


def _tables():
    from importlib import resources
    if not hasattr(_tables, 'd'):
        _tables.d = json.loads(resources.files('occulta').joinpath('stellar_tables.json').read_text())
    return _tables.d


def ld_coeffs(teff, logg, band='TESS'):
    """Quadratic limb-darkening (u1, u2) of the nearest grid point (Teff in 250-K, log g in 0.5-dex units): Claret
    (2017) for TESS, Claret & Bloemen (2011) for Kepler and K2; ATLAS above 3500 K and PHOENIX below, solar
    metallicity. Unknown Teff or log g: 5750 K, 4.5."""
    g = np.array(_tables()['claret']['TESS' if band == 'TESS' else 'Kepler'])
    teff = teff if teff and np.isfinite(teff) else 5750.0
    logg = logg if logg and np.isfinite(logg) else 4.5
    i = int(np.argmin(((g[:, 0] - teff) / 250) ** 2 + ((g[:, 1] - logg) / 0.5) ** 2))
    return float(g[i, 2]), float(g[i, 3])


def ms_dwarfs():
    """(Teff, Bp-Rp, R, M, rho) of the mean dwarf sequence of Pecaut & Mamajek (2013), rho in g/cm^3."""
    m = np.array(_tables()['mamajek'])
    return m[:, 0], m[:, 1], m[:, 2], m[:, 3], 1.40980 * m[:, 3] / m[:, 2] ** 3


def ms_teff(bp_rp):
    """Dwarf Teff at a Gaia Bp-Rp colour (no reddening correction), None outside the sequence."""
    teff, c = ms_dwarfs()[:2]
    if bp_rp is None or not np.isfinite(bp_rp) or not c.min() <= bp_rp <= c.max():
        return None
    o = np.argsort(c)
    return float(np.interp(bp_rp, c[o], teff[o]))


def ms_rho(teff):
    """Mean density (g/cm^3) of a dwarf of this Teff."""
    t, _, _, _, rho = ms_dwarfs()
    return float(np.exp(np.interp(teff, t, np.log(rho))))


def ms_radius(rho):
    """Radius (R_sun) of the dwarf of this mean density, from the sequence ordered by density; NaN off the sequence."""
    _, _, r, _, d = ms_dwarfs()
    o = np.argsort(d)
    x = np.log(np.asarray(rho, float))
    return np.where((x >= np.log(d.min())) & (x <= np.log(d.max())), np.interp(x, np.log(d[o]), r[o]), np.nan)


def rho_from_ar(ar, period_d):
    """rho_LC = 3 pi (a/R*)^3 / (G P^2) in g/cm^3."""
    return 3 * np.pi * np.asarray(ar) ** 3 / (G_CGS * (period_d * 86400.0) ** 2)


def ar_from_rho(rho, period_d):
    return (G_CGS * rho * (period_d * 86400.0) ** 2 / (3 * np.pi)) ** (1 / 3)


def _models(tf, cad_s, period, u, fac=None):
    """One batman model per cadence, integrated over the exposure (supersampled to about 4 min). fac ({cadence:
    step factor}) reuses the integration step of an earlier model, which batman otherwise tunes on every call."""
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = 0.0, period, 0.1, 20.0, 90.0, 0.0, 90.0, list(u), 'quadratic'
    out = []
    for c in np.unique(cad_s):
        idx = np.flatnonzero(cad_s == c)
        ss = max(1, int(round(c / 240.0)))
        m = batman.TransitModel(p, tf[idx], supersample_factor=ss, exp_time=c / 86400.0,
                                fac=None if fac is None else fac.get(c))
        out.append((idx, m, p, c))
    return out


def _cadence(t):
    """Exposure (s) of each point: the shorter gap to a neighbour, snapped to the nearest mission cadence."""
    d = np.diff(t) * 86400.0
    g = np.minimum(np.r_[d[:1], d], np.r_[d, d[-1:]]) if len(d) else np.array([1800.0])
    c = np.array(CADENCES_S)
    return c[np.argmin(np.abs(np.log(np.maximum(g, 1.0))[:, None] - np.log(c)[None, :]), axis=1)]


def transit_model(theta, models, n, period):
    """Flux of the folded light curve for theta = (k, log10 a/R*, b, dt0, f0)."""
    k, lar, b, dt0, f0 = theta
    ar = 10 ** lar
    y = np.empty(n)
    for idx, m, p, _ in models:
        p.rp, p.a, p.t0 = k, ar, dt0
        p.inc = float(np.degrees(np.arccos(min(b / ar, 1.0))))
        y[idx] = m.light_curve(p)
    return y - 1 + f0


def fit_coverage(tf, dur_h):
    """(ok, note): FIT_COVER_MIN of FIT_BINS bins across the transit hold data, and FIT_OOT_MIN points lie on
    each side of it."""
    h = dur_h / 48.0
    inn = tf[np.abs(tf) < h]
    filled = len(np.unique(np.floor((inn + h) / (2 * h) * FIT_BINS).clip(0, FIT_BINS - 1))) / FIT_BINS
    pre, post = np.sum(tf < -1.2 * h), np.sum(tf > 1.2 * h)
    ok = filled >= FIT_COVER_MIN and min(pre, post) >= FIT_OOT_MIN
    return ok, f'{100 * filled:.0f} per cent of the transit covered, {pre}/{post} points before/after'


def transit_fit(t, f, period, t0, dur_h, depth_ppm, u=(0.4, 0.25), n_boot=N_BOOT, seed=0):
    """Least-squares fit of a limb-darkened transit to the folded points (t, f = relative flux - 1), with bootstrap
    errors. Returns None when the transit is not covered, else a dict with the best fit and the bootstrap samples of
    k, a/R*, b, rho_LC (g/cm^3) and T14 (h)."""
    from scipy.optimize import least_squares
    t, f = np.asarray(t, float), np.asarray(f, float)
    o = np.argsort(t)
    t, f = t[o], f[o]
    cad = _cadence(t)
    tf = (t - t0 + 0.5 * period) % period - 0.5 * period
    ok, cnote = fit_coverage(tf, dur_h)
    if not ok:
        return {'ok': False, 'note': cnote}
    out = np.abs(tf) > 0.75 * dur_h / 24
    sig = np.empty(len(f))
    for c in np.unique(cad):                     # white noise of each cadence from its out-of-transit scatter
        m = cad == c
        ref = f[m & out] if (m & out).sum() >= 10 else f[m]
        sig[m] = max(1.4826 * np.median(np.abs(ref - np.median(ref))), 1e-6)
    models = _models(tf, cad, period, u)
    fac = {c: m.fac for _, m, _, c in models}
    T = dur_h / 24
    k0 = np.sqrt(max(depth_ppm, 10.0) / 1e6)
    lo = [1e-3, np.log10(1.05), 0.0, -0.25 * T, -0.05]
    hi = [1.0, 3.0, 1.2, 0.25 * T, 0.05]

    best = None
    for b0 in (0.1, 0.5, 0.8, 0.95):
        kk = min(0.99, k0 * (1.3 if b0 > 0.9 else 1.0))
        ar0 = np.sqrt(max((1 + kk) ** 2 - b0 ** 2, 0.05)) / np.sin(min(np.pi * T / period, 1.5))
        x0 = np.clip([kk, np.log10(max(ar0, 1.1)), b0, 0.0, 0.0], np.add(lo, 1e-9), np.subtract(hi, 1e-9))
        x = _boot_one(f, sig, models, x0, lo, hi, T, period)
        chi = float(np.sum(((f - transit_model(x, models, len(f), period)) / sig) ** 2))
        if best is None or chi < best[1]:
            best = (x, chi)
    # each resampling is fitted twice, from the best fit and from a random impact parameter along the valley of equal
    # duration (b and a/R_star trade off at long cadence and low SNR), and keeps the better fit, so that the bootstrap
    # spread includes that degeneracy instead of staying at the first minimum
    rng = np.random.default_rng(seed)
    xb = best[0]
    t14b = float(t14_h(xb[0], 10 ** xb[1], xb[2], period)) / 24
    samp = []
    for _ in range(n_boot):
        idx = np.sort(rng.integers(0, len(f), len(f)))
        mods = _models(tf[idx], cad[idx], period, u, fac)
        b1 = rng.uniform(0, min(1.0 + xb[0], hi[2]) - 1e-3)
        ar1 = np.sqrt(max((1 + xb[0]) ** 2 - b1 ** 2, 1e-3)) / np.sin(min(np.pi * t14b / period, 1.5))
        x1 = np.clip([xb[0], np.log10(max(ar1, 1.06)), b1, xb[3], xb[4]], np.add(lo, 1e-9), np.subtract(hi, 1e-9))
        fits = [_boot_one(f[idx], sig[idx], mods, x0, lo, hi, T, period) for x0 in (xb, x1)]
        chi = [np.sum(((f[idx] - transit_model(x, mods, len(idx), period)) / sig[idx]) ** 2) for x in fits]
        samp.append(fits[int(np.argmin(chi))])
    s = np.array(samp)
    k, ar, b = s[:, 0], 10 ** s[:, 1], s[:, 2]
    return {'ok': True, 'note': cnote, 'k': k, 'ar': ar, 'b': b, 'rho': rho_from_ar(ar, period),
            't14': t14_h(k, ar, b, period), 'best': {'k': xb[0], 'ar': 10 ** xb[1], 'b': xb[2],
                                                     'rho': float(rho_from_ar(10 ** xb[1], period)),
                                                     'chi2': best[1], 'n': len(f)}}


def _boot_one(f, sig, models, x0, lo, hi, T, period):
    from scipy.optimize import least_squares
    r = least_squares(lambda th: (f - transit_model(th, models, len(f), period)) / sig, x0, bounds=(lo, hi),
                      x_scale=[0.05, 0.1, 0.1, 0.01 * T + 1e-3, 1e-4], ftol=1e-6, xtol=1e-6)
    return r.x


def t14_h(k, ar, b, period):
    """First-to-fourth-contact duration (h) on a circular orbit (Winn 2010, eq. 14)."""
    k, ar, b = np.asarray(k), np.asarray(ar), np.asarray(b)
    inc = np.arccos(np.clip(b / ar, 0, 1))
    x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (ar * np.sin(inc))
    return period / np.pi * np.arcsin(np.clip(x, 0, 1)) * 24


def class_from_rp(rp_rj):
    """Class when at least CONF of the companion radii (R_J) fall in one class interval, else None."""
    rp = np.asarray(rp_rj, float)
    rp = rp[np.isfinite(rp)]
    if not len(rp):
        return None
    for t, m in (('planet', rp < PLANET_MAX_RJ), ('brown dwarf or low-mass star', (rp >= PLANET_MAX_RJ) & (rp <= BD_MAX_RJ)),
                 ('eclipsing binary', rp > BD_MAX_RJ)):
        if np.mean(m) >= CONF and len(rp) == len(np.asarray(rp_rj)):
            return t
    return None


def _ci(x):
    lo, med, hi = np.percentile(x, [2.5, 50, 97.5])
    return f'{med:.3g} ({lo:.3g}-{hi:.3g})'


def shape_long(fit, rho_cat, rstar):
    """(class or None, note) for a signal too long for its catalogued host. Grazing in CONF of the fits: an eclipse,
    since a grazing chord is shorter than a central one and makes the duration even less possible. rho_LC below
    RHO_EB_RATIO x rho_cat in CONF of the fits: no orbit with e < ECC_MAX around the host gives the duration, so the
    eclipsed star is a larger, lower-density star (blended or evolved). A flat-bottomed fit whose rho_LC agrees with
    rho_cat: the box duration was wrong, and the class follows from Rp = k R_star."""
    k, b, rho = fit['k'], fit['b'], fit['rho']
    note = f"shape fit: rho_LC {_ci(rho)} g/cc against catalogue {rho_cat:.3g}, b {_ci(b)}, k {_ci(k)}"
    graze = np.mean((b > GRAZE_B) & (b + k > 1))
    if graze >= CONF:
        return 'eclipsing binary', f'grazing V-shaped eclipse too long for the host ({100 * graze:.0f} per cent of fits); {note}'
    low = np.mean(rho < RHO_EB_RATIO * rho_cat)
    if low >= CONF:
        return 'eclipsing binary', (f'eclipsed star less dense than the host by more than an orbit with e < {ECC_MAX} '
                                    f'allows (rho_LC < {RHO_EB_RATIO:.3f} rho_cat in {100 * low:.0f} per cent of fits); {note}')
    lo, hi = np.percentile(rho, [100 * (1 - CONF) / 2, 100 * (1 + CONF) / 2])
    flat = np.mean(b + k < 1)
    if flat >= CONF and lo <= rho_cat * RHO_CAT_FAC and hi >= rho_cat / RHO_CAT_FAC and rstar:
        rp = k * rstar / RJ_RSUN
        t = class_from_rp(rp)
        if t:
            return t, f'flat-bottomed fit consistent with the catalogue density, Rp {_ci(rp)} R_J; {note}'
        return None, f'shape fit consistent with the host but Rp {_ci(rp)} R_J spans two classes; {note}'
    return None, f'shape fit inconclusive; {note}'


def shape_norad(fit, teff, tsrc):
    """(class or None, note, rstar) for a host without a catalogue radius. The fitted density is compared with that of
    a dwarf of the host's Teff; when they agree within MS_RHO_FAC (both within the bootstrap interval), the host is
    taken as a dwarf with the radius of the dwarf sequence at rho_LC, and Rp = k R_star gives the class."""
    k, rho = fit['k'], fit['rho']
    note = f"shape fit: rho_LC {_ci(rho)} g/cc, b {_ci(fit['b'])}, k {_ci(k)}"
    if teff is None:
        return None, f'no Teff or Gaia colour for a main-sequence radius; {note}', None
    rms = ms_rho(teff)
    lo, hi = np.percentile(rho, [100 * (1 - CONF) / 2, 100 * (1 + CONF) / 2])
    note += f'; dwarf of Teff {teff:.0f} K ({tsrc}) has rho {rms:.3g} g/cc'
    if hi < rms / MS_RHO_FAC or lo > rms * MS_RHO_FAC:
        return None, f'rho_LC inconsistent with a dwarf of this Teff; {note}', None
    rs = ms_radius(rho)
    rp = k * rs / RJ_RSUN
    t = class_from_rp(rp)
    rmed = float(np.nanmedian(rs))
    if t:
        return t, f'Rp {_ci(rp[np.isfinite(rp)])} R_J from main-sequence R_star {rmed:.2f} R_sun at rho_LC; {note}', rmed
    return None, f'main-sequence R_star {rmed:.2f} R_sun gives Rp {_ci(rp[np.isfinite(rp)]) if np.isfinite(rp).any() else "off the sequence"} R_J; {note}', None


def load_windows(tic):
    """Transit windows kept by the vetting: (t, f, period, t0) or None."""
    p = survey.FIND / sources.name_of(tic).replace(' ', '') / 'transit_windows.npz'
    if not p.exists():
        return None
    z = np.load(p)
    return z['t'], z['f'], float(z['period']), float(z['t0'])


def stellar_teff(con, tic):
    r = con.execute('SELECT teff FROM tic WHERE tic=?' if sources.mission_of(tic) == 'TESS'
                    else 'SELECT teff FROM star WHERE sid=?', (tic,)).fetchone()
    return r[0] if r and r[0] and r[0] > 0 else None


def shape_evolved(fit, rsamp, seed=2):
    """(class or None, note) for an evolved host with a Gaia parallax radius but no mass: the fitted density must lie
    within what GIANT_MASS allows at that radius (an eclipse of a still larger star below it), and Rp = k R_star from
    the bootstrap k times the radius samples gives the class."""
    rng = np.random.default_rng(seed)
    r = rng.choice(np.asarray(rsamp), len(fit['k']))
    note = f"shape fit: rho_LC {_ci(fit['rho'])} g/cc, b {_ci(fit['b'])}, k {_ci(fit['k'])}"
    rho_lo = 1.40980 * GIANT_MASS[0] / r ** 3
    low = np.mean(fit['rho'] < RHO_EB_RATIO * rho_lo)
    if low >= CONF:
        return 'eclipsing binary', (f'eclipsed star less dense than an evolved host of {GIANT_MASS[0]}-{GIANT_MASS[1]} M_sun '
                                    f'allows with e < {ECC_MAX} ({100 * low:.0f} per cent of fits); {note}')
    rp = fit['k'] * r / RJ_RSUN
    t = class_from_rp(rp)
    if t:
        return t, f'Rp {_ci(rp)} R_J = k x Gaia parallax radius; {note}'
    return None, f'Rp {_ci(rp)} R_J (k x Gaia parallax radius) spans two classes; {note}'


def shape_test(con, tic, snr, dur_h, depth_ppm, g, rs=None, m1=None, rsamp=None):
    """(class or None, note) from the transit-shape fit of one undetermined signal: with the host density when it
    has one, with the Gaia parallax radius (rsamp) of an evolved host of unknown mass, and from the dwarf sequence
    when it has no radius."""
    w = load_windows(tic)
    if w is None:
        return None, 'shape fit: no transit windows'
    if not snr or snr < FIT_SNR_MIN:
        return None, f'shape fit: SNR {snr} below {FIT_SNR_MIN}'
    t, f, period, t0 = w
    band = 'TESS' if sources.mission_of(tic) == 'TESS' else 'Kepler'
    rho_cat = stellar_density(rs, m1)
    teff, tsrc, _ = host_teff(con, tic, g)
    if rho_cat:
        logg = 4.438 + np.log10(m1 / rs ** 2)
    elif rsamp is not None:
        logg = 4.438 + np.log10(np.sqrt(np.prod(GIANT_MASS)) / np.median(rsamp) ** 2)
    elif teff:
        tt, _, r, m, _ = ms_dwarfs()
        logg = 4.438 + np.log10(np.interp(teff, tt, m) / np.interp(teff, tt, r) ** 2)
    else:
        logg = None
    u = ld_coeffs(teff, logg, band)
    fit = transit_fit(t, f, period, t0, dur_h, depth_ppm or 0, u=u)
    if not fit['ok']:
        return None, f"shape fit: {fit['note']} (too little)"
    lnote = f' (u {u[0]:.2f}, {u[1]:.2f})'
    if rho_cat:
        c, n = shape_long(fit, rho_cat, rs)
    elif rsamp is not None:
        c, n = shape_evolved(fit, rsamp)
    else:
        c, n, _ = shape_norad(fit, teff, tsrc)
    return c, n + lnote


def stellar_mass(con, tic):
    r = con.execute('SELECT mass FROM tic WHERE tic=?' if sources.mission_of(tic) == 'TESS'
                    else 'SELECT mass FROM star WHERE sid=?', (tic,)).fetchone()
    return r[0] if r and r[0] and r[0] > 0 else None


# ---- host radius from the Gaia DR3 parallax (Stefan-Boltzmann law), for hosts without a catalogue radius
PLX_SNR_MIN = 5.0                   # parallax / error: below this the inverse parallax is a biased distance (Bailer-Jones 2015)
PLX_ZP_MAS = -0.017                 # mean Gaia (E)DR3 parallax zero-point, subtracted (Lindegren et al. 2021)
MBOL_SUN = 4.74                     # IAU 2015 Resolution B2 (Mamajek et al. 2015)
TEFF_SIG_CAT, TEFF_SIG_COL = 150.0, 250.0   # Teff error (K): catalogue or GSP-Phot; dwarf-sequence colour (dwarf/giant offset)
BC_SIG, BC_SIG_EDGE = 0.10, 0.30    # BC_G error inside the 3300-8000 K validity range of the fit, and when clipped to it
AG_SIG = 0.10                       # extinction error (mag) when GSP-Phot gives A_G
N_MC = 2000
# BC_G(Teff) = sum a_i (Teff - 5772 K)^i (Andrae et al. 2018, eq. 7 and Table 8; also used by Gaia DR3 FLAME)
BCG_HOT = (6.000e-02, 6.731e-05, -6.647e-08, 2.859e-11, -7.197e-15)     # 4000-8000 K
BCG_COOL = (1.749e+00, 1.977e-03, 3.737e-07, -8.966e-11, -4.183e-14)    # 3300-4000 K
EVOLVED_R_FAC = 1.5                 # radius above 1.5 x the dwarf radius at that Teff: past the main sequence (subgiant/giant)
GIANT_MASS = (0.8, 3.0)             # mass range of red giants in the APOKASC-2 asteroseismic sample (Pinsonneault et al. 2018)
SPEC_LOGG_GIANT = 3.5               # spectroscopic log g below this marks an evolved star (Huber et al. 2016, EPIC classes)
SPEC_LOGG_SIG_MAX = 0.3             # spectroscopic log g used only if its error is below this (dex)


def bc_g(teff):
    """Gaia G-band bolometric correction at Teff, clipped to the 3300-8000 K range of the Andrae et al. (2018) fit."""
    t = np.clip(np.asarray(teff, float), 3300.0, 8000.0)
    x = t - 5772.0
    hot = sum(a * x ** i for i, a in enumerate(BCG_HOT))
    cool = sum(a * x ** i for i, a in enumerate(BCG_COOL))
    return np.where(t >= 4000.0, hot, cool)


def host_teff(con, tic, g):
    """(Teff, source, error) of the host: catalogue, Gaia DR3 GSP-Phot, or Gaia Bp-Rp on the dwarf sequence."""
    teff = stellar_teff(con, tic)
    if teff:
        return teff, 'catalogue', TEFF_SIG_CAT
    if g and g.get('teff'):
        return g['teff'], 'Gaia DR3 GSP-Phot', TEFF_SIG_CAT
    c = ms_teff(g.get('bp_rp')) if g else None
    if c:
        return c, 'Gaia Bp-Rp on the dwarf sequence', TEFF_SIG_COL
    return None, '', None


def gaia_radius(g, teff, tsig, n=N_MC, seed=1):
    """(radius samples in R_sun or None, note) from L = 4 pi R^2 sigma Teff^4 with the luminosity from the Gaia G
    magnitude, the zero-point-corrected parallax, the GSP-Phot extinction (none when missing, so R is then a lower
    limit) and BC_G(Teff); Monte Carlo over parallax, Teff, A_G and BC_G errors."""
    if not g or teff is None:
        return None, 'no Gaia parallax radius: no Teff or Gaia colour'
    plx, e, gm = g.get('plx'), g.get('plx_e'), g.get('gmag')
    if plx is None or not e or gm is None:
        return None, 'no Gaia parallax radius: no parallax or G magnitude'
    p = plx - PLX_ZP_MAS
    if p / e < PLX_SNR_MIN:
        return None, f'no Gaia parallax radius: parallax/error {p / e:.1f} below {PLX_SNR_MIN:g}'
    rng = np.random.default_rng(seed)
    ps = rng.normal(p, e, n)
    ts = rng.normal(teff, tsig, n)
    ok = (ps > 0) & (ts > 2000)
    ps, ts = ps[ok], ts[ok]
    ag = g.get('ag')
    ags = np.clip(rng.normal(ag, AG_SIG, len(ps)), 0, None) if ag is not None else np.zeros(len(ps))
    inside = 3300.0 <= teff <= 8000.0
    bc = bc_g(ts) + rng.normal(0, BC_SIG if inside else BC_SIG_EDGE, len(ps))
    mbol = gm + 5 * np.log10(ps) - 10 - ags + bc
    r = np.sqrt(10 ** (-0.4 * (mbol - MBOL_SUN))) * (5772.0 / ts) ** 2
    note = (f'Gaia DR3 parallax radius {_ci(r)} R_sun (plx {p:.3f}+-{e:.3f} mas, G {gm:.2f}, '
            f'A_G {ag if ag is not None else 0:.2f}{"" if ag is not None else " assumed: R is a lower limit"}, '
            f'Teff {teff:.0f}+-{tsig:.0f} K{"" if inside else ", BC_G clipped at the fit range"})')
    return r, note


def host_state(rs, teff):
    """'dwarf', 'evolved' or None: CONF of the radius samples below or above EVOLVED_R_FAC x the dwarf radius at Teff."""
    if rs is None or teff is None:
        return None
    t, _, r, _, _ = ms_dwarfs()
    rms = float(np.interp(teff, t, r))
    up = np.mean(np.asarray(rs) > EVOLVED_R_FAC * rms)
    return 'evolved' if up >= CONF else 'dwarf' if up <= 1 - CONF else None


def orbit_inside(rs, period, m1=None):
    """Note when CONF of the radius samples exceed the semi-major axis at the period for the largest mass allowed
    (m1, else the top of GIANT_MASS): a companion would orbit inside the star (Kepler's third law), so the signal
    cannot be an eclipse of, or a transit across, this host. '' otherwise."""
    if rs is None or not period:
        return ''
    m = m1 or GIANT_MASS[1]
    a_rsun = (m * (period / 365.25) ** 2) ** (1 / 3) * 215.03     # 1 au = 215.03 R_sun
    f = np.mean(np.asarray(rs) > a_rsun)
    if f < CONF:
        return ''
    return (f'orbit inside the host: a = {a_rsun:.0f} R_sun at P {period:.2f} d for {m:.1f} M_sun, '
            f'R_star {_ci(rs)} R_sun')


def spec_logg(pos):
    """{tic: (log g, error, survey)} from APOGEE DR17 (III/286) or GALAH DR3 (J/MNRAS/506/150), nearest entry within
    2 arcsec, read from the CfA VizieR mirror."""
    import requests
    from concurrent.futures import ThreadPoolExecutor
    url = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
    cats = (('III/286/catalog', 'logg,e_logg', 'APOGEE DR17'), ('J/MNRAS/506/150/stars', 'logg', 'GALAH DR3'))

    def one(item):
        t, (ra, dec) = item
        for src, cols, name in cats:
            try:
                txt = requests.get(url, params={'-source': src, '-c': f'{ra} {dec}', '-c.rs': 2, '-out': cols,
                                                '-sort': '_r', '-out.max': 1}, timeout=60).text
            except Exception:
                continue
            for line in txt.splitlines():
                p = line.split('\t')
                if line.startswith('#') or not p[0].strip():
                    continue
                try:
                    lg = float(p[0])
                    el = float(p[1]) if len(p) > 1 and p[1].strip() else 0.2   # GALAH DR3 typical log g error
                except ValueError:
                    continue
                return t, (lg, el, name)
        return t, None

    with ThreadPoolExecutor(8) as ex:
        return {t: v for t, v in ex.map(one, pos.items()) if v}


def spec_mass(logg, rstar):
    """Mass (M_sun) from a spectroscopic log g and a radius: M = g R^2 / G, log g_sun = 4.438."""
    if logg is None or not rstar:
        return None
    return float(10 ** (logg - 4.438) * rstar ** 2)


# ---- eccentric orbit from a secondary eclipse (Winn 2010, eqs. 33-34; Kipping 2010)
def ecc_from_secondary(phase, dur_ratio):
    """(e cos w, e sin w, e) from the secondary phase (primary at 0) and T_sec / T_pri, to first order in e."""
    ecw = np.pi / 2 * (phase - 0.5)
    esw = (dur_ratio - 1) / (dur_ratio + 1)
    return ecw, esw, float(min(np.hypot(ecw, esw), 0.95))


def rho_ecc_factor(esw, e):
    """rho_circ / rho_true = ((1 + e sin w) / sqrt(1 - e^2))^3: a circular fit to an eccentric transit misreads the
    stellar density by this factor (Kipping 2010; Dawson & Johnson 2012)."""
    return ((1 + esw) / np.sqrt(1 - e ** 2)) ** 3


# ---- secondary eclipse at any phase and phase curve, from the full light curve (occulta classify --lc)
SEC_FAP = 0.01                      # false-alarm probability of the strongest dip over all trial phases and widths
SEC_WIDTHS = (0.5, 0.75, 1.0, 1.5, 2.0)     # trial secondary durations in units of the primary T14 (e sin w up to ~0.33)
SEC_EXCL = 1.0                      # phases within this many T14 of the primary are not searched
SEC_MIN_ORBITS = 3                  # a secondary must be covered in at least three orbits, each with its own baseline
SEC_TOP = 50                        # the 50 strongest stacked trials are retested orbit by orbit
SEC_SIGMA_MIN = 10.0                # floor on the detection threshold: the same scan on the inverted light curves of the
                                    # 216 searchable ready candidates of 12 Nov 2026 never passed 9.0 sigma (99th
                                    # percentile 8.0), so a dip above 10 sigma is not noise or systematics at the 1/216
                                    # level (inverted-light-curve calibration as in Kepler DR25, Thompson et al. 2018)
ELL_P_MAX = 1.0                     # ellipsoidal test only below this period (d): the 2 d running-median detrending
                                    # keeps a P/2 harmonic only when P/2 is at most a quarter of its window
ELL_ALPHA = 1.37                    # ellipsoidal coefficient 0.15 (15 + u)(1 + tau) / (3 - u), u 0.6, tau 0.4
                                    # (Morris & Naftilan 1993; Shporer 2017)
ELL_SIGMA = 5.0                     # ellipsoidal amplitude must exceed this many (red-noise-inflated) sigma


def _planck_ratio(t1, t2, lam_m=8.0e-7):
    """B_lambda(t1) / B_lambda(t2) at the centre of the TESS and Kepler bands (~0.8 and ~0.65 micron)."""
    c2 = 1.4388e-2 / lam_m
    return np.expm1(c2 / t2) / np.expm1(c2 / t1)


def occultation_max(k, ar, teff):
    """Deepest occultation (secondary) a planet can make: full reflection (geometric albedo 1) plus thermal emission
    at the hottest dayside temperature, no redistribution and zero Bond albedo, T_d = Teff sqrt(R/a) (2/3)^(1/4)
    (Cowan & Agol 2011). Anything deeper is the eclipse of a luminous companion."""
    td = teff * np.sqrt(1 / ar) * (2 / 3) ** 0.25
    return k ** 2 * (1 / ar ** 2 + _planck_ratio(td, teff))


def _beta(res, n_in):
    """Red-noise factor (Pont et al. 2006): scatter of means of n_in points over the white-noise expectation."""
    n_in = max(int(n_in), 1)
    nb = len(res) // n_in
    if nb < 5 or n_in < 2:
        return 1.0
    m = res[:nb * n_in].reshape(nb, n_in).mean(axis=1)
    return max(1.0, float(np.std(m) / (np.std(res) / np.sqrt(n_in))))


def secondary_scan(T, Y, P, t0, t14_h, depth_ppm, teff=5750.0):
    """Strongest box dip away from the primary over trial phases (step T14/4) and widths SEC_WIDTHS x T14. The local
    baseline flanks each box; errors carry a red-noise factor. Detection: false-alarm probability below SEC_FAP after
    all trials (Bonferroni), both halves of the data above SEC_HALF_SIGMA. A detected dip deeper than a planet's
    largest occultation is a stellar eclipse."""
    from scipy.stats import norm
    T, Y = np.asarray(T, float), np.asarray(Y, float)
    t14 = t14_h / 24
    ph = ((T - t0) / P) % 1
    prim = (np.minimum(ph, 1 - ph) * P) < SEC_EXCL * t14
    base = Y[~prim]
    sig0 = 1.4826 * np.median(np.abs(base - np.median(base)))
    beta = _beta(base - np.median(base), max(2, int(np.sum(np.abs(((T - t0) / P + 0.5) % 1 - 0.5) * P < t14 / 2)
                                                    / max(1, len(np.unique(np.round((T - t0) / P)))))))
    half = T < np.median(T)
    orbit = np.floor((T - t0) / P).astype(int)
    lo, hi = (SEC_EXCL * t14) / P, 1 - (SEC_EXCL * t14) / P

    def search(Ys):
        """Best (sigma, phase, width, depth, error, half sigmas, orbits, chi2) of the dips of Ys, and the trial count."""
        trials, n_trials = [], 0
        for wf in SEC_WIDTHS:
            w = wf * t14
            for c in np.arange(lo + w / (2 * P), hi - w / (2 * P), t14 / (4 * P)):
                n_trials += 1
                d = np.abs(ph - c) * P
                inn = (d < 0.35 * w) & ~prim
                out = (d > 0.75 * w) & (d < 0.75 * w + 0.5) & ~prim
                if inn.sum() < 4 or out.sum() < 10:
                    continue
                dep = np.median(Ys[out]) - np.mean(Ys[inn])
                trials.append((dep / (beta * sig0 / np.sqrt(inn.sum())), c, wf, inn, out))

        def per_orbit(inn, out):
            """Depth of each orbit with its own in-box and flanking points: a real eclipse recurs in every orbit, a
            single event (systematic, flare, a neighbour's lone eclipse) does not."""
            ds, es = [], []
            for n in np.unique(orbit[inn]):
                i2, o2 = inn & (orbit == n), out & (orbit == n)
                if i2.sum() >= 3 and o2.sum() >= 6:
                    ds.append(np.median(Ys[o2]) - np.mean(Ys[i2])); es.append(beta * sig0 / np.sqrt(i2.sum()))
            return np.array(ds), np.array(es)

        best = None
        for _, c, wf, inn, out in sorted(trials, key=lambda x: -x[0])[:SEC_TOP]:
            ds, es = per_orbit(inn, out)
            if len(ds) < SEC_MIN_ORBITS:
                continue
            wt = 1 / es ** 2
            dep = float(np.sum(wt * ds) / np.sum(wt))
            chi = float(np.sum(wt * (ds - dep) ** 2) / (len(ds) - 1))
            err = float(np.sqrt(max(chi, 1.0) / np.sum(wt)))     # scatter between orbits inflates the error
            s = dep / err
            if best is None or s > best[0]:
                hs = []
                for hm in (half, ~half):
                    dd, ee = per_orbit(inn & hm, out & hm)
                    if len(dd):
                        w2 = 1 / ee ** 2
                        hs.append(float(np.sum(w2 * dd) / np.sqrt(np.sum(w2))))
                best = (s, c, wf, dep, err, hs, len(ds), chi)
        return best, n_trials, bool(trials)

    best, n_trials, any_ = search(Y)
    if best is None:
        return {'ok': False, 'note': 'secondary scan: no phase covered' if not any_ else
                f'secondary scan: fewer than {SEC_MIN_ORBITS} orbits cover any trial phase'}
    # the same search on the inverted light curve: its strongest "dip" (a brightening) measures how often noise and
    # systematics alone pass, as with the inverted light curves of Kepler DR25 (Thompson et al. 2018)
    inv = search(-Y)[0]
    s_inv = float(inv[0]) if inv else 0.0
    s, c, wf, dep, err, hs, n_orb, chi = best
    thr = max(float(norm.isf(SEC_FAP / max(n_trials, 1))), SEC_SIGMA_MIN)
    k = np.sqrt(max(depth_ppm, 1.0) / 1e6)
    ar = P / (np.pi * t14) * (1 + k)
    dmax = float(occultation_max(k, ar, teff or 5750.0))
    det = bool(s > thr and s > s_inv and len(hs) == 2 and min(hs) > 2.0)
    ecw, esw, e = ecc_from_secondary(c, wf)
    return {'ok': True, 'sigma': float(s), 'sigma_inv': s_inv, 'threshold': thr, 'phase': float(c),
            'width_ratio': float(wf), 'depth_ppm': float(dep * 1e6), 'err_ppm': float(err * 1e6),
            'halves': [float(h) for h in hs], 'detected': det, 'stellar': bool(det and dep > dmax),
            'occ_max_ppm': dmax * 1e6, 'beta': beta, 'n_trials': n_trials, 'n_orbits': n_orb, 'chi2_orbits': chi,
            'ecosw': float(ecw), 'esinw': float(esw), 'ecc': e}


def phase_curve(T, Y, P, t0, t14_h, depth_ppm, m1=1.0):
    """Ellipsoidal, reflection and beaming harmonics out of eclipse, y = c + a1 cos(2 pi phi) + b1 sin(2 pi phi) +
    a2 cos(4 pi phi), for P <= ELL_P_MAX. The ellipsoidal amplitude -a2 gives q = A / (ELL_ALPHA (R/a)^3) at sin i ~ 1
    (Shporer 2017, eq. 6); a companion above MSTAR_MIN is a star."""
    if P > ELL_P_MAX:
        return {'ok': False, 'note': f'phase curve: P {P:.2f} d above {ELL_P_MAX:g} d (detrending removes it)'}
    T, Y = np.asarray(T, float), np.asarray(Y, float)
    t14 = t14_h / 24
    ph = ((T - t0) / P) % 1
    m = (np.minimum(ph, 1 - ph) * P > t14) & (np.abs(ph - 0.5) * P > t14)
    if m.sum() < 50:
        return {'ok': False, 'note': 'phase curve: too few points'}
    x = 2 * np.pi * ph[m]
    A = np.c_[np.ones(m.sum()), np.cos(x), np.sin(x), np.cos(2 * x)]
    coef, *_ = np.linalg.lstsq(A, Y[m], rcond=None)
    res = Y[m] - A @ coef
    cov = np.linalg.inv(A.T @ A) * np.var(res) * _beta(res, max(2, int(0.1 * P / np.median(np.diff(T)))) ) ** 2
    a_ell, e_ell = float(-coef[3]), float(np.sqrt(cov[3, 3]))
    k = np.sqrt(max(depth_ppm, 1.0) / 1e6)
    ar = P / (np.pi * t14) * (1 + k)
    q = a_ell / (ELL_ALPHA * ar ** -3)
    m2 = q * (m1 or 1.0)
    det = a_ell > ELL_SIGMA * e_ell
    return {'ok': True, 'ell_ppm': a_ell * 1e6, 'ell_err_ppm': e_ell * 1e6, 'refl_ppm': float(-coef[1] * 1e6),
            'beam_ppm': float(coef[2] * 1e6), 'q': float(q), 'm2': float(m2), 'detected': bool(det),
            'stellar': bool(det and m2 > MSTAR_MIN)}


def lc_checks(tics=None, status='ready', keep=False):
    """Download the full light curves of the candidates, run secondary_scan and phase_curve on the detrended
    series of the vetting, and store both in summary.json (keys secondary_any, phase_curve). Nothing else changes."""
    import shutil
    from . import planets
    con = survey.db()
    rows = con.execute('SELECT tic, period, t0, dur_h, depth_ppm FROM candidate WHERE status=?', (status,)).fetchall()
    if tics:
        rows = [r for r in rows if r[0] in set(tics)]
    for tic, P, t0, dur, dep in rows:
        p = survey.FIND / sources.name_of(tic).replace(' ', '') / 'summary.json'
        if not p.exists():
            print(tic, 'no summary.json'); continue
        s = json.loads(p.read_text())
        P, t0 = s.get('period_d') or P, s.get('t0_bjd') or t0
        t14, dep = s.get('t14_h') or dur, s.get('depth_ppm') or dep
        work = survey.CACHE / f'cls_{tic}'
        work.mkdir(parents=True, exist_ok=True)
        try:
            cached = work / 'series.npz'
            if cached.exists():                  # detrended series of an earlier --keep run
                z = np.load(cached)
                T, Y = z['T'], z['Y']
            else:
                data = {}
                for author, unit, path in planets.fetch_all(tic, work):
                    ser = planets.series(author, path)
                    if ser is not None:
                        data.setdefault(unit, {})[author] = ser
                main = [(d.get('QLP') or d.get('TESS-SPOC') or d.get('SPOC') or d.get('Kepler') or d.get('K2SFF'))
                        for d in data.values()]
                main = [m for m in main if m is not None]
                if not main:
                    raise RuntimeError('no light curves')
                T = np.concatenate([m[0] for m in main]); Y = np.concatenate([m[1] for m in main])
                o = np.argsort(T); T, Y = T[o], Y[o]
                if keep:
                    np.savez_compressed(cached, T=T, Y=Y)
            t0b = t0 - planets.BTJD if t0 > 2e6 else t0
            teff = stellar_teff(con, tic) or 5750.0
            s['secondary_any'] = secondary_scan(T, Y, P, t0b, t14, dep, teff)
            s['phase_curve'] = phase_curve(T, Y, P, t0b, t14, dep, stellar_mass(con, tic) or 1.0)
            print(tic, 'secondary', round(s['secondary_any'].get('sigma', np.nan), 1), s['secondary_any'].get('detected'),
                  'ell', s['phase_curve'].get('detected'))
        except Exception as e:
            s['secondary_any'] = {'ok': False, 'note': f'secondary scan failed: {e}'}
            print(tic, 'failed', e)
        p.write_text(json.dumps(s, indent=1, default=float))
        if not keep:
            shutil.rmtree(work, ignore_errors=True)


def lc_tests(tic):
    """(secondary_any, phase_curve) dicts from the vetting summary, {} when lc_checks has not run."""
    p = survey.FIND / sources.name_of(tic).replace(' ', '') / 'summary.json'
    if not p.exists():
        return {}, {}
    s = json.loads(p.read_text())
    return s.get('secondary_any') or {}, s.get('phase_curve') or {}


def lc_test(sa, pc):
    """(is_binary, note) from the any-phase secondary and the ellipsoidal variation."""
    notes = []
    if sa.get('ok'):
        # results stored by an older run are held to the present floor
        sa = {**sa, 'detected': bool(sa.get('detected') and sa['sigma'] > SEC_SIGMA_MIN),
              'threshold': max(sa['threshold'], SEC_SIGMA_MIN)}
        sa['stellar'] = bool(sa['detected'] and sa['depth_ppm'] > sa['occ_max_ppm'])
        if sa.get('stellar'):
            return True, (f"secondary eclipse at phase {sa['phase']:.3f} ({sa['sigma']:.1f} sigma, threshold "
                          f"{sa['threshold']:.1f}; {sa['depth_ppm']:.0f} ppm, above the {sa['occ_max_ppm']:.0f} ppm a planet "
                          f"can occult; e cos w {sa['ecosw']:+.2f}, e sin w {sa['esinw']:+.2f})")
        if sa.get('detected'):
            notes.append(f"occultation-depth dip at phase {sa['phase']:.3f} ({sa['depth_ppm']:.0f} ppm, {sa['sigma']:.1f} sigma)")
        else:
            notes.append(f"no secondary at any phase (best {sa['sigma']:.1f} sigma, threshold {sa['threshold']:.1f})")
    if pc.get('ok'):
        if pc.get('stellar'):
            return True, (f"ellipsoidal variation {pc['ell_ppm']:.0f}+-{pc['ell_err_ppm']:.0f} ppm, companion "
                          f"{pc['m2']:.2f} M_sun")
        notes.append(f"ellipsoidal {pc['ell_ppm']:.0f}+-{pc['ell_err_ppm']:.0f} ppm (M2 < {MSTAR_MIN} M_sun or not detected)")
    return False, '; '.join(notes)


def run(status='ready', tics=None, gaia=True, external=True):
    con = survey.db()
    if 'objtype' not in [r[1] for r in con.execute('PRAGMA table_info(candidate)')]:
        con.execute('ALTER TABLE candidate ADD COLUMN objtype TEXT')
    rows = con.execute('SELECT tic, depth_ppm, notes, period, t0, dur_h, snr FROM candidate WHERE status=?',
                       (status,)).fetchall()
    if tics:
        rows = [r for r in rows if r[0] in set(tics)]
    fill_tic(con, [r[0] for r in rows if stellar_radius(con, r[0]) is None])
    pos = positions(con, [r[0] for r in rows])
    gd = gaia_dr3(con, [r[0] for r in rows]) if gaia else {}
    if gaia and pos and not gd:
        # without Gaia the orbit, RV and eclipsing-binary tests vanish and binaries fall back to planets: write nothing
        raise RuntimeError('Gaia DR3 cross-match failed or matched no host; no class written (retry, or --no-gaia)')
    ap = apogee(pos) if external else {}
    lm = lamost(pos) if external else {}
    zl = ztf_all(pos, [(r[0], r[1]) for r in rows]) if external else {}
    sl = spec_logg(pos) if external else {}
    print(f'Gaia {len(gd)}, APOGEE {len(ap)}, LAMOST {len(lm)}, ZTF {len(zl)}, spectroscopic log g {len(sl)} hosts matched')
    counts = {}
    for tic, dep, notes, period, t0, dur, snr in rows:
        mission = sources.mission_of(tic)
        zeb, znote, zd = ztf_test(zl.get(tic), mission, period, t0, dur, dep or 0)
        rs, m1, src = stellar_radius(con, tic), stellar_mass(con, tic), 'TIC v8.2'
        g = gd.get(tic) or {}
        hnotes, rsamp = [], None
        if rs is None and g.get('rad'):          # no TIC radius: Gaia DR3 FLAME or GSP-Phot
            rs, src = g['rad'], 'Gaia DR3'
            m1 = m1 or g.get('mass')
        teff, _, tsig = host_teff(con, tic, g)
        if rs is None:                           # nor FLAME/GSP-Phot: luminosity radius from the parallax
            rsamp, rnote = gaia_radius(g, teff, tsig)
            hnotes.append(rnote)
            if rsamp is not None:
                rs, src = float(np.median(rsamp)), 'Gaia DR3 parallax'
                state = host_state(rsamp, teff)
                hnotes.append(f'host {state or "dwarf or evolved (undecided)"} at Teff {teff:.0f} K')
                if state == 'dwarf' and not m1:
                    tt, _, _, mm, _ = ms_dwarfs()
                    m1 = float(np.interp(teff, tt, mm))
        sp = sl.get(tic)
        if sp and sp[1] <= SPEC_LOGG_SIG_MAX and rs:
            ms = spec_mass(sp[0], rs)
            hnotes.append(f'{sp[2]} log g {sp[0]:.2f}+-{sp[1]:.2f}, mass g R^2/G {ms:.2f} M_sun')
            if not m1 or sp[0] < SPEC_LOGG_GIANT:   # TIC masses come from dwarf relations: wrong for evolved hosts
                m1 = ms
        if rsamp is not None and m1:
            rsamp = None                         # mass known: the density tests below apply as for catalogue hosts
        dk = zd if dep and zd > ZTF_DILUTION * dep else dep
        if dep and zd > ZTF_DILUTION * dep:      # ZTF's 1 arcsec seeing resolves contaminants that TESS pixels blend
            t, why = label(zd, rs, *vet_tests(tic), src=src)
            why += f' (ZTF undiluted depth, TESS/Kepler depth {dep:.0f} ppm)'
        else:
            t, why = label(dep, rs, *vet_tests(tic), src=src)
        inside = ''
        if src == 'Gaia DR3 parallax' and dk:
            r_all = gaia_radius(g, teff, tsig)[0]
            inside = orbit_inside(r_all, period, m1)
            if inside:                           # no companion can orbit this star at the period
                t, why = 'unclassified', f'{inside}: the signal is not from this host; was {t}: {why}'
            elif t != 'eclipsing binary':        # the class must hold over the radius errors
                tc = class_from_rp(np.sqrt(dk / 1e6) * r_all / RJ_RSUN)
                if tc is None:
                    t, why = 'unclassified', f'radius error spans two classes: {why}'
                else:
                    t = tc
        if hnotes:
            why += '; ' + '; '.join(hnotes)
        ceb, cnote = contact_test(dur, period)
        if ceb and t != 'eclipsing binary':
            t, why = 'eclipsing binary', f'{cnote}; was {t}: {why}'
        # a planet or brown-dwarf size rests on the catalogued host; a duration the host density cannot produce means
        # that host is not the eclipsed star (or its parameters are wrong), so the size is undetermined
        long_, dnote = duration_test(dur, period, stellar_density(rs, m1), dk or 0)
        sa, pc = lc_tests(tic)
        if long_ and sa.get('detected') and not sa.get('stellar'):
            # an occultation-depth secondary fixes e cos w and e sin w; a circular reading of an eccentric transit
            # gives the density rho_true x rho_ecc_factor, so the duration is retested against that density
            f_e = rho_ecc_factor(sa['esinw'], sa['ecc'])
            long_, dnote = duration_test(dur, period, stellar_density(rs, m1) * f_e, dk or 0)
            dnote += f" (eccentric: e {sa['ecc']:.2f}, density factor {f_e:.2f} from the secondary)"
        if long_ and t in ('planet', 'brown dwarf or low-mass star'):
            beb, bnote = blend_test(vet_blends(tic))
            if beb:                              # the eclipsed star is a resolved neighbour, not the host
                t, why = 'eclipsing binary', f'{bnote}; too long for the host: {dnote}; was {t}: {why}'
            else:
                t, why = 'unclassified', f'too long for the host: {dnote}; was {t}: {why}'
        else:
            why += '; ' + ('too long for the host: ' if long_ else '') + dnote
        if t == 'unclassified' and why != 'no transit depth' and not inside:   # radius missing or too long: shape
            use = long_ or (rs and m1 and src != 'TIC v8.2')
            st, snote = shape_test(con, tic, snr, dur, dk or 0, g, rs if use else None, m1 if use else None,
                                   rsamp=rsamp)
            t, why = (st, f'{snote}; was unclassified: {why}') if st else (t, f'{why}; {snote}')
        tests = [gaia_test(gd.get(tic), period, m1), gaia_rv(gd.get(tic), period, m1),
                 apogee_test(ap.get(tic), period, m1), lamost_test(lm.get(tic), period, m1), (zeb, znote),
                 lc_test(sa, pc)]
        for eb, note in tests:
            if eb:
                t, why = 'eclipsing binary', note + '; ' + why
            elif note:
                why += '; ' + note
        kept = strip_class(notes)
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
    ap.add_argument('--lc', action='store_true', help='only download the light curves and store the any-phase '
                    'secondary and phase-curve tests in summary.json (no class change)')
    ap.add_argument('--keep', action='store_true', help='with --lc: keep the downloaded light curves')
    ap.add_argument('tic', type=int, nargs='*')
    a = ap.parse_args(argv)
    if a.lc:
        lc_checks(a.tic, a.status, keep=a.keep)
        return
    run(a.status, a.tic, gaia=not a.no_gaia, external=not a.no_external)
