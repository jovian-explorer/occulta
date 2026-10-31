"""Single transits and aperiodic dippers: clean dips that the survey saw once in all the data of a star, with the shortest
period that the observed time span still allows and the period that the transit duration implies through the stellar
density; stars with several clean dips that do not match one another in depth and duration (aperiodic dippers); and a
shape test that separates a symmetric transit from the sharp ingress and slow egress of a disintegrating planet's tail.

    occulta singles scan [--snr S]   every star in the database, from the dips that the survey already measured (table
                                     event); nothing is downloaded.
    occulta singles search TIC ...   download every light curve of the stars, measure the shape of each clean dip, and
                                     classify it as a transit or an asymmetric dip.
    occulta singles list             the table single, strongest first.

Results go to the table single and to findings/singles/ (singles.csv, report.md). The survey runs the scan after each
pass when "search_singles": true is set in config.json.
"""
import argparse, csv, json, time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from . import lightcurves as lcs
from .planets import trapezoid_ss

G, RSUN, MSUN, REARTH = 6.674e-11, 6.957e8, 1.989e30, 6.371e6
SNR_SINGLE = 9.0        # a lone dip needs more than the 7 sigma of the event search, as no repeat confirms it
SAME_D = 0.5            # dips closer than this (two products of one sector) are one dip
B_HI = 0.8              # the upper end of the duration period assumes this impact parameter
P_MAX_D = 3000.0        # longest period tried for the coverage limit
DCHI2_TAIL = 25.0       # the tail profile (as many parameters as the trapezoid) must fit better by this much
COLS = ('tic', 'sector', 't0', 'dur_h', 'depth_ppm', 'snr', 'kind', 'p_min_d', 'p_dur_d', 'p_dur_lo', 'p_dur_hi',
        'rp_earth', 'asym_sigma', 'notes')


# ---------------------------------------------------------------- period from one transit
def period_from_duration(t14_h, k, rstar, mstar, b=0.0):
    """Period (d) of a circular orbit whose transit lasts t14_h hours, for radius ratio k and impact parameter b:
    T14 = P R* sqrt((1+k)^2 - b^2) / (pi a) with Kepler's third law gives P = G M* pi T14^3 / (4 R*^3 L^3),
    L = sqrt((1+k)^2 - b^2) (Seager & Mallen-Ornelas 2003, small-angle form)."""
    L = np.sqrt(np.maximum((1 + k) ** 2 - b ** 2, 1e-6))
    T = t14_h * 3600.0
    return G * mstar * MSUN * np.pi * T ** 3 / (4 * (rstar * RSUN) ** 3 * L ** 3) / 86400.0


def min_period(t0, dur_h, windows, p_lo=None):
    """Shortest period for which no other transit falls inside an observed window [(start, end)] (BTJD), each window
    shrunk by half the duration so that a partly covered transit is not counted. Gaps inside a window are treated as
    observed, so the limit can be too long when a repeat fell into a data gap."""
    half = dur_h / 48.0
    W = np.array([(a + half, b - half) for a, b in windows if b - a > 2 * half], float).reshape(-1, 2)
    p_lo = p_lo or max(2 * dur_h / 24.0, 0.3)
    P = np.geomspace(p_lo, P_MAX_D, 40000)
    lo = np.ceil((W[None, :, 0] - t0) / P[:, None])
    hi = np.floor((W[None, :, 1] - t0) / P[:, None])
    n = np.maximum(hi - lo + 1, 0) - ((lo <= 0) & (hi >= 0))        # repeats other than the dip itself
    ok = np.all(n <= 0, axis=1)
    return float(P[np.argmax(ok)]) if ok.any() else P_MAX_D


# ---------------------------------------------------------------- shape
def tail(t, ta, w, d, tau):
    """Dip of a dusty tail: linear ingress over w to the depth d, then exponential recovery with time constant tau
    (the profile shape of the disintegrating planets, Rappaport et al. 2012)."""
    x = t - ta
    ramp = np.clip(x / max(w, 1e-6), 0, 1)
    rec = np.exp(-np.clip(x - w, 0, None) / max(tau, 1e-6))
    return -d * np.where(x < w, ramp, rec)


def shape(t, y, t0, dur_h, sig):
    """Symmetric trapezoid against the tail profile on the window around one dip; returns (asym_sigma, dchi2), dchi2 > 0
    when the tail profile fits better."""
    from .rogues import mirror_asymmetry
    d = dur_h / 24.0
    m = np.abs(t - t0) < 1.5 * d + 0.25
    t, y = t[m], y[m]
    if t.size < 20:
        return np.nan, np.nan
    cad = float(np.median(np.diff(t)))
    dep = max(-np.min(np.convolve(y, np.ones(5) / 5, 'same')), 3 * sig)
    trap = lambda p: (y - p[4] - trapezoid_ss(t, p[0], p[1], p[2], p[3] * p[2] / 2, cad)) / sig
    ft = least_squares(trap, [t0, dep, d, 0.3, 0], bounds=([t0 - d, 0, d / 4, 0.02, -1], [t0 + d, 1, 3 * d, 1, 1]))
    tl = lambda p: (y - p[4] - tail(t, p[0], p[1], p[2], p[3])) / sig
    best = None
    for frac in (0.1, 0.3):
        f = least_squares(tl, [t0 - d / 2, frac * d, dep, d / 2, 0],
                          bounds=([t0 - 2 * d, cad / 2, 0, cad / 2, -1], [t0 + d, 2 * d, 1, 3 * d, 1]))
        best = f if best is None or f.cost < best.cost else best
    dchi2 = 2 * (ft.cost - best.cost)
    tc = ft.x[0]
    return mirror_asymmetry(t, -(y - ft.x[4]), tc, 1.5 * d, sig), float(dchi2)


# ---------------------------------------------------------------- survey tables
def _star(con, tic):
    r = con.execute('SELECT rad, mass FROM tic WHERE tic=?', (tic,)).fetchone()
    rad = r[0] if r and r[0] else None
    mass = r[1] if r and r[1] else None
    return rad, mass


def classify(con, snr=SNR_SINGLE, inv_ratio=1.0):
    """Lone and mutually inconsistent clean dips from the table event: [(row dict)]. Stars that the survey already made
    candidates of, and dips that match a known object or eclipsing binary, are left out."""
    from .survey import known_match
    q = ('SELECT e.tic, e.sector, e.t0, e.dur_h, e.depth_ppm, e.snr, ss.inv_snr FROM event e JOIN star_sector ss '
         'USING(tic, sector) WHERE e.snr >= ? AND e.flags = "" AND NOT EXISTS (SELECT 1 FROM candidate c WHERE c.tic=e.tic) '
         'ORDER BY e.tic, e.snr DESC')
    by = {}
    for r in con.execute(q, (min(snr, 7.0),)):
        if r[5] >= inv_ratio * (r[6] or 0) and known_match(con, r[0], r[2]) not in ('match', 'eb'):
            by.setdefault(r[0], []).append(r)
    out = []
    for tic, evs in by.items():
        uniq = []
        for e in evs:                                          # sorted by decreasing S/N
            if all(abs(e[2] - u[2]) >= SAME_D for u in uniq):
                uniq.append(e)
        if uniq[0][5] < snr:
            continue
        windows = [(a, b) for a, b in con.execute('SELECT t_start, t_end FROM star_sector WHERE tic=? AND t_start IS NOT NULL',
                                                  (tic,))]
        rad, mass = _star(con, tic)
        notes = [] if rad and mass else ['no stellar radius or mass in the catalogue: the Sun assumed']
        rad, mass = rad or 1.0, mass or 1.0
        kind = 'single' if len(uniq) == 1 else 'aperiodic'
        for e in (uniq[:1] if kind == 'single' else uniq):
            k = np.sqrt(max(e[4], 0) / 1e6)
            rp = k * rad * RSUN / REARTH
            row = dict(zip(COLS, (tic, e[1], e[2], e[3], e[4], e[5], kind, None, None, None, None, rp, None, '')))
            if kind == 'single':
                row['p_min_d'] = min_period(e[2], e[3], windows)
                row['p_dur_lo'] = float(period_from_duration(e[3], k, rad, mass, 0.0))
                row['p_dur_hi'] = float(period_from_duration(e[3], k, rad, mass, B_HI))
                row['p_dur_d'] = float(period_from_duration(e[3], k, rad, mass, 0.5))
                n = list(notes)
                if row['p_dur_hi'] < row['p_min_d']:
                    n.append('the duration implies a shorter period than the coverage allows: an eccentric orbit, a '
                             'grazing transit, or a wrong stellar radius')
                if rp > 25:
                    n.append(f'radius {rp:.0f} Re: a stellar companion is more likely than a planet')
                row['notes'] = '; '.join(n)
            else:
                row['notes'] = f'{len(uniq)} dips of different depth or duration'
            out.append(row)
    return out


def store(con, rows, findings):
    for r in rows:
        con.execute(f'INSERT OR REPLACE INTO single VALUES ({",".join("?" * len(COLS))})', tuple(r[c] for c in COLS))
    con.commit()
    write(con, findings)


def write(con, findings):
    d = Path(findings) / 'singles'
    d.mkdir(parents=True, exist_ok=True)
    rows = [dict(zip(COLS, r)) for r in con.execute(f'SELECT {",".join(COLS)} FROM single ORDER BY kind, snr DESC')]
    with open(d / 'singles.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, COLS)
        w.writeheader()
        w.writerows(rows)
    fmt = lambda x, f: '' if x is None else format(x, f)
    L = ['# Single transits and aperiodic dips', '',
         f'{sum(r["kind"] == "single" for r in rows)} single transits, '
         f'{len({r["tic"] for r in rows if r["kind"] == "aperiodic"})} stars with aperiodic dips, '
         f'{sum(r["kind"] == "dipper" for r in rows)} asymmetric (tail-shaped) dips. Updated {time.strftime("%Y-%m-%d %H:%M")}.', '',
         '| TIC | t0 (BTJD) | duration (h) | depth (ppm) | S/N | kind | P min (d) | P from duration (d) | Rp (Re) | notes |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        pd = '' if r['p_dur_d'] is None else f"{r['p_dur_d']:.0f} ({r['p_dur_lo']:.0f}-{r['p_dur_hi']:.0f})"
        L.append(f"| {r['tic']} | {r['t0']:.3f} | {r['dur_h']:.1f} | {r['depth_ppm']:.0f} | {r['snr']:.1f} | {r['kind']} | "
                 f"{fmt(r['p_min_d'], '.1f')} | {pd} | {fmt(r['rp_earth'], '.1f')} | {r['notes'] or ''} |")
    (d / 'report.md').write_text('\n'.join(L) + '\n', encoding='utf-8')


def run_survey(con, findings, snr=SNR_SINGLE):
    from .survey import CFG
    rows = classify(con, snr, CFG.get('inv_ratio', 1.0))
    store(con, rows, findings)
    return len(rows)


# ---------------------------------------------------------------- light curves
def search_lc(t, f, events):
    """Shape of each clean dip in one light curve: [(t0, asym_sigma, dchi2, kind)]."""
    o = np.argsort(t)
    t, f = t[o], f[o]
    y, _ = lcs.detrend(t, f)
    ok = np.isfinite(y)
    t, y = t[ok], y[ok]
    sig = lcs.mad_std(y)
    out = []
    for e in events:
        a, dc = shape(t, y, e['t0'], e['dur_h'], sig)
        kind = 'dipper' if dc >= DCHI2_TAIL else 'transit'
        out.append((e['t0'], a, dc, kind))
    return out


def search_star(con, tic, findings):
    """Download every light curve of a star, run the dip search and the shape test; dips classified as tail-shaped are
    stored with kind 'dipper', the others update the asymmetry of a stored single."""
    import shutil
    from . import planets, rogues
    work = Path(findings).parent / 'data' / 'cache' / f'single_{tic}'
    work.mkdir(parents=True, exist_ok=True)
    res = []
    try:
        for author, num, path in planets.fetch_all(tic, work):
            unit = num % 1000 + rogues.UNIT_OFFSET[author]
            fmt = planets.AUTHORS[author]
            r = lcs.process(path, tic, unit, fmt, do_bls=False)
            if 'error' in r:
                continue
            clean = [e for e in r['events'] if e['flags'] == '']
            t, f, *_ = lcs.read_lc(path, fmt)
            for t0, a, dc, kind in search_lc(t, f, clean):
                e = next(x for x in clean if x['t0'] == t0)
                res.append((unit, e, a, dc, kind))
                if kind == 'dipper':
                    rad, _ = _star(con, tic)
                    rp = np.sqrt(max(e['depth_ppm'], 0) / 1e6) * (rad or 1.0) * RSUN / REARTH
                    con.execute(f'INSERT OR REPLACE INTO single VALUES ({",".join("?" * len(COLS))})',
                                (tic, unit, t0, e['dur_h'], e['depth_ppm'], e['snr'], 'dipper', None, None, None, None,
                                 rp, a, f'tail profile better by delta chi2 {dc:.0f}'))
                else:
                    con.execute('UPDATE single SET asym_sigma=? WHERE tic=? AND ABS(t0-?) < ?', (a, tic, t0, SAME_D))
        con.commit()
        write(con, findings)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return res


def main(argv=None):
    from . import paths, survey
    ap = argparse.ArgumentParser(prog='occulta singles', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='sub', required=True)
    s = sp.add_parser('scan', help='single transits and aperiodic dips from the dips already in the database')
    s.add_argument('--snr', type=float, default=SNR_SINGLE, help=f'S/N of a lone dip (default {SNR_SINGLE})')
    s = sp.add_parser('search', help='download the light curves of stars and test the shape of each clean dip')
    s.add_argument('tic', nargs='+', type=int)
    sp.add_parser('list', help='the table single')
    a = ap.parse_args(argv)
    paths.init()
    con = survey.db()
    if a.sub == 'scan':
        n = run_survey(con, survey.FIND, a.snr)
        print(f'{n} single or aperiodic dips; see {Path(survey.FIND) / "singles" / "report.md"}')
    elif a.sub == 'search':
        for tic in a.tic:
            for unit, e, asym, dc, kind in search_star(con, tic, survey.FIND):
                print(tic, unit, f"t0 {e['t0']:.3f}", f"{e['dur_h']:.1f} h", f"{e['depth_ppm']:.0f} ppm",
                      f'asymmetry {asym:.1f} sigma', f'delta chi2 {dc:.0f}', kind)
    else:
        for r in con.execute(f'SELECT {",".join(COLS)} FROM single ORDER BY snr DESC'):
            print(json.dumps(dict(zip(COLS, r)), default=float))
