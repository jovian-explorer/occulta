"""Candidate vetting: re-fetch every sector of one star, resolve the period, fit the
transit shape, run the standard false-positive tests, write a report with plots,
keep only a compact transit-window light curve and purge the downloads."""
import json, shutil, time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from . import lightcurves as lcsearch
from . import sources as sources

BTJD = 2457000.0
G = 6.674e-11
RSUN, MSUN, REARTH = 6.957e8, 1.989e30, 6.371e6


AUTHORS = {'QLP': 'qlp', 'SPOC': 'spoc', 'TESS-SPOC': 'spoc', 'Kepler': 'kepler', 'K2SFF': 'k2sff'}
P_MIN = 0.5            # shortest period proposed from a pair of dips (d)
MAX_PAIR_DIPS = 25     # only the strongest dips propose periods
MAX_PANELS = 24       # individual transits drawn in transits.png
P_HALF_MIN = 0.2       # shortest P/2 trialled for a BLS period (d); the BLS grid starts at 0.5 d, so true periods below that appear only doubled
REFINE_SPAN, REFINE_TRIALS = 30.0, 241  # BLS period refit: +-30 steps of P dur / B, 4 trials per step, 9 epochs each (2,169 folds); LHS 3844 b needed more than 12 steps
DETECT_SNR = 7.0            # per-transit S/N of the candidate depth needed before an undetected epoch counts as absent
STACK_SNR = 3.0             # the stack of all covered epochs of a period must show the transit at this S/N
HALF_ODD_EVEN_SIGMA = 3.0  # P/2 is kept only if its two parities agree in depth within this many sigma
V_MIN_CADENCES = 8.0       # the V-shape cut is used only if the transit spans this many median exposures: a rounded (limb-darkened) or timing-smeared transit of fewer cadences always fits as a V
V_PREFER_SIGMA = 3.0       # V-shape cut: the V (ingress = half the transit) must beat a flat-bottomed transit by this many sigma (delta chi2 = sigma^2) after cadence integration
SEC_HALF_SIGMA = 2.0       # a secondary eclipse rejects the candidate only if each half of the data (early, late) also shows it at this significance
VET_LIMIT_S = 900      # a period search that takes longer marks the candidate vet_error instead of blocking the run


class VetTimeout(Exception):
    pass


def fetch_all(sid, workdir):
    """Download every long-baseline light curve of the star: [(author, unit, path)]."""
    import re, requests, lightkurve as lk
    res = lk.search_lightcurve(sources.name_of(sid))
    out = []
    tab = res.table
    for i in range(len(tab)):
        author = str(tab['author'][i])
        exp = float(tab['exptime'][i])
        if author not in AUTHORS or exp < 100 or (author == 'Kepler' and exp < 1000):
            continue
        num = int(re.findall(r'\d+', str(tab['mission'][i]))[-1])
        unit = 3000 + num if author == 'Kepler' else 2000 + num if author == 'K2SFF' else num
        uri = str(tab['dataURI'][i])
        path = Path(workdir) / uri.split('/')[-1]
        try:
            r = requests.get('https://mast.stsci.edu/api/v0.1/Download/file', params={'uri': uri}, timeout=300)
            if r.status_code != 200 or r.content[:6] != b'SIMPLE':
                continue
            path.write_bytes(r.content)
        except Exception:
            continue
        out.append((author, unit, str(path)))
    return out


def series(author, path):
    """(t [BTJD], detrended flux) for one light curve, using the same detrending as the search."""
    t, f, *_ = lcsearch.read_lc(path, AUTHORS[author])
    o = np.argsort(t); t, f = t[o], f[o]
    if t.size < 100:
        return None
    y, _ = lcsearch.detrend(t, f)
    m = np.isfinite(y)
    return t[m], y[m]


def get_star(con, sid):
    """(ra, dec, mag, teff, rad, mass, dist, gaia, contratio) from the local tables, else the archives."""
    m = sources.mission_of(sid)
    if m == 'TESS':
        r = con.execute('SELECT ra, dec, tmag, teff, rad, mass, dist, gaia, contratio FROM tic WHERE tic=?', (sid,)).fetchone()
        if r:
            return r
        try:
            from astroquery.mast import Catalogs
            t = Catalogs.query_criteria(catalog='Tic', ID=sid)[0]
            f = lambda k: float(t[k]) if np.isfinite(t[k]) else None
            return (f('ra'), f('dec'), f('Tmag'), f('Teff'), f('rad'), f('mass'), f('d'), str(t['GAIA']), f('contratio'))
        except Exception:
            return None
    r = con.execute('SELECT ra, dec, mag, teff, rad, mass, dist FROM star WHERE sid=?', (sid,)).fetchone()
    if r:
        return (*r, None, None)
    try:
        import csv, io, requests
        q = (f"select ra, dec, k2_kepmag, k2_teff, k2_rad, k2_mass, k2_dist from k2targets where epic_number={sid - sources.K2_OFF}"
             if m == 'K2' else
             f"select ra, dec, kepmag, teff, radius, mass, dist from q1_q17_dr25_stellar where kepid={sid - sources.KEP_OFF}")
        rows = list(csv.DictReader(io.StringIO(requests.get('https://exoplanetarchive.ipac.caltech.edu/TAP/sync',
                                                          params={'query': q, 'format': 'csv'}, timeout=120).text)))
        v = [float(x) if x not in ('', None) else None for x in rows[0].values()]
        return (*v, None, None)
    except Exception:
        return None


def trapezoid(t, t0, depth, t14, t12):
    x = np.abs(t - t0)
    t12 = np.clip(t12, 1e-4, t14 / 2)
    y = np.zeros_like(t)
    full = x <= t14 / 2 - t12
    ing = (x > t14 / 2 - t12) & (x < t14 / 2)
    y[full] = -depth
    y[ing] = -depth * (t14 / 2 - x[ing]) / t12
    return y


def trapezoid_ss(t, t0, depth, t14, t12, exp):
    """Trapezoid integrated over the exposure time `exp` (supersampled), so that a transit of a few cadences is compared fairly."""
    exp = np.broadcast_to(np.asarray(exp, float), np.shape(t))
    off = (np.arange(5) - 2) / 5
    y = sum(trapezoid(t + o * exp, t0, depth, t14, t12) for o in off) / len(off)
    return np.where(exp > t14 / 40, y, trapezoid(t, t0, depth, t14, t12))


def _win(t, y, tc, half):
    a, b = np.searchsorted(t, [tc - half, tc + half])
    return t[a:b], y[a:b]


def box_depth(t, y, tc, dur_d):
    t, y = _win(t, y, tc, dur_d * 0.75 + 0.55)
    if t.size < 3:
        return None
    cad = np.median(np.diff(t))
    # central 70 % of the transit, but never narrower than the cadence allows (30-min Kepler/K2 data)
    half = min(max(dur_d * 0.35, 0.75 * cad), dur_d * 0.45)
    inn = np.abs(t - tc) < half
    out = (np.abs(t - tc) > dur_d * 0.75) & (np.abs(t - tc) < dur_d * 0.75 + 0.5)
    if inn.sum() < 2 or out.sum() < 10:
        return None
    sig = lcsearch.mad_std(y[out])
    d = np.median(y[out]) - np.mean(y[inn])
    return d, sig / np.sqrt(inn.sum()), inn.sum()


def coverage(t, tc, dur_d):
    """Fraction of the transit window that has data."""
    dt = np.median(np.diff(t[:2000])) if t.size > 1 else 0.01
    w, _ = _win(t, t, tc, dur_d / 2)
    n = w.size
    return min(1.0, n * dt / dur_d)


def refine_period(T, Y, P, T0, dur_d):
    """Scan a narrow grid of periods around a BLS period and keep the one whose folded light curve is deepest in transit.
    Over a baseline B a period error dP moves the last transit by (B/P) dP, so the grid spans +-REFINE_SPAN steps of
    P dur / B (one transit duration over the baseline) and is sampled at a quarter of that step, with nine epoch shifts of up to 0.75
    durations. The fold is anchored at the BLS epoch T0, where the BLS phase is right by construction, and B is the longest distance
    from T0 to either end of the data. Returns (P, t0), t0 being the epoch nearest to the middle of the data."""
    B = max(T[-1] - T0, T0 - T[0])
    if B < 20 * P:
        return P, T0
    res = P * dur_d / B
    Ps = P + res * np.linspace(-REFINE_SPAN, REFINE_SPAN, REFINE_TRIALS)
    base = np.median(Y)
    best = (-np.inf, P, T0)
    for q in Ps:
        for dt in dur_d * np.linspace(-0.75, 0.75, 9):      # the epoch moves together with the period
            x = ((T - T0 - dt) / q + 0.5) % 1 - 0.5
            inn = np.abs(x * q) < 0.35 * dur_d
            if inn.sum() >= 3 and base - np.mean(Y[inn]) > best[0]:
                best = (base - np.mean(Y[inn]), float(q), float(T0 + dt))
    q, t0 = best[1], best[2]
    return q, float(t0 + np.round((0.5 * (T[0] + T[-1]) - t0) / q) * q)


def gaia_checks(ra, dec, gaia_id, depth):
    out = {}
    try:
        from astroquery.gaia import Gaia
        q = f"""SELECT source_id, ra, dec, phot_g_mean_mag, ruwe, non_single_star, phot_variable_flag,
               DISTANCE(POINT({ra},{dec}), POINT(ra,dec))*3600 AS sep
               FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({ra},{dec},{60/3600}))"""
        tab = Gaia.launch_job(q).get_results()
        tab.sort('sep')
        # the TIC's own Gaia identifier names the target (a high-proper-motion star can lie farther away than a neighbour)
        mine = [i for i, sid in enumerate(tab['source_id']) if str(sid) == str(gaia_id)]
        tgt = tab[mine[0]] if mine else tab[0]
        out['ruwe'] = float(tgt['ruwe']) if tgt['ruwe'] is not None else None
        out['non_single_star'] = int(tgt['non_single_star'] or 0)
        out['variable'] = str(tgt['phot_variable_flag'])
        g0 = float(tgt['phot_g_mean_mag'])
        blends = []
        for r in tab:
            if r['source_id'] == tgt['source_id']:
                continue
            if r['phot_g_mean_mag'] is None or not np.isfinite(r['phot_g_mean_mag']):
                continue
            frac = 10 ** (-0.4 * (float(r['phot_g_mean_mag']) - g0))
            frac = frac / (1 + frac)
            need = depth / frac
            if need < 1.0:
                blends.append({'sep_arcsec': round(float(r['sep']), 1), 'dG': round(float(r['phot_g_mean_mag']) - g0, 2),
                               'eclipse_needed': round(need, 3)})
        out['possible_blends'] = blends
        eb = Gaia.launch_job(f"SELECT source_id FROM gaiadr3.vari_eclipsing_binary WHERE source_id={int(tgt['source_id'])}").get_results()
        out['gaia_eb'] = len(eb) > 0
    except Exception as e:
        out['error'] = repr(e)[:150]
    return out


def vet(con, tic, findings):
    cand = con.execute('SELECT kind, period, t0, dur_h, depth_ppm FROM candidate WHERE tic=?', (tic,)).fetchone()
    kind, P0, T0, dur_h, depth_ppm = cand
    snr0 = con.execute('SELECT snr FROM candidate WHERE tic=?', (tic,)).fetchone()[0] or 0
    star = get_star(con, tic)
    work = Path(findings).parent / 'data' / 'cache' / f'vet_{tic}'
    work.mkdir(parents=True, exist_ok=True)
    outdir = Path(findings) / sources.name_of(tic).replace(' ', '')
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        lcs = fetch_all(tic, work)
        data = {}
        for author, sector, path in lcs:
            s = series(author, path)
            if s is not None:
                data.setdefault(sector, {})[author] = s
        if not data:
            raise RuntimeError('no light curves')
        # one series per sector for the analysis: QLP (the search data), SPOC as the independent cross-check
        main = {s: (d.get('QLP') or d.get('TESS-SPOC') or d.get('SPOC') or d.get('Kepler') or d.get('K2SFF'))
                for s, d in data.items()}
        T = np.concatenate([v[0] for v in main.values()])
        Y = np.concatenate([v[1] for v in main.values()])
        # exposure time of every point: the median cadence of its sector (QLP is 30, 10 and 3.3 min in different sectors)
        EXPO = np.concatenate([np.full(len(v[0]), np.median(np.diff(v[0])) if len(v[0]) > 1 else 0.0) for v in main.values()])
        o = np.argsort(T); T, Y, EXPO = T[o], Y[o], EXPO[o]

        # ---- individual dips matching the candidate's depth (gap-edge and scattered-light dips excluded)
        events = []
        for s, (t, y) in sorted(main.items()):
            tg, dt, (yg,) = lcsearch.to_grid(t, y)
            cut = np.where(np.diff(t) > lcsearch.GAP_D)[0]
            segs_t = list(zip(np.r_[t[0], t[cut + 1]], np.r_[t[cut], t[-1]]))
            for e in lcsearch.single_events(tg, dt, yg, lcsearch.mad_std(yg)):
                c = lcsearch.characterise(e, tg, dt, yg, yg * 0, yg * 0, yg * 0, segs_t)
                # a star vetted by hand has no candidate depth yet: every dip enters, and the pairing below matches depths
                if (not depth_ppm or 0.5 < c['depth_ppm'] / depth_ppm < 2.0) and c['snr'] >= 7:
                    events.append((s, c['t0'], c['depth_ppm'], c['snr'], c['dur_h'], c['flags']))
        dur_d = (dur_h or 4) / 24
        notes = []

        # ---- period: every pair of similar dips proposes P = gap/n; keep the periods that the most dips share
        tol = max(0.06, dur_d / 4)
        start = time.monotonic()
        def check_time():
            if time.monotonic() - start > VET_LIMIT_S:
                raise VetTimeout(f'period search stopped after {VET_LIMIT_S} s ({len(events)} dips); vet by hand')
        ranked, halved, anchors = [], None, {}
        if kind == 'bls' and P0:
            # a BLS period error of 1e-4 d smears the stacked dips over a multi-year baseline: refit it before the epoch test
            Pr, Tr = refine_period(T, Y, P0, T0, dur_d)
            on = lambda ev, P_, T_: [e for e in ev if abs((((e[1] - T_) / P_ + 0.5) % 1 - 0.5) * P_) < tol]
            sup = on(events, Pr, Tr)
            if len(sup) < 2:
                sup = [(None, Tr, depth_ppm, snr0, dur_h, '')]
            ranked = [(Pr, sup)]; anchors[Pr] = Tr
            if P0 / 2 >= P_HALF_MIN:     # BLS often locks on 2 x P for short, faint planets: P/2 competes with P
                Ph, Th = refine_period(T, Y, P0 / 2, T0, dur_d)
                half = on(events, Ph, Th)
                ranked.append((Ph, half if len(half) >= 2 else [(None, Th, depth_ppm, snr0, dur_h, '')]))
                halved = Ph; anchors[Ph] = Th
        else:
            # the strongest clean dips only; for each pair, every period gap/n down to P_MIN is tested on all of them
            # at once, so that dips spread over years cost seconds, not hours
            ev = sorted((e for e in events if e[5] in ('', 'edge')), key=lambda e: -e[3])[:MAX_PAIR_DIPS]
            tv, dv = np.array([e[1] for e in ev]), np.array([e[2] for e in ev])
            trials = {}
            for i in range(len(ev)):
                for j in range(i + 1, len(ev)):
                    a_, b_ = ev[i], ev[j]
                    if not (0.65 < a_[2] / b_[2] < 1.55) or 'edge' in a_[5] + b_[5] and a_[5] == b_[5] == 'edge':
                        continue
                    g = abs(b_[1] - a_[1])
                    if g < P_MIN:
                        continue
                    check_time()
                    Ps = g / np.arange(1, int(g / P_MIN) + 1)
                    x = (((tv - a_[1]) / Ps[:, None] + 0.5) % 1 - 0.5) * Ps[:, None]
                    sup = (np.abs(x) < tol) & (0.65 < dv / a_[2]) & (dv / a_[2] < 1.55)
                    for P, row, n in zip(Ps.tolist(), sup, sup.sum(axis=1).tolist()):
                        key = round(P, 3)
                        if n >= 2 and (key not in trials or n > trials[key][1]):
                            trials[key] = (P, n, row)
            if trials:
                from scipy.stats import binom
                n_ev = max(1, len(events) - 2)
                v = list(trials.values())
                Ps, ns = np.array([x[0] for x in v]), np.array([x[1] for x in v])
                # -log10 chance that the extra dips line up with this period by accident
                sig = np.where(ns > 2, -np.log10(np.maximum(binom.sf(ns - 3, n_ev, np.minimum(1.0, 2 * tol / Ps)), 1e-300)), 0.0)
                ranked = [(v[k][0], [e for e, m in zip(ev, v[k][2]) if m]) for k in np.lexsort((-Ps, -sig))[:300]]
        if not ranked:
            raise RuntimeError(f'fewer than two consistent dips on re-analysis: {[(e[0], round(e[1], 3), round(e[2]), e[5]) for e in events]}')
        alive, tried = [], []
        dt_T = np.median(np.diff(T[:2000]))
        for P, sup in ranked:
            check_time()
            is_half = P == halved
            tt = np.array(sorted(e[1] for e in sup))
            ep = np.round((tt - tt[0]) / P)
            if P in anchors:      # the refit ephemeris is better than a line through a few matched, partly chance, dips
                t0 = anchors[P]
            else:
                P, t0 = np.polyfit(ep, tt, 1) if np.ptp(ep) > 0 else (P, tt[0])
            exp = float(np.median([e[2] for e in sup])) / 1e6
            n0, n1 = int(np.floor((T[0] - t0) / P)) - 1, int(np.ceil((T[-1] - t0) / P)) + 1
            ok, seen, misses, sw, sdw = True, 0, [], 0.0, 0.0
            par = {0: [], 1: []}
            tcs = t0 + np.arange(n0, n1 + 1) * P
            cov = (np.searchsorted(T, tcs + dur_d / 2) - np.searchsorted(T, tcs - dur_d / 2)) * dt_T / dur_d  # coverage()
            for tc in tcs[cov >= 0.6]:
                bd = box_depth(T, Y, tc, dur_d)
                if bd is None:
                    continue
                d, e, _ = bd
                seen += 1
                sw += 1 / e ** 2; sdw += d / e ** 2
                par[int(round((tc - t0) / P)) % 2].append((d, e))
                # a predicted transit is missing only if no dip is seen at all, the expected one would be obvious, and the
                # data are good enough to have detected it on its own (per-transit S/N of the candidate depth above the
                # dip-detection threshold); on a faint star single transits are below that limit and cannot be absent
                if exp / e >= DETECT_SNR and d / e < 2 and (exp - d) / e > 4:
                    ok = False; misses.append(round(tc, 3)); break
            if ok and sw > 0 and sdw / np.sqrt(sw) < STACK_SNR:
                # faint star: no single epoch is decisive, but the stack of all covered epochs shows no transit at all
                # (the candidate depth is no yardstick here: it is set by the strongest, noise-boosted dips)
                ok = False; misses.append('stack flat')
            if ok and is_half and min(len(v) for v in par.values()) >= 1:
                # P/2 must be one planet: the in-between epochs show a depth consistent with the others
                (d0, e0), (d1, e1) = [(np.average([d for d, _ in v], weights=[1 / e ** 2 for _, e in v]),
                                       1 / np.sqrt(sum(1 / e ** 2 for _, e in v))) for v in par.values()]
                if abs(d0 - d1) / np.hypot(e0, e1) > HALF_ODD_EVEN_SIGMA:
                    ok = False; misses.append('P/2 depths differ')
            tried.append((round(P, 4), len(sup), seen, misses))
            if ok:
                alive.append((P, t0, seen, sdw / np.sqrt(sw) if sw else 0.0, exp))
        if not alive:
            raise RuntimeError('every allowed period predicts transits that are absent in covered data; dips '
                               f'{[(e[0], round(e[1], 3), round(e[2]), e[5]) for e in events]}; tried {tried[:12]}')
        # the true period stacks the most signal over every covered epoch it predicts
        alive = [a for a in alive if np.isfinite(a[3])]
        top = max((a[3] for a in alive), default=0.0)
        if top <= 0:
            raise RuntimeError('no allowed period stacks a positive transit signal over its covered epochs')
        alive = [a for a in alive if a[3] >= 0.9 * top]
        P, t0, seen, _, exp = min(alive, key=lambda a: a[0])  # shortest period the data allow
        depth_ppm = exp * 1e6
        others = sorted({round(a[0], 4) for a in alive if abs(a[0] / P - round(a[0] / P)) > 0.01})
        if others:
            notes.append('periods that are not multiples of the chosen one also survive: ' + ', '.join(map(str, others[:8])))

        # ---- folded trapezoid fit
        ph = ((T - t0) / P + 0.5) % 1 - 0.5
        x = ph * P
        w = np.abs(x) < 3 * dur_d
        lo_b = [-dur_d / 2, 0, dur_d / 4, 1e-4]
        hi_b = [dur_d / 2, 1.0, dur_d * 3, dur_d * 1.5]
        p0 = np.clip([0, depth_ppm / 1e6, dur_d, dur_d / 6], np.array(lo_b) + 1e-9, np.array(hi_b) - 1e-9)
        fit = least_squares(lambda p: Y[w] - trapezoid_ss(x[w], p[0], p[1], p[2], p[3], EXPO[w]), p0, bounds=(lo_b, hi_b))
        dt0, dep, t14, t12 = fit.x
        t0 += dt0
        resid = Y[w] - trapezoid_ss(x[w], dt0, dep, t14, t12, EXPO[w])
        sig = lcsearch.mad_std(resid)
        # shape test: refit with the ingress fixed at a V (t12 = t14/2) and at a flat bottom (t12 = 0.1 t14), both integrated over
        # the cadence; the V is a reason to reject only if the data prefer it over the flat bottom
        def chi2_fixed(frac):
            lo2, hi2 = [-dur_d / 2, 0, dur_d / 4], [dur_d / 2, 1.0, dur_d * 3]
            q0 = np.clip([dt0, dep, t14], np.array(lo2) + 1e-9, np.array(hi2) - 1e-9)
            f2 = least_squares(lambda q: Y[w] - trapezoid_ss(x[w], q[0], q[1], q[2], frac * q[2] / 2, EXPO[w]), q0, bounds=(lo2, hi2))
            return np.sum(f2.fun ** 2) / sig ** 2
        v_pref = chi2_fixed(0.2) - chi2_fixed(1.0)   # > 0: the V fits better than the flat bottom
        n_in = np.sum(np.abs(x[w] - dt0) < t14 / 2)
        snr = dep / (sig / np.sqrt(max(n_in, 1)))
        # odd / even and secondary
        ep = np.round((T - t0) / P).astype(int)
        def d_at(mask, centre):
            inn = mask & (np.abs(((T - centre) / P + 0.5) % 1 - 0.5) * P < t14 * 0.35)
            out = mask & (np.abs(((T - centre) / P + 0.5) % 1 - 0.5) * P > t14 * 0.75) & \
                  (np.abs(((T - centre) / P + 0.5) % 1 - 0.5) * P < t14 * 0.75 + 0.5)
            if inn.sum() < 4 or out.sum() < 10:
                return np.nan, np.nan
            return np.median(Y[out]) - np.mean(Y[inn]), lcsearch.mad_std(Y[out]) / np.sqrt(inn.sum())
        # odd/even from fully covered transits only, at least two of each parity
        per = {0: [], 1: []}
        for n in sorted(set(ep)):
            tc = t0 + n * P
            if coverage(T, tc, t14) >= 0.8:
                bd = box_depth(T, Y, tc, t14)
                if bd:
                    per[n % 2].append(bd[:2])
        if len(per[0]) >= 2 and len(per[1]) >= 2:
            def wm(v):
                w = np.array([1 / e ** 2 for _, e in v]); d = np.array([d for d, _ in v])
                return np.sum(w * d) / np.sum(w), 1 / np.sqrt(np.sum(w))
            (d0, e0), (d1, e1) = wm(per[0]), wm(per[1])
            oe_sig = abs(d0 - d1) / np.hypot(e0, e1)
        else:
            oe_sig = np.nan
        # every fully covered epoch with baseline on both sides: is the transit seen in each on its own?
        epochs = []
        for n in sorted(set(ep)):
            tc = t0 + n * P
            if coverage(T, tc, t14) < 0.8:
                continue
            a_, b_ = np.searchsorted(T, [tc - t14 - 0.4, tc + t14 + 0.4])
            tw = T[a_:b_]
            if np.sum(tw < tc - 0.6 * t14) < 5 or np.sum(tw > tc + 0.6 * t14) < 5:
                continue        # one-sided baseline: at a gap edge, cannot be judged alone
            bd = box_depth(T, Y, tc, t14)
            if bd:
                yw = Y[a_:b_]
                pre, post = yw[(tw < tc - 0.6 * t14)], yw[(tw > tc + 0.6 * t14)]
                step = abs(np.median(pre) - np.median(post))
                step_e = np.hypot(lcsearch.mad_std(pre) / np.sqrt(pre.size), lcsearch.mad_std(post) / np.sqrt(post.size))
                epochs.append((int(n), bd[0], bd[1], step, step_e))
        sec, sec_e = d_at(np.ones_like(T, bool), t0 + P / 2)
        sec_sig = sec / sec_e if np.isfinite(sec) else np.nan
        # a real secondary eclipse is in all the data; one that sits in a single subset (e.g. the provisional QLP
        # timestamps of the newest sectors) is not: check the early and the late half separately
        half_sig = []
        for hm in (T < np.median(T), T >= np.median(T)):
            hs, he = d_at(hm, t0 + P / 2)
            if np.isfinite(hs):
                half_sig.append(hs / he)
        sec_consistent = not half_sig or min(half_sig) > SEC_HALF_SIGMA

        # ---- physical parameters
        ra, dec, tmag, teff, rad, mass, dist, gaia, contr = star if star else [None] * 9
        rad = rad or np.nan; mass = mass or np.nan
        # the TIC contamination ratio c = (flux of the neighbours) / (flux of the target) dilutes the transit, so that
        # depth_true = depth_obs (1 + c) (Stassun et al. 2019); a missing or negative ratio counts as no dilution
        c_ = contr if contr is not None and np.isfinite(contr) and contr > 0 else 0.0
        rp_obs = np.sqrt(dep) * rad * RSUN / REARTH
        rp = np.sqrt(dep * (1 + c_)) * rad * RSUN / REARTH
        a_au = ((mass * (P / 365.25) ** 2) ** (1 / 3)) if np.isfinite(mass) else np.nan
        rho_star = mass * MSUN / (4 / 3 * np.pi * (rad * RSUN) ** 3) / 1000 if np.isfinite(rad) else np.nan  # g/cc
        k = np.sqrt(dep)
        ar = P / (np.pi * t14) * (1 + k)  # central transit
        rho_transit = 3 * np.pi * ar ** 3 / (G * (P * 86400) ** 2) / 1000
        teq = teff * np.sqrt(rad * RSUN / (2 * a_au * 1.496e11)) * (1 - 0.3) ** 0.25 if np.isfinite(a_au) and teff else np.nan
        shape = t12 / (t14 / 2)

        gaia_res = gaia_checks(ra, dec, gaia, dep) if ra is not None else {}

        # ---- verdict
        reasons = []
        if rp > 24: reasons.append(f'radius {rp:.1f} Re is stellar-size')
        if dep > 0.25: reasons.append(f'depth {dep*100:.0f}% is too deep for a planet on any main-sequence star')
        if not np.isfinite(rp) and dep > 0.05: reasons.append(f'depth {dep*100:.1f}% with no stellar radius to rule out a stellar eclipse')
        resolved = t14 / np.median(EXPO[w]) >= V_MIN_CADENCES
        if shape > 0.8 and resolved and v_pref >= V_PREFER_SIGMA ** 2: reasons.append(f'V-shaped (ingress fraction {shape:.2f}, V preferred over flat bottom at {np.sqrt(v_pref):.1f} sigma)')
        elif shape > 0.8:
            notes.append(f'ingress fraction {shape:.2f}, V preferred over a flat bottom at {np.sqrt(max(v_pref, 0)):.1f} sigma, {t14 / np.median(EXPO[w]):.1f} cadences in transit: '
                         + ('the V does not beat a flat bottom' if resolved else 'ingress not resolved') + ', V-shape cut not applied')
        if np.isfinite(oe_sig) and oe_sig > 3: reasons.append(f'odd/even depth differ at {oe_sig:.1f} sigma')
        if np.isfinite(sec_sig) and sec_sig > 3:
            if sec_consistent: reasons.append(f'secondary eclipse at {sec_sig:.1f} sigma')
            else: notes.append(f'secondary at {sec_sig:.1f} sigma not confirmed in both halves of the data ({", ".join(f"{h:.1f}" for h in half_sig)} sigma): not used')
        if gaia_res.get('gaia_eb'): reasons.append('Gaia DR3 eclipsing binary')
        if gaia_res.get('ruwe') and gaia_res['ruwe'] > 1.4: reasons.append(f'Gaia RUWE {gaia_res["ruwe"]:.2f} (unresolved companion)')
        if snr < 10: reasons.append(f'folded SNR only {snr:.1f}')
        steps = sum(st > 0.5 * dep and st > 3 * se for _, _, _, st, se in epochs)
        if epochs and steps >= max(1, len(epochs) / 2):
            reasons.append(f'baseline differs before and after the dip in {steps} of {len(epochs)} epochs: a step or ramp, not a transit')
        # only epochs whose own data could show the fitted depth (per-transit S/N >= DETECT_SNR) can be judged alone; on a
        # faint star the stack (folded S/N above) carries the evidence instead
        epochs = [x[:3] for x in epochs]
        powered = [x for x in epochs if dep / x[2] >= DETECT_SNR]
        clear = [(n, d, e) for n, d, e in powered if d / e >= 5 and 0.4 < d / dep < 2.5]
        if len(powered) >= 2 and len(clear) < 2:
            reasons.append(f'only {len(clear)} of {len(powered)} fully covered epochs show the transit on its own')
        elif powered and sum(d / e >= 3 for _, d, e in powered) < 0.6 * len(powered):
            reasons.append(f'transit absent in {sum(d / e < 3 for _, d, e in powered)} of {len(powered)} fully covered epochs')
        if len(clear) >= 2:
            d_ = np.array([d for _, d, _ in clear]); e_ = np.array([e for _, _, e in clear])
            m_ = np.sum(d_ / e_ ** 2) / np.sum(1 / e_ ** 2)
            chi = np.sum(((d_ - m_) / e_) ** 2)
            from scipy.stats import chi2 as _chi2
            if _chi2.sf(chi, len(d_) - 1) < 1e-6 and d_.max() / d_.min() > 2.5:
                reasons.append(f'depth changes between epochs ({d_.min()*1e6:.0f}-{d_.max()*1e6:.0f} ppm): aperture-dependent dilution of a blend')
        if np.isfinite(rho_star) and rho_transit < 0.05 * rho_star:
            reasons.append(f'transit far too long for the star (density {rho_transit:.3f} vs {rho_star:.2f} g/cc)')
        for (otic, oP, oT0) in con.execute('SELECT tic, period, t0 FROM candidate WHERE tic != ? AND period > 0', (tic,)).fetchall():
            if abs(oP / P - 1) < 0.003 and oT0 is not None and abs((((oT0 - t0) / P + 0.5) % 1 - 0.5) * P) < 0.3:
                reasons.append(f'same period and epoch as {sources.name_of(otic)}: sector-synchronised systematic')
                break
        if np.isfinite(rho_star) and 0.05 * rho_star <= rho_transit < 0.2 * rho_star:
            reasons.append(f'transit too long for the star at this period (density {rho_transit:.2f} vs {rho_star:.2f} g/cc)')
        if not np.isfinite(rho_star) and rho_transit < 0.05:
            reasons.append(f'transit duration implies a stellar density of only {rho_transit:.3f} g/cc (giant host or systematic)')
        # coherent stellar variability faster than the transit (pulsator): strongest out-of-transit oscillation vs depth
        oot = np.abs(T - (t0 + np.round((T - t0) / P) * P)) > t14
        if oot.sum() > 200 and 1.5 / t14 < 50:
            from astropy.timeseries import LombScargle
            ls = LombScargle(T[oot], Y[oot])
            fq, pw = ls.autopower(minimum_frequency=1.5 / t14, maximum_frequency=50, samples_per_peak=5)
            fb = fq[np.argmax(pw)]
            mdl = ls.model(np.linspace(0, 1 / fb, 50), fb)
            amp = (mdl.max() - mdl.min()) / 2
            if amp > 0.5 * dep:
                reasons.append(f'star oscillates with semi-amplitude {amp*1e6:.0f} ppm at {24/fb:.2f} h, comparable to the {dep*1e6:.0f} ppm dip: pulsator, not a transit')
        # dips of the same depth recurring away from the ephemeris (rotational/spot or other variability)
        on_eph = [e for e in events if abs((((e[1] - t0) / P + 0.5) % 1 - 0.5) * P) <= t14]
        n_off = len(events) - len(on_eph)
        if n_off >= 10 and n_off >= 3 * max(1, len(on_eph)):
            reasons.append(f'{n_off} dips of similar depth occur away from the ephemeris (against {len(on_eph)} on it): stellar variability, not a transit')
        if len(epochs) < 2:
            reasons.append(f'only {len(epochs)} transit(s) covered with baseline on both sides; cannot be confirmed')
        if any(99 <= s <= 104 for s in main):
            notes.append('QLP Sectors 99-104 await timestamp correction; mid-times from them are provisional')
        status = 'rejected' if reasons else 'ready'
        score = float(snr) * (0.5 if notes else 1.0)

        # ---- keep the relevant data: transit windows only
        n_all = np.round((T - t0) / P)
        keep = np.abs(T - (t0 + n_all * P)) < 1.0
        np.savez_compressed(outdir / 'transit_windows.npz', t=T[keep], f=Y[keep], period=P, t0=t0)

        # ---- plots: each covered transit, then the fold
        tr = sorted({int(n) for n in n_all[np.abs(T - (t0 + n_all * P)) < t14 / 2]})
        shown = tr[:MAX_PANELS]   # stars observed for years have hundreds of transits; the fold shows them all
        ncol = min(4, max(1, len(shown)))
        nrow = int(np.ceil(len(shown) / ncol)) + 1
        fig = plt.figure(figsize=(4 * ncol, 3 * nrow))
        for i, n in enumerate(shown):
            ax = fig.add_subplot(nrow, ncol, i + 1)
            tc = t0 + n * P
            m = np.abs(T - tc) < 1.0
            ax.plot((T[m] - tc) * 24, Y[m] * 1e6, 'o', ms=3, color='0.15')
            ax.axvspan(-t14 * 12, t14 * 12, color='C0', alpha=0.12)
            ax.set_title(f'epoch {n}  BTJD {tc:.3f}', fontsize=9)
            ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('ppm')
        ax = fig.add_subplot(nrow, 1, nrow)
        m = np.abs(x) < 3 * dur_d
        ax.plot(x[m] * 24, Y[m] * 1e6, 'o', ms=2.5, color='0.15')
        xs = np.linspace(-3 * dur_d, 3 * dur_d, 600)
        ax.plot(xs * 24, trapezoid(xs, dt0, dep, t14, t12) * 1e6, color='C3', lw=2)
        ax.set_xlabel('hours from mid-transit (folded)'); ax.set_ylabel('ppm')
        ax.set_title(f'{sources.name_of(tic)}  P = {P:.5f} d  depth {dep*1e6:.0f} ppm  T14 {t14*24:.2f} h  Rp {rp:.2f} Re', fontsize=10)
        fig.tight_layout(); fig.savefig(outdir / 'transits.png', dpi=110); plt.close(fig)

        summary = {'tic': tic, 'status': status, 'reasons': reasons, 'notes': notes, 'period_d': P,
                   't0_bjd': t0 + BTJD, 'depth_ppm': dep * 1e6, 't14_h': t14 * 24, 'ingress_frac': shape,
                   'snr': snr, 'n_transits_seen': len(tr), 'sectors': sorted(main), 'rp_earth': rp, 'rp_earth_uncorrected': rp_obs, 'a_au': a_au,
                   'teq_k': teq, 'rho_star': rho_star, 'rho_transit_b0': rho_transit, 'odd_even_sigma': oe_sig,
                   'secondary_sigma': sec_sig, 'gaia': gaia_res,
                   'epochs': [{'epoch': n, 'depth_ppm': d * 1e6, 'err_ppm': e * 1e6} for n, d, e in epochs],
                   'star': {'tmag': tmag, 'teff': teff, 'rad': rad, 'mass': mass, 'dist_pc': dist, 'contratio': contr,
                            'ra': ra, 'dec': dec},
                   'matched_dips': [{'sector': e[0], 't0_btjd': e[1], 'depth_ppm': e[2], 'snr': e[3]} for e in events],
                   'surviving_periods': [a[0] for a in alive], 'vetted': time.strftime('%Y-%m-%d %H:%M')}
        (outdir / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
        (outdir / 'report.md').write_text(report_md(summary), encoding='utf-8')
        con.execute('UPDATE candidate SET status=?, period=?, t0=?, depth_ppm=?, dur_h=?, snr=?, score=?, notes=?, updated=? '
                    'WHERE tic=?', (status, P, t0, dep * 1e6, t14 * 24, snr, score,
                                    '; '.join(reasons + notes), time.strftime('%Y-%m-%d %H:%M'), tic))
        con.commit()
    except VetTimeout as e:
        con.execute("UPDATE candidate SET status='vet_error', notes=?, updated=? WHERE tic=?",
                    (str(e), time.strftime('%Y-%m-%d %H:%M'), tic))
        con.commit()
        (outdir / 'report.md').write_text(f'# TIC {tic}\n\nNot vetted: {e}\n', encoding='utf-8')
    except RuntimeError as e:
        con.execute("UPDATE candidate SET status='rejected', notes=?, updated=? WHERE tic=?",
                    (str(e), time.strftime('%Y-%m-%d %H:%M'), tic))
        con.commit()
        (outdir / 'report.md').write_text(f'# TIC {tic}\n\nRejected during vetting: {e}\n', encoding='utf-8')
    finally:
        shutil.rmtree(work, ignore_errors=True)


def report_md(s):
    st = s['star']
    L = [f"# TIC {s['tic']} — {s['status'].upper()}", '',
         f"Period {s['period_d']:.5f} d, mid-transit {s['t0_bjd']:.4f} BJD, depth {s['depth_ppm']:.0f} ppm, "
         f"duration {s['t14_h']:.2f} h, folded SNR {s['snr']:.1f}, {s['n_transits_seen']} transits seen in sectors "
         f"{', '.join(map(str, s['sectors']))}.", '',
         f"Host: Tmag {st['tmag']}, Teff {st['teff']} K, R {st['rad']} Rsun, M {st['mass']} Msun, d {st['dist_pc']} pc, "
         f"TIC contamination ratio {st['contratio']}.", '',
         f"Planet: Rp {s['rp_earth']:.2f} Re after the dilution correction ({s['rp_earth_uncorrected']:.2f} Re without it), a {s['a_au']:.4f} AU, Teq {s['teq_k']:.0f} K.", '',
         f"Tests: ingress fraction {s['ingress_frac']:.2f}, odd/even {s['odd_even_sigma']:.1f} sigma, "
         f"secondary {s['secondary_sigma']:.1f} sigma, stellar density from TIC {s['rho_star']:.2f} g/cc versus "
         f"{s['rho_transit_b0']:.2f} g/cc from the transit (central-transit assumption).", '',
         f"Gaia: {json.dumps(s['gaia'])}", '']
    if s['reasons']:
        L += ['Rejection reasons: ' + '; '.join(s['reasons']), '']
    if s['notes']:
        L += ['Notes: ' + '; '.join(s['notes']), '']
    L += ['![transits](transits.png)', '']
    return '\n'.join(L)
