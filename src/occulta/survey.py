"""Transit hunt over TESS (FFI and 2-min), K2 and Kepler: stream light curves of every star type,
search them, keep only the detections, purge the raw files, and promote
cross-sector signals that no catalogue knows about to candidates.

    occulta survey run --minutes 55 --mode both   # one budgeted, resumable work session
    occulta survey status                         # progress and candidates
    occulta planets vet TIC                       # vet one star by hand
"""
import argparse, csv, io, json, os, shutil, sqlite3, sys, time, traceback
import concurrent.futures as cf
from pathlib import Path

from . import paths

import numpy as np
import requests

from . import sources as sources

ROOT = paths.HOME
DATA, RAW, CACHE = ROOT / 'data', ROOT / 'data' / 'raw', ROOT / 'data' / 'cache'
FIND, LOGS = ROOT / 'findings', ROOT / 'logs'
DB = DATA / 'hunt.db'
CFG = paths.config()
LOCK = DATA / 'run.lock'
BTJD = 2457000.0


def log(*a):
    msg = time.strftime('%Y-%m-%d %H:%M:%S ') + ' '.join(str(x) for x in a)
    print(msg, flush=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    with open(LOGS / 'survey.log', 'a', encoding='utf-8') as fh:
        fh.write(msg + '\n')


def db(path=None):
    con = sqlite3.connect(path or DB, timeout=60)
    con.executescript("""
    CREATE TABLE IF NOT EXISTS tic(tic INTEGER PRIMARY KEY, ra REAL, dec REAL, tmag REAL, teff REAL, rad REAL,
        mass REAL, dist REAL, gaia TEXT, contratio REAL);
    CREATE TABLE IF NOT EXISTS tic_band(band INTEGER PRIMARY KEY, n INTEGER);
    CREATE TABLE IF NOT EXISTS known(tic INTEGER, src TEXT, period REAL, epoch REAL, dur_h REAL, depth_ppm REAL);
    CREATE INDEX IF NOT EXISTS known_tic ON known(tic);
    CREATE TABLE IF NOT EXISTS sector(sector INTEGER PRIMARY KEY, status TEXT, n_targets INTEGER, version TEXT,
        started TEXT, finished TEXT);
    CREATE TABLE IF NOT EXISTS target(sector INTEGER, tic INTEGER, status INTEGER DEFAULT 0, url TEXT, PRIMARY KEY(sector, tic));
    CREATE TABLE IF NOT EXISTS job(job TEXT PRIMARY KEY, status TEXT, n_targets INTEGER, started TEXT, finished TEXT);
    CREATE TABLE IF NOT EXISTS star(sid INTEGER PRIMARY KEY, mag REAL, teff REAL, rad REAL, mass REAL, ra REAL, dec REAL, dist REAL);
    CREATE TABLE IF NOT EXISTS star_sector(tic INTEGER, sector INTEGER, cam INTEGER, ccd INTEGER, npts INTEGER,
        t_start REAL, t_end REAL, noise_ppm REAL, noise1h_ppm REAL, inv_snr REAL, PRIMARY KEY(tic, sector));
    CREATE TABLE IF NOT EXISTS event(tic INTEGER, sector INTEGER, cam INTEGER, ccd INTEGER, t0 REAL, dur_h REAL,
        depth_ppm REAL, snr REAL, n_in INTEGER, vratio REAL, flags TEXT);
    CREATE INDEX IF NOT EXISTS event_tic ON event(tic);
    CREATE TABLE IF NOT EXISTS bls(tic INTEGER, sector INTEGER, period REAL, t0 REAL, dur_h REAL, depth_ppm REAL,
        snr REAL, sde REAL, ntr INTEGER, oddeven_sig REAL, sec_depth_ppm REAL, PRIMARY KEY(tic, sector));
    CREATE TABLE IF NOT EXISTS candidate(tic INTEGER PRIMARY KEY, kind TEXT, period REAL, t0 REAL, dur_h REAL,
        depth_ppm REAL, snr REAL, n_events INTEGER, sectors TEXT, status TEXT, score REAL, notes TEXT,
        created TEXT, updated TEXT);
    """)
    return con


def dir_gb(p):
    return sum(f.stat().st_size for f in Path(p).rglob('*') if f.is_file()) / 1e9


# ---------------------------------------------------------------- catalogues
TIERS = [10.0, 11.0, 12.0]


def band_key(tmag, band):
    return band if tmag == 10.0 else int(round(tmag * 10)) * 1000 + band + 90


def build_catalog(con, deadline, tmag):
    """TIC stars with TESS magnitude in the tier (previous tier, tmag], all sky, in 1-degree bands (resumable)."""
    from astroquery.mast import Catalogs
    lo = max([-5.0] + [t for t in TIERS if t < tmag])
    done = {r[0] for r in con.execute('SELECT band FROM tic_band')}
    for band in range(-90, 90):
        key = band_key(tmag, band)
        if key in done:
            continue
        if time.time() > deadline:
            return False
        for attempt in range(3):
            try:
                t = Catalogs.query_criteria(catalog='Tic', Tmag=[lo + 1e-4 if lo > -5 else -5, tmag],
                                            dec=[band, band + 1], objType='STAR')
                break
            except Exception as e:
                log('TIC band', band, 'retry', attempt, repr(e)[:120]); time.sleep(20)
        else:
            return False
        rows = []
        for r in t:
            if str(r['disposition']) in ('SPLIT', 'DUPLICATE', 'ARTIFACT'):
                continue
            f = lambda k: float(r[k]) if np.isfinite(r[k]) else None
            rows.append((int(r['ID']), f('ra'), f('dec'), f('Tmag'), f('Teff'), f('rad'), f('mass'), f('d'),
                         str(r['GAIA']), f('contratio')))
        con.executemany('INSERT OR REPLACE INTO tic VALUES (?,?,?,?,?,?,?,?,?,?)', rows)
        con.execute('INSERT INTO tic_band VALUES (?,?)', (key, len(rows)))
        con.commit()
        log(f'TIC tier {tmag} band', band, len(rows), 'stars')
    return True


def refresh_known(con):
    """TOIs, CTOIs, confirmed planets, Kepler KOIs, K2 candidates and EB catalogues; refreshed weekly.
    Epochs are stored as full BJD; star ids use the sources.py namespace."""
    stamp = CACHE / 'known.stamp'
    if stamp.exists() and time.time() - stamp.stat().st_mtime < 7 * 86400:
        return
    rows = []
    def num(x):
        try: return float(x)
        except Exception: return None
    def tap(q):
        r = requests.get('https://exoplanetarchive.ipac.caltech.edu/TAP/sync', params={'query': q, 'format': 'csv'}, timeout=600)
        return list(csv.DictReader(io.StringIO(r.text)))
    toi = requests.get('https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv', timeout=300).text
    for r in csv.DictReader(io.StringIO(toi)):
        rows.append((int(r['TIC ID']), 'TOI ' + r['TOI'] + ' ' + r['TFOPWG Disposition'], num(r['Period (days)']),
                     num(r['Epoch (BJD)']), num(r['Duration (hours)']), num(r['Depth (ppm)'])))
    ctoi = requests.get('https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv', timeout=300).text
    for r in csv.DictReader(io.StringIO(ctoi)):
        rows.append((int(r['TIC ID']), 'CTOI ' + r['CTOI'], num(r['Period (days)']), num(r['Transit Epoch (BJD)']),
                     num(r['Duration (hrs)']), num(r['Depth ppm'])))
    for r in tap("select tic_id,pl_name,pl_orbper,pl_tranmid,pl_trandur from pscomppars where tic_id is not null"):
        try: tid = int(str(r['tic_id']).replace('TIC', '').strip())
        except Exception: continue
        rows.append((tid, 'Confirmed ' + r['pl_name'], num(r['pl_orbper']), num(r['pl_tranmid']), num(r['pl_trandur']), None))
    for r in tap("select kepid,kepoi_name,koi_disposition,koi_period,koi_time0bk,koi_duration,koi_depth from cumulative"):
        e = num(r['koi_time0bk'])
        disp = r['koi_disposition'].replace('FALSE POSITIVE', 'FP')
        rows.append((sources.KEP_OFF + int(r['kepid']), f"KOI {r['kepoi_name']} {disp}", num(r['koi_period']),
                     e + 2454833.0 if e else None, num(r['koi_duration']), num(r['koi_depth'])))
    for r in tap("select epic_hostname,pl_name,disposition,pl_orbper,pl_tranmid,pl_trandur from k2pandc"):
        try: eid = int(str(r['epic_hostname']).replace('EPIC', '').strip())
        except Exception: continue
        disp = r['disposition'].replace('FALSE POSITIVE', 'FP')
        rows.append((sources.K2_OFF + eid, f"K2 {r['pl_name']} {disp}", num(r['pl_orbper']), num(r['pl_tranmid']),
                     num(r['pl_trandur']), None))
    from astroquery.vizier import Vizier
    for cat, col, off in (('J/ApJS/258/16/tess-ebs', 'TIC', 0), ('J/AJ/151/68/catalog', 'KIC', sources.KEP_OFF)):
        try:
            eb = Vizier(row_limit=-1).get_catalogs(cat)[0]
            pcol = next(c for c in ('Per', 'Period', 'P') if c in eb.colnames)
            ecol = next((c for c in ('BJD0', 'T0', 'Epoch') if c in eb.colnames), None)
            for r in eb:
                e = num(r[ecol]) if ecol else None
                if e is not None and e < 2e6:
                    e += 2400000.0 if e > 50000 else 2454833.0
                rows.append((off + int(r[col]), 'EB catalogue', num(r[pcol]), e, None, None))
        except Exception as e:
            log('EB catalogue fetch failed', cat, repr(e)[:120])
    con.execute('DELETE FROM known')
    con.executemany('INSERT INTO known VALUES (?,?,?,?,?,?)', rows)
    con.commit()
    stamp.touch()
    log('known objects refreshed:', len(rows))


# ---------------------------------------------------------------- jobs
def jobs():
    """The search plan, in priority order: (job id, mission, tier, unit)."""
    for step in CFG['plan']:
        m, tier = step['mission'], step.get('tmag')
        for u in sources.units(m):
            yield f'{m}:{tier}:{u}', m, tier, u


def next_job(con):
    status = dict(con.execute('SELECT job, status FROM job'))
    for jid, m, tier, u in jobs():
        if status.get(jid) in (None, 'targets', 'running'):
            return jid, m, tier, u
    return None


def prepare_job(con, jid, m, tier, u, deadline):
    if m == 'QLP' and tier is not None and not build_catalog(con, deadline, tier):
        return False
    kep = None
    if m == 'Kepler':
        kep = [int(r[0]) - sources.KEP_OFF for r in con.execute('SELECT sid FROM star WHERE sid >= ?', (sources.KEP_OFF,))]
        if not kep:
            rows = sources.kepler_targets()
            val = lambda r, k: float(r[k]) if r.get(k) not in ('', None) else None
            con.executemany('INSERT OR REPLACE INTO star VALUES (?,?,?,?,?,?,?,?)',
                            [(sources.KEP_OFF + int(r['kepid']), *(val(r, k) for k in
                              ('kepmag', 'teff', 'radius', 'mass', 'ra', 'dec', 'dist'))) for r in rows])
            con.commit()
            kep = [int(r['kepid']) for r in rows]
    tg = sources.unit_targets(u, kep)
    if tg is None and m == 'ZTF':
        con.execute("INSERT OR REPLACE INTO job VALUES (?,?,?,?,?)", (jid, 'missing', 0, None, None))
        con.commit(); log(jid, 'ZTF has no target list; search a position with `occulta ztf RA DEC`'); return True
    if tg is None:
        con.execute("INSERT OR REPLACE INTO job VALUES (?,?,?,?,?)", (jid, 'missing', 0, None, None))
        con.commit(); log(jid, 'no data at the archive'); return True
    if m == 'QLP':
        con.execute('CREATE TEMP TABLE IF NOT EXISTS ids(tic INTEGER PRIMARY KEY)')
        con.execute('DELETE FROM ids')
        con.executemany('INSERT OR IGNORE INTO ids VALUES (?)', ((i,) for i in tg))
        if tier is None:   # final pass: everything QLP has that earlier passes did not take
            sel = [x[0] for x in con.execute(
                'SELECT ids.tic FROM ids WHERE ids.tic NOT IN (SELECT tic FROM target WHERE sector=?)', (u,))]
        else:
            lo = max([-5.0] + [t for t in TIERS if t < tier])
            sel = [x[0] for x in con.execute('SELECT ids.tic FROM ids JOIN tic USING(tic) WHERE tic.tmag > ? '
                                              'AND tic.tmag <= ? ORDER BY tic.tmag', (lo, tier))]
        tg = {t: None for t in sel}
    con.executemany('INSERT OR IGNORE INTO target(sector, tic, url) VALUES (?,?,?)', ((u, t, v) for t, v in tg.items()))
    con.execute("INSERT OR REPLACE INTO job VALUES (?,?,?,?,?)", (jid, 'targets', len(tg), time.strftime('%Y-%m-%d %H:%M'), None))
    con.commit()
    log(jid, f'({sources.unit_name(u)}):', len(tg), 'light curves queued')
    return True


SESSION = requests.Session()


def fetch(u, sid, url):
    out = RAW / f'u{u}_{sid}.fits'
    urls = [url] if url else sources.url_for(u, sid)
    for v in urls:
        for attempt in range(3):
            try:
                r = SESSION.get(v, timeout=180)
                if r.status_code == 200 and r.content[:6] == b'SIMPLE':
                    out.write_bytes(r.content)
                    return str(out)
                if r.status_code in (403, 404):
                    break
            except Exception:
                time.sleep(5)
    return None


def store(con, res):
    if 'error' in res:
        con.execute('UPDATE target SET status=2 WHERE sector=? AND tic=?', (res['sector'], res['tic']))
        return
    t, s = res['tic'], res['sector']
    con.execute('INSERT OR REPLACE INTO star_sector VALUES (?,?,?,?,?,?,?,?,?,?)',
                (t, s, res['cam'], res['ccd'], res['npts'], res['t_start'], res['t_end'], res['noise_ppm'],
                 res['noise1h_ppm'], res['inv_snr']))
    con.execute('DELETE FROM event WHERE tic=? AND sector=?', (t, s))
    con.executemany('INSERT INTO event VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                    [(t, s, res['cam'], res['ccd'], e['t0'], e['dur_h'], e['depth_ppm'], e['snr'], e['n_in'],
                      e['vratio'], e['flags']) for e in res['events']])
    b = res['bls']
    if b and b['sde'] >= CFG['bls_sde'] and b['snr'] >= CFG['bls_snr']:
        con.execute('INSERT OR REPLACE INTO bls VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                    (t, s, b['period'], b['t0'], b['dur_h'], b['depth_ppm'], b['snr'], b['sde'], b['ntr'],
                     b['oddeven_sig'], b['sec_depth_ppm']))
    con.execute('UPDATE target SET status=1 WHERE sector=? AND tic=?', (s, t))


def work_job(con, jid, m, u, deadline):
    from . import lightcurves as lcsearch
    fmt = sources.fmt_of_unit(u)
    do_bls = m != 'Kepler'     # Kepler's periodic signals were exhausted by its own pipeline; we look for single events
    con.execute("UPDATE job SET status='running' WHERE job=?", (jid,)); con.commit()
    todo = con.execute('SELECT tic, url FROM target WHERE sector=? AND status=0', (u,)).fetchall()
    log(jid, len(todo), 'light curves to go')
    n = 0
    with cf.ThreadPoolExecutor(CFG['download_threads']) as dl, cf.ProcessPoolExecutor(CFG['workers']) as pp:
        pending_dl, pending_pp = {}, {}
        it = iter(todo)
        exhausted = False
        while True:
            while (not exhausted and len(pending_dl) < CFG['download_threads'] * 2
                   and len(pending_pp) < CFG['workers'] * 6 and time.time() < deadline
                   and dir_gb(RAW) < CFG['raw_cap_gb']):
                nxt = next(it, None)
                if nxt is None:
                    exhausted = True; break
                pending_dl[dl.submit(fetch, u, nxt[0], nxt[1])] = nxt[0]
            if not pending_dl and not pending_pp:
                break
            done, _ = cf.wait(list(pending_dl) + list(pending_pp), timeout=30, return_when=cf.FIRST_COMPLETED)
            for f in done:
                if f in pending_dl:
                    sid = pending_dl.pop(f)
                    path = f.result()
                    if path is None:   # not observed in this unit, or missing at the archive
                        con.execute('UPDATE target SET status=3 WHERE sector=? AND tic=?', (u, sid))
                        continue
                    pending_pp[pp.submit(lcsearch.process, path, sid, u, fmt, do_bls)] = path
                else:
                    path = pending_pp.pop(f)
                    try:
                        store(con, f.result())
                    except Exception:
                        log('store failed', traceback.format_exc()[-300:])
                    try: os.remove(path)          # purge raw data as soon as it is searched
                    except OSError: pass
                    n += 1
                    if n % 500 == 0:
                        con.commit(); log(jid, n, 'processed this session')
        con.commit()
    left = con.execute('SELECT COUNT(*) FROM target WHERE sector=? AND status=0', (u,)).fetchone()[0]
    if left == 0:
        mark_systematics(con, u)
        con.execute("UPDATE job SET status='done', finished=? WHERE job=?", (time.strftime('%Y-%m-%d %H:%M'), jid))
        con.commit()
        log(jid, 'complete')
    return left == 0


def mark_systematics(con, s):
    """Events shared by many stars on the same CCD at the same time are instrumental."""
    ev = con.execute('SELECT rowid, cam, ccd, t0 FROM event WHERE sector=?', (s,)).fetchall()
    by = {}
    for rid, cam, ccd, t0 in ev:
        by.setdefault((cam, ccd), []).append((t0, rid))
    bad = []
    for k, lst in by.items():
        lst.sort()
        ts = np.array([x[0] for x in lst])
        for i, (t0, rid) in enumerate(lst):
            n = np.searchsorted(ts, t0 + 0.1) - np.searchsorted(ts, t0 - 0.1)
            if n >= CFG['common_mode_n']:
                bad.append(rid)
    for rid in bad:
        con.execute("UPDATE event SET flags = CASE WHEN flags='' THEN 'common' ELSE flags||',common' END WHERE rowid=?", (rid,))
    con.commit()
    log('sector', s, 'common-mode events flagged:', len(bad))


# ---------------------------------------------------------------- candidates
def known_match(con, tic, t0_btjd, tol_d=0.3):
    """'host' if the star carries a known object, 'match' if this event is one of its transits, 'eb' if a catalogued EB."""
    rows = con.execute('SELECT src, period, epoch FROM known WHERE tic=?', (tic,)).fetchall()
    if not rows:
        return None
    for src, P, E in rows:
        if src.startswith('EB catalogue') or src.split()[-1:] in (['FP'], ['EB'], ['FA']):
            return 'eb'
        if P and E and P > 0:
            e = E - BTJD if E > 2e6 else E
            ph = ((t0_btjd - e) / P + 0.5) % 1 - 0.5
            if abs(ph * P) < tol_d:
                return 'match'
    return 'host'


def promote(con):
    """Cross-sector stitching: stars with repeated, mutually consistent, clean dips -> candidates."""
    c = CFG
    stars = [r[0] for r in con.execute(
        'SELECT tic FROM event WHERE snr >= ? AND (flags = "" OR flags = "edge") GROUP BY tic HAVING COUNT(*) >= 2',
        (c['event_snr'],))]
    new = 0
    for tic in stars:
        if con.execute('SELECT 1 FROM candidate WHERE tic=?', (tic,)).fetchone():
            continue
        ev = con.execute('SELECT e.sector, e.t0, e.dur_h, e.depth_ppm, e.snr, e.flags, ss.inv_snr FROM event e '
                         'JOIN star_sector ss USING(tic, sector) WHERE e.tic=? AND e.snr >= ? AND (e.flags="" OR e.flags="edge")',
                         (tic, c['event_snr'])).fetchall()
        ev = [e for e in ev if e[4] >= c['inv_ratio'] * e[6] and known_match(con, tic, e[1]) not in ('match', 'eb')]
        clean = [e for e in ev if e[5] == '']
        if len(ev) < 2 or len(clean) < 1:
            continue
        # group events of similar depth and duration
        ev.sort(key=lambda e: -e[4])
        ref = ev[0]
        grp = [e for e in ev if 0.6 < e[3] / ref[3] < 1.6 and 0.5 <= e[2] / ref[2] <= 2.0]
        # the same transit seen in two light curves (QLP and SPOC of one sector) counts once: keep the stronger
        uniq = []
        for e in grp:                                   # grp is sorted by decreasing SNR
            if all(abs(e[1] - u[1]) >= 0.5 for u in uniq):
                uniq.append(e)
        grp = uniq
        if len(grp) < 2:
            continue
        times = sorted(e[1] for e in grp)
        snr = float(np.sqrt(sum(e[4] ** 2 for e in grp)))
        host = known_match(con, tic, times[0])
        con.execute('INSERT INTO candidate VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                    (tic, 'multi', None, times[0], float(np.median([e[2] for e in grp])),
                     float(np.median([e[3] for e in grp])), snr, len(grp), ','.join(str(e[0]) for e in grp),
                     'new', None, ('known host; ' if host == 'host' else ''), time.strftime('%Y-%m-%d %H:%M'), None))
        new += 1
    # short-period periodic signals from per-sector BLS
    for tic, s, P, t0, dur, dep, snr, sde, ntr, oe in con.execute(
            'SELECT tic, sector, period, t0, dur_h, depth_ppm, snr, sde, ntr, oddeven_sig FROM bls '
            'WHERE ntr >= 3 AND oddeven_sig < 3 AND snr >= ? AND sde >= ?', (c['bls_snr'], c['bls_sde'])).fetchall():
        if con.execute('SELECT 1 FROM candidate WHERE tic=?', (tic,)).fetchone():
            continue
        if known_match(con, tic, t0) in ('match', 'eb'):
            continue
        con.execute('INSERT INTO candidate VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                    (tic, 'bls', P, t0, dur, dep, snr, ntr, str(s), 'new', None, '', time.strftime('%Y-%m-%d %H:%M'), None))
        new += 1
    con.commit()
    if new:
        log('promoted', new, 'new candidates')
    return new


def vet_new(con, deadline):
    from . import planets as vet
    for (tic,) in con.execute("SELECT tic FROM candidate WHERE status='new' ORDER BY snr DESC").fetchall():
        if time.time() > deadline:
            break
        try:
            vet.vet(con, tic, FIND)
        except Exception as e:
            log('vet failed', tic, traceback.format_exc()[-400:])
            # archive or network hiccups are retried on the next run; anything else is a bug to look at
            transient = isinstance(e, (requests.exceptions.RequestException, ConnectionError, TimeoutError)) or                 type(e).__name__ in ('RemoteServiceError', 'HTTPError', 'ReadTimeout')
            con.execute("UPDATE candidate SET status=?, updated=? WHERE tic=?",
                        ('new' if transient else 'vet_error', time.strftime('%Y-%m-%d %H:%M'), tic))
            con.commit()


def analyse_binaries(deadline):
    """Stellar signals rejected as planets: cross-match with variable-star catalogues, then analyse every
    uncatalogued one as an eclipsing binary (binaries/eb/<TIC>.json/.png)."""
    from .binaries import crossmatch as new_variables, analysis as eb_vet
    from astropy.table import Table
    try:
        new_variables.main()
    except SystemExit as e:
        log("binary analysis: catalogue cross-match failed", e); return
    eb_vet.OUT.mkdir(parents=True, exist_ok=True)
    for row in Table.read(paths.BIN / 'new_variables_isolated.csv'):
        if time.time() > deadline:
            break
        if not (eb_vet.OUT / f'{int(row["tic"])}.json').exists():
            eb_vet.one(row)


# ---------------------------------------------------------------- entry points
def run(minutes, mode='both'):
    """One budgeted, resumable work session. mode: 'planets' searches and vets planet candidates only; 'binaries' and
    'both' also analyse every stellar signal that no catalogue lists as an eclipsing binary."""
    if LOCK.exists() and time.time() - LOCK.stat().st_mtime < (minutes + 30) * 60:
        try:
            pid = int(LOCK.read_text())
            import psutil  # noqa
            if psutil.pid_exists(pid):
                log('another run is active, exiting'); return
        except Exception:
            log('another run is active, exiting'); return
    LOCK.write_text(str(os.getpid()))
    deadline = time.time() + minutes * 60
    try:
        for d in (RAW, CACHE, FIND, LOGS):
            d.mkdir(parents=True, exist_ok=True)
        for f in RAW.glob('*'):          # leftovers of an interrupted run
            f.unlink()
        if dir_gb(ROOT) > CFG['disk_cap_gb']:
            log('disk cap reached, not downloading'); return
        con = db()
        try:
            refresh_known(con)
        except requests.exceptions.RequestException as e:
            if not con.execute('SELECT 1 FROM known LIMIT 1').fetchone():
                log('the lists of known planets and binaries could not be downloaded, so that no signal can be judged '
                    'new; the archive may be busy, so try again later:', repr(e)[:200]); return
            log('the lists of known objects could not be refreshed; the previous lists are used:', repr(e)[:200])
        vet_new(con, deadline - 600)      # finish any outstanding vetting first
        while time.time() < deadline - 120:
            nj = next_job(con)
            if nj is None:
                log('the whole plan has been searched'); break
            jid, m, tier, u = nj
            if con.execute('SELECT status FROM job WHERE job=?', (jid,)).fetchone() is None:
                if not prepare_job(con, jid, m, tier, u, deadline - 300):
                    log('catalogue still building'); break
                continue
            work_job(con, jid, m, u, deadline - 240)
            promote(con)
        promote(con)
        vet_new(con, deadline - 300)
        if mode in ('binaries', 'both'):
            try:
                analyse_binaries(deadline - 60)
            except Exception:
                log("binary analysis failed", traceback.format_exc()[-300:])
        status(con)
        try:
            from . import dashboard as dashboard; dashboard.write()
        except Exception:
            log('dashboard failed', traceback.format_exc()[-200:])
    finally:
        LOCK.unlink(missing_ok=True)


def status(con=None):
    con = con or db()
    q = lambda sql, *a: con.execute(sql, a).fetchall()
    bands = q('SELECT COUNT(*), SUM(n) FROM tic_band')[0]
    log(f'catalogue bands {bands[0]}, stars {bands[1] or 0}')
    for s, st, n, fin in q('SELECT job, status, n_targets, finished FROM job ORDER BY started'):
        u = int(s.split(':')[-1])
        d = q('SELECT COUNT(*) FROM target WHERE sector=? AND status>0', u)[0][0]
        log(f'  sector {s}: {st} {d}/{n} {fin or ""}')
    log('events', q('SELECT COUNT(*) FROM event')[0][0], '| bls hits', q('SELECT COUNT(*) FROM bls')[0][0])
    for r in q('SELECT status, COUNT(*) FROM candidate GROUP BY status'):
        log('  candidates', r)
    log(f'disk used {dir_gb(ROOT):.2f} GB of {CFG["disk_cap_gb"]} GB')


def vet_one(tic):
    """Vet one star by hand, adding it as a candidate when the survey has not promoted it."""
    from . import planets as vet
    con = db()
    if not con.execute('SELECT 1 FROM candidate WHERE tic=?', (tic,)).fetchone():
        con.execute("INSERT INTO candidate(tic, kind, status, created) VALUES (?,?,?,?)",
                    (tic, 'manual', 'new', time.strftime('%Y-%m-%d %H:%M'))); con.commit()
    vet.vet(con, tic, FIND)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--minutes', type=float, default=55)
    sp.add_parser('status')
    v = sp.add_parser('vet'); v.add_argument('tic', type=int)
    a = ap.parse_args()
    if a.cmd == 'run':
        run(a.minutes)
    elif a.cmd == 'status':
        status()
        from . import dashboard as dashboard; print(dashboard.write())
    elif a.cmd == 'vet':
        from . import planets as vet
        con = db()
        if not con.execute('SELECT 1 FROM candidate WHERE tic=?', (a.tic,)).fetchone():
            con.execute("INSERT INTO candidate(tic, kind, status, created) VALUES (?,?,?,?)",
                        (a.tic, 'manual', 'new', time.strftime('%Y-%m-%d %H:%M'))); con.commit()
        vet.vet(con, a.tic, FIND)
