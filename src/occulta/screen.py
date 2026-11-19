"""Final false-detection screen for the ready set.

planets.vet() and classify.py already test radius, depth, V-shape, odd/even, secondary eclipse, Gaia eclipsing-binary
flags, SNR, baseline steps, depth change, sector-synchronised matches and pulsators. This module adds the checks that
the ready set was never put through, using only what the database already holds (the per-sector event table, the
TIC row and the notes written by the earlier vetting), so it runs in seconds and needs no light-curve download:

  FEW      too few transits: two events in one sector below SNR 30, or three events in one sector below SNR 15
  PROV     every event lies in QLP sectors 99-104, whose time stamps still await correction
  GAP      most events touch an orbit gap or a momentum-dump edge (flag 'edge'), or the period is a TESS orbit alias
           (13.7 d, 27.4 d) and at least one event touches an edge
  CENTROID most events carry a centroid shift flag
  BKG      most events coincide with a background excursion (scattered light)
  DEPTH    event depths differ by more than a factor 3 between sectors, so the signal is not one astrophysical depth
  VSHAPE   most events are V-shaped while the object is labelled a planet
  CROSS    three or more ready stars share one period and one epoch (a common systematic, not three coincidences)
  BLEND    a neighbour contributes more flux than the target (TIC contamination ratio >= 1)
  RADIUS   (relabel, not reject) a planet whose radius, corrected for dilution, exceeds 2 R_J
  ODDEVEN, SECOND  (not for eclipsing binaries) repeats of the stored odd/even and secondary statistics against the vetting thresholds

reject(row, events, star, siblings) returns a list of (code, detail); an empty list keeps the candidate."""
import re

ORBIT = 13.70            # TESS orbit, days
PROV = range(99, 105)    # QLP sectors with provisional time stamps


def matched(events, p, t0, dur_h):
    """Events of the star that fall on the ephemeris (within a transit duration or 2 per cent of the period)."""
    tol = max(dur_h / 24.0, 0.02 * p)
    out = []
    for sec, t, depth, flags in events:
        ph = (t - t0) / p
        if abs(ph - round(ph)) * p <= tol:
            out.append((sec, depth, set((flags or '').split(','))))
    return out


def _frac(ev, flag):
    return sum(flag in f for _, _, f in ev) / len(ev) if ev else 0.0


def _num(pat, text):
    m = re.search(pat, text or '')
    return float(m.group(1)) if m else None


def reject(row, events, star, siblings=0):
    """row: dict(tic, kind, period, t0, dur_h, snr, n_events, sectors, objtype, notes). events: [(sector, t0, depth_ppm, flags)].
    star: dict(contratio) or {}. siblings: number of other ready stars on the same ephemeris."""
    out = []
    p, notes, label = row['period'], row.get('notes') or '', row.get('objtype')
    ev = matched(events, p, row['t0'], row['dur_h'])
    secs = {int(x) for x in str(row.get('sectors') or '').split(',') if x.strip().isdigit()} or {s for s, _, _ in ev}
    nev = row['n_events']
    tess = row['tic'] < 10 ** 10
    if len(secs) <= 1 and ((nev == 2 and row['snr'] < 30) or (nev == 3 and row['snr'] < 15)):
        out.append(('FEW', f'{nev} events in {max(len(secs), 1)} sector, SNR {row["snr"]:.1f}'))
    if tess and secs and all(s in PROV for s in secs):
        out.append(('PROV', f'all events in provisional QLP sector(s) {sorted(secs)}'))
    if len(ev) >= 2 and _frac(ev, 'edge') >= 0.5:
        out.append(('GAP', f'{_frac(ev, "edge"):.0%} of {len(ev)} events touch a data gap or momentum-dump edge'))
    elif tess and any(abs(p - k * ORBIT) / (k * ORBIT) < 0.03 for k in (1, 2)) and ev and _frac(ev, 'edge') > 0:
        out.append(('GAP', f'period {p:.2f} d is a TESS orbit alias and an event touches an edge'))
    if len(ev) >= 2 and (_frac(ev, 'centroid_x') + _frac(ev, 'centroid_y') >= 1.0 or
                         sum(('centroid_x' in f or 'centroid_y' in f) for _, _, f in ev) / len(ev) >= 0.5):
        out.append(('CENTROID', f'centroid shift in {sum(("centroid_x" in f or "centroid_y" in f) for _, _, f in ev)} of {len(ev)} events'))
    if len(ev) >= 2 and _frac(ev, 'bkg') >= 0.5:
        out.append(('BKG', f'{_frac(ev, "bkg"):.0%} of {len(ev)} events coincide with a background excursion'))
    bysec = {}
    for s, d, _ in ev:
        bysec.setdefault(s, []).append(d)
    med = [sorted(v)[len(v) // 2] for v in bysec.values() if v]
    if len(med) >= 2 and min(med) > 0 and max(med) / min(med) > 3.0:
        out.append(('DEPTH', f'sector depths differ by {max(med) / min(med):.1f}x ({min(med):.0f}-{max(med):.0f} ppm)'))
    if label == 'planet' and len(ev) >= 2 and _frac(ev, 'vshape') >= 0.5:
        out.append(('VSHAPE', f'{_frac(ev, "vshape"):.0%} of events are V-shaped for a planet label'))
    if siblings >= 2:
        out.append(('CROSS', f'{siblings + 1} ready stars share this period and epoch'))
    cr = (star or {}).get('contratio')
    if cr is not None and cr >= 1.0:
        out.append(('BLEND', f'TIC contamination ratio {cr:.2f}: a neighbour is brighter than the target'))
    oe = _num(r'odd/even ([\d.]+) sigma', notes)
    if label != 'eclipsing binary' and oe is not None and oe > 3.0:       # an EB's odd/even difference only doubles its period
        out.append(('ODDEVEN', f'odd/even depths differ at {oe:.1f} sigma'))
    sec = _num(r'secondary (-?[\d.]+) sigma', notes)
    if label == 'planet' and sec is not None and sec > 3.0:
        out.append(('SECOND', f'secondary eclipse at {sec:.1f} sigma in a planet candidate'))
    best = _num(r'best ([\d.]+) sigma, threshold ([\d.]+)', notes)
    thr = _num(r'threshold ([\d.]+)', notes)
    if label == 'planet' and best is not None and thr and best >= thr:
        out.append(('SECOND', f'best secondary at any phase {best:.1f} sigma, above threshold {thr:.1f}'))
    return out


def relabel(row, star):
    """New object type after the dilution-corrected radius test, or None. A 'planet' larger than 2 R_J once the TIC
    contamination is removed is a companion of unknown nature, not a planet."""
    if row.get('objtype') != 'planet':
        return None
    rp = _num(r'Rp ([\d.]+)(?: \([\d.\-]+\))? R_J', row.get('notes'))
    if rp is None:
        return None
    cr = (star or {}).get('contratio') or 0.0
    corrected = rp * (1.0 + max(cr, 0.0)) ** 0.5
    return ('unclassified', f'Rp {rp:.2f} R_J becomes {corrected:.2f} R_J after dilution, above 2 R_J') if corrected > 2.0 else None


def siblings_of(rows, tol=0.002):
    """{tic: number of other rows with the same period (within tol) and epoch (within one duration)}"""
    rows = sorted(rows, key=lambda r: r['period'])
    out = {r['tic']: 0 for r in rows}
    for i, a in enumerate(rows):
        for b in rows[i + 1:]:
            if b['period'] > a['period'] * (1 + tol):
                break
            dph = abs((a['t0'] - b['t0']) / a['period'])
            if min(dph % 1, 1 - dph % 1) * a['period'] <= max(a['dur_h'], b['dur_h']) / 24.0:
                out[a['tic']] += 1
                out[b['tic']] += 1
    return out


def main(argv=None):
    """occulta screen [--apply]. Without --apply it prints the count per reason and writes nothing. With --apply (copy
    hunt.db first) it takes run.lock, writes the reason codes into the notes of each rejected candidate (status
    'rejected'), marks the candidates found by 'occulta recovered' with status 'recovered', and relabels planets that
    are over-sized after the dilution correction as 'unclassified'."""
    import argparse
    import collections
    import csv
    import os
    import sqlite3
    import time
    from . import paths, recovered
    ap = argparse.ArgumentParser(prog='occulta screen')
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--matches', default=str(paths.HOME / 'xmatch' / 'matches.csv'))
    a = ap.parse_args(argv)
    lock = paths.HOME / 'data' / 'run.lock'
    if a.apply:
        if lock.exists() and time.time() - lock.stat().st_mtime < 3600:
            raise SystemExit(f'{lock} is held by another run; try again later')
        lock.write_text(str(os.getpid()))
    con = sqlite3.connect(paths.HOME / 'data' / 'hunt.db')
    try:
        cols = 'tic kind period t0 dur_h snr n_events sectors objtype notes'.split()
        rows = [dict(zip(cols, r)) for r in con.execute(f"SELECT {','.join(cols)} FROM candidate WHERE status='ready'")]
        with open(a.matches, encoding='utf8') as f:
            m = list(csv.DictReader(f))
        for r in m:
            r['tic'] = int(r['tic'])
            r['sep'] = float(r['sep']) if r['sep'] else None
        found = {t: v for t, v in recovered.decide(m).items() if v['verdict'] == 'recovered'}
        sib = siblings_of([r for r in rows if r['tic'] not in found])
        now = time.strftime('%Y-%m-%d %H:%M')
        n = collections.Counter()
        for r in rows:
            note = ''
            status = 'ready'
            if r['tic'] in found:
                x = found[r['tic']]['match']
                status = 'recovered'
                note = f"recovered: {x['cat']} {x['cat_name']} ({x['ref']}), catalogue period {x['cat_period']}, relation {x['relation']}"
                n['recovered'] += 1
            else:
                ev = con.execute('SELECT sector, t0, depth_ppm, flags FROM event WHERE tic=?', (r['tic'],)).fetchall()
                st = con.execute('SELECT contratio FROM tic WHERE tic=?', (r['tic'],)).fetchone()
                st = {'contratio': st[0]} if st else {}
                why = reject(r, ev, st, sib[r['tic']])
                if why:
                    status = 'rejected'
                    note = 'screen ' + '; '.join(f'{c}: {d}' for c, d in why)
                    n['rejected'] += 1
                    for c, _ in why:
                        n[c] += 1
                else:
                    n['kept'] += 1
                    new = relabel(r, st)
                    if new:
                        n['relabelled'] += 1
                        if a.apply:
                            con.execute('UPDATE candidate SET objtype=? WHERE tic=?', (new[0], r['tic']))
                        note = 'screen RADIUS: ' + new[1]
            if a.apply and note:
                con.execute('UPDATE candidate SET status=?, notes=?, updated=? WHERE tic=?',
                            (status, f"{r['notes'] or ''} | {note}", now, r['tic']))
        if a.apply:
            con.commit()
        print(dict(n), '(applied)' if a.apply else '(dry run)')
    finally:
        con.close()
        if a.apply:
            lock.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
