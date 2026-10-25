"""Rebuild findings/objects/CATALOGUE_B (table.md rows, catalogue.csv, catalogue.json) from the current eb_vet results.
Ephemerides (P, T0) come from eb/<TIC>.json; the secondary phase is re-measured per sector from individually fitted
eclipses (eclipses.py, relative to same-sector primaries) and averaged over sectors with inverse-variance
weights; depths, names and notes are kept from the existing catalogue.json. Stars listed in MOVED have their own
reports and are dropped. Usage: occulta binaries catalogue"""
import json, subprocess, sys
from pathlib import Path

from .. import paths

import numpy as np

B = paths.BIN
ROOT = paths.HOME
CAT = ROOT / 'findings' / 'objects' / 'CATALOGUE_B'
# Decisions made by hand while the catalogue is reviewed live in binaries/catalogue_overrides.json:
#   moved     stars that received their own reports and leave the catalogue;
#   reassign  eclipses that the difference images (pixels.py) place on another star, so that the row is relabelled
#             to the true source, whose eclipse depth is the depth that star alone must have;
#   dropped   eclipses that belong to a known binary or to a star fainter than the scope of the catalogue;
#   notes     remarks added to the row of a star.
_OV = json.loads((B / 'catalogue_overrides.json').read_text()) if (B / 'catalogue_overrides.json').exists() else {}
MOVED = set(_OV.get('moved', []))
REASSIGN = {int(k): v for k, v in _OV.get('reassign', {}).items()}
DROPPED = {int(k): v for k, v in _OV.get('dropped', {}).items()}
NOTES = {int(k): v for k, v in _OV.get('notes', {}).items()}


def sector_phases(tic, P, T0, ph2, w_h):
    """Weighted mean secondary phase over sectors, from eclipses.py; None when no secondary is fitted."""
    out = subprocess.run([sys.executable, '-W', 'ignore', '-m', 'occulta.binaries.eclipses', str(tic), repr(P), repr(T0),
                          repr(ph2), repr(w_h)], cwd=ROOT, capture_output=True, text=True).stdout
    ph, er = [], []
    for line in out.splitlines()[1:]:
        tok = line.split()
        if len(tok) < 7 or 'nan' in tok[-1]:
            continue
        v, e = tok[-1].split('±')
        ph.append(float(v)); er.append(max(float(e), 2e-4))   # single-eclipse sectors have zero scatter: floor
    if not ph:
        return None, None, 0
    ph, er = np.array(ph), np.array(er)
    W = 1 / er ** 2
    m = float(np.sum(W * ph) / W.sum())
    # error: formal, inflated by the scatter between sectors
    chi = np.sum(W * (ph - m) ** 2) / max(1, ph.size - 1)
    return m, float(np.sqrt(max(1.0, chi) / W.sum())), int(ph.size)


def main(argv=None):
    # key by the TESS target whose light curve was analysed (differs from 'tic' for reassigned rows)
    old = {int(r.get('target_tic') or r['tic']): r for r in json.loads((CAT / 'catalogue.json').read_text(encoding='utf-8'))}
    rows = []
    for tic, r in old.items():
        if tic in MOVED or tic in DROPPED:
            continue
        if r.get('fixed_ephemeris'):
            # ephemeris measured by hand (e.g. when data gaps defeat the automatic analysis); kept as listed
            rows.append(dict(r)); continue
        e = json.loads((B / 'eb' / f'{tic}.json').read_text())
        eph = e.get('ephemeris') or {}
        P, T0 = e['period'], e['t0_btjd']
        dT = ((T0 + (r['ph2'] * P if r.get('swap') and r['ph2'] else 0) - r['T0'] + P / 2) % P) - P / 2
        ph2 = (e.get('secondary') or {}).get('phase') or (r['ph2'] - dT / P if r['ph2'] else 0.5)
        m, s, n = sector_phases(tic, P, T0, ph2, max(1.0, min(e.get('width_h') or 3, P * 24 / 6)))
        if r['ph2'] is None:                                           # no secondary detected before: keep it so
            m = s = None
        if not eph:
            # no timing ephemeris (too few eclipses outside Sectors 99-104): epoch from the individually fitted eclipses
            import pandas as pd
            c = pd.read_csv(B / 'eb' / f'{tic}_eclipses.csv')
            c = c[(c.type == 'P') & ~c.unit.between(99, 104)]
            if len(c):
                T0 = float(np.median(c.t_btjd - c.cycle * P))
        if r.get('swap') and m is not None:                            # the deeper eclipse is listed as primary
            # its epoch is taken directly from the individually fitted deeper ('S') eclipses, since the phase of a
            # single eclipse per sector is too uncertain to transfer an epoch
            import pandas as pd
            c = pd.read_csv(B / 'eb' / f'{tic}_eclipses.csv')
            c = c[(c.type == 'S') & ~c.unit.between(99, 104)]
            guess = T0 + m * P
            T0d = float(np.median(c.t_btjd - np.round((c.t_btjd - guess) / P) * P)) if len(c) else guess
            m = 1 - (((T0d - T0) / P) % 1)
            T0 = T0d
        if tic in NOTES and NOTES[tic] not in (r.get('sec_note') or ''):
            r = dict(r, sec_note=('; '.join(x for x in (r.get('sec_note'), NOTES[tic]) if x)))
        if tic in REASSIGN:
            r = dict(r, **REASSIGN[tic], target_tic=tic)
        rows.append(dict(r, P=P, Perr=eph.get('P_err_d'), T0=T0, ph2=m, ph2_err=s, ph2_nsec=n,
                         ecw=None if m is None else float(np.pi / 2 * (m - 0.5))))
        print(tic, r['name'], 'P %.7f' % P, 'T0 %.5f (moved %.1f min)' % (T0, dT * 1440),
              'ph2 old %s new %s ± %s (%d sectors)' % (r['ph2'] and round(r['ph2'], 4), m and round(m, 4), s and round(s, 4), n), flush=True)
    rows.sort(key=lambda x: x['P'])
    (CAT / 'catalogue.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')
    with open(CAT / 'catalogue.csv', 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('tic,name,ra,dec,tmag,teff,period_d,period_err_d,t0_btjd,primary_depth_pct,secondary_depth_pct,'
                 'secondary_phase,secondary_phase_err,note\n')
        for x in rows:
            fh.write(','.join('' if v is None else str(v) for v in (
                x['tic'], x['name'], x['ra'], x['dec'], x['tmag'], x['teff'], x['P'], x['Perr'], x['T0'], x['D1'],
                x['D2'], x['ph2'], x['ph2_err'], x['sec_note'])) + '\n')
    lines = ['| # | TIC | Name | RA, Dec (deg) | T | Teff (K) | P (d) | T0 (BTJD) | Primary depth (%) | '
             'Secondary depth (%) | Secondary phase | Note |', '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for i, x in enumerate(rows, 1):
        teff = '-' if not x['teff'] or x['teff'] != x['teff'] else '%d' % x['teff']
        per = '%.6f ± %.6f' % (x['P'], x['Perr']) if x['Perr'] else '%.6f' % x['P']
        ph = '-' if x['ph2'] is None else '%.4f ± %.4f' % (x['ph2'], x['ph2_err'])
        d2 = '-' if x['D2'] is None else '%.2f' % x['D2']
        lines.append(f"| {i} | {x['tic']} | {x['name']} | {x['ra']:.4f}, {x['dec']:.4f} | {x['tmag']:.2f} | {teff} | "
                     f"{per} | {x['T0']:.4f} | {x['D1']:.2f} | {d2} | {ph} | {x['sec_note']} |")
    (CAT / 'table.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
