"""Collect everything known about one new eclipsing binary into findings/objects/<TIC>/:
params.json (ephemeris, eclipse depths and widths, orbital constraints, star, Gaia, catalogues), the diagnostic
figure, and a parameter table in Markdown (table.md) to be quoted by the hand-written report.md.

Orbital constraints: for an eclipsing orbit seen close to edge-on, the secondary phase and the ratio of the
eclipse widths give e cos(w) and e sin(w) (Sterne 1939; Kopal 1959):
  phase_2 - phase_1 = 1/2 + (1/pi) [ arctan(e cos w / sqrt(1 - e^2)) + e cos w sqrt(1 - e^2) / (1 - e^2 sin^2 w) ]
  width_2 / width_1 = (1 + e sin w) / (1 - e sin w)     (widths taken at half depth, FWHM)
which are solved together for e and w.
Usage: occulta binaries report <TIC>"""
import json, shutil, sys
from pathlib import Path

from .. import paths

import numpy as np
from scipy.optimize import brentq

ROOT = paths.HOME


def phase_sep(ecw, esw):
    e2 = ecw ** 2 + esw ** 2
    s2 = esw ** 2 / e2 if e2 > 0 else 0.0
    q = np.sqrt(1 - e2)
    return 0.5 + (np.arctan(ecw / q) + ecw * q / (1 - e2 * s2)) / np.pi


def orbit(phase2, width_ratio):
    esw = (width_ratio - 1) / (width_ratio + 1)
    f = lambda ecw: phase_sep(ecw, esw) - phase2
    lim = np.sqrt(max(0.0, 0.99 ** 2 - esw ** 2))
    try:
        ecw = brentq(f, -lim, lim)
    except ValueError:
        return None
    e = float(np.hypot(ecw, esw))
    return {'e_cos_w': float(ecw), 'e_sin_w': float(esw), 'e': e, 'omega_deg': float(np.degrees(np.arctan2(esw, ecw)) % 360)}


def build(tic):
    from . import catalogues as catalogues
    eb = json.loads((paths.BIN / 'eb' / f'{tic}.json').read_text())
    s = json.loads((ROOT / 'findings' / f'TIC{tic}' / 'summary.json').read_text())
    st = s['star']
    out = ROOT / 'findings' / 'objects' / f'TIC{tic}'
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(paths.BIN / 'eb' / f'{tic}.png', out / 'eclipses.png')
    p = {'tic': tic, 'star': st, 'gaia': eb.get('gaia', {}), 'eb': eb}
    sec = eb.get('secondary')
    if sec and sec.get('fwhm_h'):
        p['orbit'] = orbit(sec['phase'], sec['fwhm_h'] / eb['fwhm_h'])
        # error on e cos w from the phase error alone (first order: e cos w ~ pi/2 (phase - 1/2))
        if p['orbit']:
            p['orbit']['e_cos_w_err'] = float(np.pi / 2 * sec['phase_err'])
    # sum of fractional radii for a circular, edge-on orbit, from the primary width: (R1 + R2)/a ~ pi width / P
    p['r_sum_circular'] = float(np.pi * eb['width_h'] / 24 / eb['period'])
    if st.get('mass') and st.get('rad'):
        a_rsun = ((st['mass'] * 2) * (eb['period'] / 365.25) ** 2) ** (1 / 3) * 215.032   # twin-mass assumption
        p['a_rsun_twin_mass'] = float(a_rsun)
    p['catalogues'] = catalogues.lookup(tic, st['ra'], st['dec'])
    (out / 'params.json').write_text(json.dumps(p, indent=1, default=float), encoding='utf-8')
    e = eb.get('ephemeris') or {}
    rows = [('TIC', tic), ('RA, Dec (deg)', f"{st['ra']:.5f}, {st['dec']:.5f}"), ('TESS magnitude', st['tmag']),
            ('Teff (K), TIC', st['teff']), ('Radius (Rsun), TIC', st['rad']), ('Distance (pc)', st['dist_pc']),
            ('Gaia DR3 RUWE', f"{(p['gaia'] or {}).get('ruwe') or float('nan'):.2f}"),
            ('Period (d)', f"{eb['period']:.7f} ± {e.get('P_err_d', float('nan')):.7f}"),
            ('Primary mid-time T0 (BTJD)', f"{eb['t0_btjd']:.5f} ± {e.get('T0_err_d', float('nan')):.5f}"),
            ('Primary depth (%)', f"{eb['primary_depth']*100:.2f} ± {eb['primary_err']*100:.2f}"),
            ('Primary width (h), 10% depth / FWHM', f"{eb['width_h']:.2f} / {eb['fwhm_h']:.2f}"),
            ('Secondary depth (%)', f"{eb['secondary_depth']*100:.2f} ± {eb['secondary_err']*100:.2f}"),
            ('Secondary phase', f"{sec['phase']:.5f} ± {sec['phase_err']:.5f}" if sec else f"{eb['secondary_phase_scan']:.3f} (scan)"),
            ('Secondary width (h), 10% depth / FWHM', f"{sec['width_h']:.2f} / {sec['fwhm_h']:.2f}" if sec and sec.get('fwhm_h') else '-'),
            ('Primary timing O-C rms (min), reduced chi2', f"{e.get('oc_rms_min', float('nan')):.2f}, {e.get('chi2_red', float('nan')):.1f}"),
            ('Eclipses covered / seen', f"{eb['eclipses_covered']} / {eb['eclipses_seen']}"),
            ('Data units (TESS sectors; 1000+ = 2-min)', ', '.join(str(u) for u in eb['units']))]
    if p.get('orbit'):
        o = p['orbit']
        rows += [('e cos ω', f"{o['e_cos_w']:.5f} ± {o['e_cos_w_err']:.5f}"), ('e sin ω (from FWHM ratio)', f"{o['e_sin_w']:.3f}"),
                 ('e (edge-on)', f"{o['e']:.3f}"), ('ω (deg)', f"{o['omega_deg']:.1f}")]
    (out / 'table.md').write_text('| Quantity | Value |\n|---|---|\n' + '\n'.join(f'| {a} | {b} |' for a, b in rows) + '\n',
                                encoding='utf-8')
    print(out)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    for t in argv:
        build(int(t))


if __name__ == '__main__':
    main()
