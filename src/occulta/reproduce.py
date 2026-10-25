"""Reproduce one detection from public data alone, blind, with the same code the survey uses.

    occulta reproduce 181149914                    # a TIC number
    occulta reproduce 181149914 --out my_run       # output folder (default reproduce_out/<TIC>)
    occulta reproduce 181149914 --mode planets     # planet vetting only

Steps, each identical to the survey: (1) every light curve of the star is downloaded from MAST; (2) each one
is searched for individual dips and periodic signals (lightcurves.process); (3) dips from different sectors are
stitched into a candidate (survey.promote) with no catalogue of known objects loaded, so nothing is set aside;
(4) the candidate is vetted as a planet (planets.vet); (5) if the signal is stellar, it is analysed as an
eclipsing binary (binaries.analysis.analyse); (6) the star is cross-matched against VSX, GCVS, the Gaia DR3
eclipsing binaries, the TESS EB catalogue, SIMBAD and the TOI/CTOI lists. Everything is written to one
folder together with the list of data files used and checksums of the code that produced it.
"""
import argparse, hashlib, json, platform, shutil, sys, time
from pathlib import Path

from . import paths

ROOT = paths.HOME
import numpy as np  # noqa: E402
from . import survey as hunt, lightcurves as lcsearch, sources, planets as vet  # noqa: E402

PKG = Path(__file__).resolve().parent


def code_hashes():
    """SHA-256 checksums of every module of the installed package, and the package version, so that a reproduction
    records exactly which code produced it."""
    from . import __version__
    h = {f.relative_to(PKG).as_posix(): hashlib.sha256(f.read_bytes()).hexdigest()[:16]
         for f in sorted(PKG.rglob('*.py'))}
    h['version'] = __version__
    return h


def main(argv=None):
    ap = argparse.ArgumentParser(prog='occulta reproduce', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tic', type=int)
    ap.add_argument('--out', default=None)
    ap.add_argument('--mode', choices=('planets', 'binaries', 'both'), default='both',
                    help='planets: stop after the planet vetting; binaries: analyse the signal as an eclipsing binary '
                         'whatever the vetting decides; both (default): analyse it as a binary when it is stellar')
    a = ap.parse_args(argv)
    paths.init()
    tic = a.tic
    out = Path(a.out or ROOT / 'reproduce_out' / str(tic)).resolve()
    if out.exists():
        shutil.rmtree(out)
    (out / 'raw').mkdir(parents=True)
    log = []

    def say(*x):
        m = time.strftime('%H:%M:%S ') + ' '.join(str(v) for v in x)
        print(m, flush=True); log.append(m)

    # 1. data
    files = vet.fetch_all(tic, out / 'raw')
    say(f'{len(files)} light curves downloaded')
    if not files:
        sys.exit('no light curves found at MAST for this star')

    # 2. search every light curve
    con = hunt.db(out / 'run.db')
    for author, unit, path in files:
        u = unit if author in ('QLP',) else (1000 + unit if author in ('SPOC', 'TESS-SPOC') else unit)
        con.execute('INSERT OR IGNORE INTO target(sector, tic, status) VALUES (?,?,0)', (u, tic))
        res = lcsearch.process(path, tic, u, vet.AUTHORS[author], do_bls=author != 'Kepler')
        hunt.store(con, res)
        n = len(res.get('events', []))
        say(f'{author} unit {unit}: {n} dips', '' if 'error' not in res else res['error'])
    con.commit()

    # 3. candidate from the stitched dips (no known-object catalogue loaded)
    hunt.promote(con)
    row = con.execute('SELECT kind, n_events, sectors FROM candidate WHERE tic=?', (tic,)).fetchone()
    if row is None:
        say('no repeating signal passed the search thresholds; nothing to vet')
        con.execute("INSERT INTO candidate(tic, kind, status, created) VALUES (?,?,?,?)",
                    (tic, 'manual', 'new', time.strftime('%Y-%m-%d %H:%M')))
        con.commit()
    else:
        say(f'candidate from the search: {row}')

    # 4. planet vetting
    vet.vet(con, tic, out)
    star_dir = out / sources.name_of(tic).replace(' ', '')
    s = json.loads((star_dir / 'summary.json').read_text()) if (star_dir / 'summary.json').exists() else {}
    say('planet vetting:', s.get('status'), '; '.join(s.get('reasons', [])))

    # 5. eclipsing-binary analysis when the signal is stellar
    eb = None
    stellar = s and s.get('status') == 'rejected' and s.get('depth_ppm', 0) > 0
    if s and s.get('period_d') and (a.mode == 'binaries' or (a.mode == 'both' and stellar)):
        from .binaries import analysis as eb_vet
        eb_vet.OUT = star_dir
        eb_vet.WORK = out / 'raw_eb'
        eb = eb_vet.analyse(tic, s['period_d'], s['t0_bjd'] - 2457000.0, s['t14_h'], s.get('gaia', {}))
        say(f"EB analysis: P {eb['period']:.5f} d, primary {eb['primary_depth']*100:.2f}%, "
            f"secondary {eb['secondary_depth']*100:.2f}% at phase "
            f"{(eb.get('secondary') or {}).get('phase', eb.get('secondary_phase_scan', float('nan'))):.4f}")

    # 6. catalogues
    from .binaries import catalogues as catalogues
    st = s.get('star', {}) if s else {}
    known = catalogues.lookup(tic, st.get('ra'), st.get('dec')) if st.get('ra') is not None else {}
    say('catalogues:', json.dumps(known))

    manifest = {'tic': tic, 'run': time.strftime('%Y-%m-%d %H:%M:%S'), 'python': platform.python_version(),
                'code_sha256': code_hashes(), 'data_files': [Path(p).name for _, _, p in files],
                'planet_vetting': s, 'eclipsing_binary': eb, 'catalogues': known}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1, default=float))
    (out / 'log.txt').write_text('\n'.join(log) + '\n')
    shutil.rmtree(out / 'raw', ignore_errors=True)
    con.close()
    say(f'done: {out}')


if __name__ == '__main__':
    main()
