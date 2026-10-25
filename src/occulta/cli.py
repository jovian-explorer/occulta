"""Command line of occulta. Every command works inside one working directory (--home, the environment variable
OCCULTA_HOME, or the current directory), which holds the database, the cache, and the findings.

    occulta init                                   create the working directory and its configuration
    occulta survey run --minutes 55 --mode both    one resumable work session of the blind search
    occulta survey status | map                    progress of the search, as text or as a map
    occulta reproduce TIC [--mode planets|binaries|both]
                                                   repeat the whole analysis of one star from public data
    occulta planets vet TIC                        vet one star as a planet host
    occulta binaries <tool> ...                    analyses of eclipsing binaries (occulta binaries --help)
"""
import argparse
import os
import sys

BINARY_TOOLS = {
    'crossmatch': ('crossmatch', 'cross-match every stellar signal with the catalogues of variable stars'),
    'analyse': ('analysis', 'analyse uncatalogued stellar signals as eclipsing binaries ([--one TIC] [--workers N])'),
    'triage': ('triage', 'one-line summary of analysed binaries: TIC ...'),
    'eclipses': ('eclipses', 'fit every eclipse individually: TIC P T0 PHASE2 [WIDTH_H]'),
    'timing': ('timing', 'eclipse timing diagram, QLP against SPOC: TIC P T0 OUT.png'),
    'pixels': ('pixels', 'which star in the TESS pixels eclipses: TIC P T0 DURATION_H RA DEC [PHASE]'),
    'pixels-all': ('pixels_all', 'the pixel test for many binaries: [TIC ...]'),
    'photometry': ('photometry', 'uniform aperture photometry of the full-frame images: TIC RA DEC [N_AP]'),
    'gaia': ('gaia', 'Gaia DR3 non-single-star solutions of the reviewed binaries'),
    'lighttime': ('lighttime', 'fit eclipse times with the light-travel time of a Gaia orbit: TIC P T0 POUT E TPERI A B F G [PNG]'),
    'lighttime-scan': ('lighttime_scan', 'scan outer periods of a light-travel-time term: TIC P T0 [POUT ...]'),
    'dips': ('dips', 'list every dip of a star, sector by sector: TIC [K]'),
    'report': ('report', 'parameter table and figure for a report: TIC ...'),
    'catalogue': ('catalogue', 'rebuild the catalogue of reviewed binaries'),
    'rerun': ('rerun', 're-analyse binaries at their adopted periods: TIC ...'),
}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    ap = argparse.ArgumentParser(prog='occulta', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--home', help='working directory (default: OCCULTA_HOME or the current directory)')
    ap.add_argument('--version', action='store_true', help='print the version and exit')
    sp = ap.add_subparsers(dest='cmd')
    sp.add_parser('init', help='create the working directory and its configuration')
    s = sp.add_parser('survey', help='the blind search')
    ss = s.add_subparsers(dest='sub', required=True)
    r = ss.add_parser('run', help='one budgeted, resumable work session')
    r.add_argument('--minutes', type=float, default=55)
    r.add_argument('--mode', choices=('planets', 'binaries', 'both'), default='both',
                   help='planets: search and vet planet candidates only; binaries or both: also analyse every '
                        'uncatalogued stellar signal as an eclipsing binary')
    ss.add_parser('status', help='progress and candidates')
    ss.add_parser('map', help='map of the plan and of the stars searched (findings/progress/progress_map.png)')
    rp = sp.add_parser('reproduce', help='repeat the analysis of one star from public data', add_help=False)
    rp.add_argument('args', nargs=argparse.REMAINDER)
    pl = sp.add_parser('planets', help='planet hunting')
    pls = pl.add_subparsers(dest='sub', required=True)
    v = pls.add_parser('vet', help='vet one star as a planet host'); v.add_argument('tic', type=int)
    b = sp.add_parser('binaries', help='eclipsing-binary hunting',
                      description='\n'.join(f'  {k:15s} {d}' for k, (_, d) in BINARY_TOOLS.items()),
                      formatter_class=argparse.RawDescriptionHelpFormatter)
    b.add_argument('tool', choices=list(BINARY_TOOLS))
    b.add_argument('args', nargs=argparse.REMAINDER)
    a = ap.parse_args(argv)

    if a.version:
        from . import __version__
        print(__version__); return
    if a.home:
        os.environ['OCCULTA_HOME'] = os.path.abspath(a.home)
    if a.cmd is None:
        ap.print_help(); return
    # the working directory is fixed only now, so that every module sees the same one
    from . import paths
    if a.cmd == 'init':
        print('working directory:', paths.init())
    elif a.cmd == 'survey':
        paths.init()
        from . import survey
        if a.sub == 'run':
            survey.run(a.minutes, a.mode)
        elif a.sub == 'status':
            survey.status()
            from . import dashboard
            print(dashboard.write())
        else:
            from . import progress
            progress.main()
    elif a.cmd == 'reproduce':
        from . import reproduce
        reproduce.main(a.args)
    elif a.cmd == 'planets':
        paths.init()
        from . import survey
        survey.vet_one(a.tic)
    elif a.cmd == 'binaries':
        import importlib
        mod = importlib.import_module(f'occulta.binaries.{BINARY_TOOLS[a.tool][0]}')
        mod.main(a.args)


if __name__ == '__main__':
    main()
