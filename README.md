<p align="center"><img src="docs/logo.png" width="180" alt="occulta logo"></p>

# occulta

occulta, Latin for hidden things, is a blind search of the photometry of the Transiting Exoplanet Survey Satellite (TESS), K2, and Kepler for transiting planets and eclipsing binaries that no published catalogue contains. It streams the light curves of every star in a chosen brightness range, searches each of them for individual dips and periodic signals, stitches dips from different sectors into candidates, and vets every candidate first as a planet and then, when the signal is stellar, as an eclipsing binary. The user chooses which of the two hunts to run: planets alone, eclipsing binaries alone, or both together.

> **occulta is public** (pre-release 0.0.1). It is ready for use and for testing; results can still change between versions. Please report problems, wrong results and suggestions on the [issue tracker](https://github.com/jovian-explorer/occulta/issues/new/choose).

## Why a blind search

Space photometry records transiting planets and eclipsing binaries whether or not anyone looks for them, so that the archives of TESS, K2, and Kepler still contain signals that the catalogues built from the same data have not listed, in particular single transits from different years, long periods, and eclipses diluted by neighbouring stars. At the same time, a signal is new only when it is absent from every catalogue, when its period is the true one, and when it belongs to the star it is attributed to, and each of these conditions has failed for some of the candidates that this search has examined. In occulta, therefore, every step that decides whether a detection is real is part of the code, documented in METHODS.md, and repeatable from public data with a single command.

## Installation

```
pip install git+https://github.com/jovian-explorer/occulta
```

or, from a copy of this repository,

```
pip install .
```

Python 3.10 or later is required. All data are downloaded from the Mikulski Archive for Space Telescopes (MAST), the Gaia archive, and VizieR as they are needed.

## Usage

Every command works inside one working directory, which holds the database of the search, the cache of light curves, and the findings; it is chosen with `--home`, with the environment variable `OCCULTA_HOME`, or, by default, as the current directory.

```
occulta init                                      # create the working directory and its configuration
occulta survey run --minutes 55 --mode both       # one resumable work session of the blind search
occulta survey status                             # progress and candidates
occulta survey map                                # map of the plan and of the stars searched
occulta reproduce 181149914 --mode both           # repeat the whole analysis of one star from public data
occulta planets vet 402710280                     # vet one star as a planet host
occulta binaries --help                           # the tools for eclipsing binaries
```

The option `--mode` decides what is hunted. With `planets`, the survey searches the light curves and vets each candidate as a planet, and stellar signals are set aside. With `binaries`, every signal that the planet vetting rejects for a stellar reason, and that no catalogue of variable stars contains, is analysed as an eclipsing binary; with `both`, the default, the two hunts run together.

### Tools for eclipsing binaries

| Command | What it does |
|---|---|
| `occulta binaries crossmatch` | cross-matches every stellar signal with VSX, GCVS, Gaia DR3, the TESS eclipsing-binary catalogues, and SIMBAD |
| `occulta binaries analyse` | measures the ephemeris, both eclipses, the secondary phase, odd and even eclipses, and centroid shifts |
| `occulta binaries triage TIC ...` | one line per binary, with the Roche-lobe test that rejects giant hosts |
| `occulta binaries eclipses TIC P T0 PHASE2` | fits every eclipse individually and measures the secondary phase sector by sector |
| `occulta binaries timing TIC P T0 OUT.png` | eclipse timing diagram, with QLP and SPOC times compared in every shared sector |
| `occulta binaries pixels TIC P T0 DUR RA DEC` | assigns the eclipses to a star with difference images of the TESS pixels |
| `occulta binaries photometry TIC RA DEC` | uniform aperture photometry of the full-frame images, for depths that are comparable between sectors |
| `occulta binaries gaia` | Gaia DR3 orbits, accelerations, and spectroscopic solutions of the reviewed binaries |
| `occulta binaries lighttime ...` | fits the eclipse times with the light-travel time of a Gaia astrometric orbit |
| `occulta binaries report TIC` | parameter table and figure for a report |

## The plan of the search

The configuration in `config.json` lists the data in the order in which they are searched, and each step runs from the newest to the oldest sector.

| Step | Data | Stars |
|---|---|---|
| 1 | TESS full-frame images, MIT Quick-Look Pipeline (QLP) | TESS magnitude ≤ 10 |
| 2 | K2, K2SFF light curves | every K2 target |
| 3 | TESS 2-minute cadence, SPOC | every 2-minute target |
| 4 | TESS QLP | 10 < T ≤ 11 |
| 5 | Kepler long cadence | every Kepler target |
| 6 | TESS QLP | 11 < T ≤ 12 |
| 7 | TESS QLP | 12 < T ≤ 13.5 |

## Tests

```
pip install .[test]
pytest                          # offline checks
OCCULTA_ONLINE=1 pytest         # also recovers two known objects from public data (about 15 minutes)
```

The recovery tests require that the pipeline find, with no prior information, the planet candidate TIC 402710280.01 at 36.340 d and the eccentric eclipsing binary CPD−37 4813 at 10.7598 d.

## Methods and citation

METHODS.md describes every step of the search and of the vetting, with references in docs/refs.bib. If you use occulta, please cite it as described in [CITATION.cff](CITATION.cff): Aggarwal, K. 2026, occulta: a blind search of TESS, K2 and Kepler photometry for uncatalogued transiting planets and eclipsing binaries.

## Licence

MIT; see LICENSE.
