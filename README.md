<p align="center"><img src="docs/logo.png" width="180" alt="occulta logo"></p>

# occulta

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23634486.svg)](https://doi.org/10.5281/zenodo.23634486)

> **occulta is public** (pre-release 0.0.1). It is ready for use and for testing; results can still change between versions. Please report problems, wrong results and suggestions on the [issue tracker](https://github.com/jovian-explorer/occulta/issues/new/choose). To cite occulta, use the DOI [10.5281/zenodo.23634486](https://doi.org/10.5281/zenodo.23634486) (all versions; see [CITATION.cff](CITATION.cff)).

occulta, Latin for hidden things, is a blind search of the photometry of the Transiting Exoplanet Survey Satellite (TESS), K2, and Kepler for transiting planets and eclipsing binaries that no published catalogue contains. Space photometry records such objects whether or not anyone looks for them, so that the archives still contain signals that the catalogues built from the same data have not listed, in particular single transits from different years, long periods, and eclipses diluted by neighbouring stars. Since a signal is new only when it is absent from every catalogue, when its period is the true one, and when it belongs to the star it is attributed to, every step that decides whether a detection is real is part of the code, documented in [METHODS.md](METHODS.md), and repeatable from public data with a single command.

## What can I do with it?

* **Hunt for planets, eclipsing binaries, or both**: occulta streams the light curves of every star in a chosen brightness range, searches each of them for individual dips and periodic signals, stitches dips from different sectors into candidates, and vets every candidate first as a planet and then, when the signal is stellar, as an eclipsing binary.
* **Check whether a signal is already known**: every candidate is cross-matched with VSX, GCVS, Gaia DR3, the TESS eclipsing-binary catalogues, the TESS Objects of Interest, and SIMBAD.
* **Find the star that is eclipsed**: difference images of the TESS pixels, modelled with the Gaia positions of every nearby star, assign the eclipses to their source.
* **Measure an eclipsing binary**: the ephemeris, both eclipses, the secondary phase sector by sector, odd and even eclipses, eclipse timing, and the light-travel time of a third body compared with Gaia astrometric orbits.
* **Repeat any result**: `occulta reproduce TIC` repeats the whole analysis of one star from public data, so that every claim can be checked by others.

## Quick start

1. Install occulta with pip (Python 3.10 or newer; details in [How to Install and Run](#-how-to-install-and-run)):
   ```
   pip install "git+https://github.com/jovian-explorer/occulta.git"
   ```
2. Create an empty working directory, enter it, and run one session of the search (the same commands work in the Windows Command Prompt and PowerShell, in macOS, and in Linux; a command given outside a working directory stops with a message instead of filling the current folder):
   ```
   mkdir occulta_work
   cd occulta_work
   occulta init
   occulta survey run --minutes 55 --mode both
   ```
3. Follow the progress with `occulta survey status` and `occulta survey map`; candidates and reports go to `findings/`.

### Worked example: WASP-18 b

To see the code recover a known planet before any search is started, repeat the analysis of WASP-18 (TIC 100100827), a hot Jupiter that Hellier et al. ([2009](https://doi.org/10.1038/nature08245)) found with an orbital period of 0.94 d. In the same working directory, run

```
occulta reproduce 100100827 --mode planets
```

occulta downloads the 23 light curves of the star from MAST, searches each of them with no information on the planet, stitches the dips into one candidate, vets it, and looks the star up in the catalogues. On an ordinary connection the run takes about 10 minutes, most of it the download, and the result is written to `reproduce_out/100100827/`. The expected values, from a run on 26 October 2026 with 226 transits in 10 TESS sectors, are:

| Quantity | occulta | Literature |
|---|---|---|
| Period | 0.94145 d | 0.94145 d (Hellier et al. 2009) |
| Depth | 9,670 ppm (0.97 %) | about 1 % |
| Duration | 2.20 h | about 2.2 h |

These numbers move in the last digit as TESS adds sectors. The report `reproduce_out/100100827/TIC100100827/report.md` marks the candidate as rejected, which is the correct outcome for the code: the occultation of this very hot planet is a real secondary eclipse (24 sigma here), which the planet vetting treats as a stellar signal, and the catalogue step lists the signal as already known (VSX, Gaia DR3, the TESS eclipsing-binary catalogue, and the TOI list as TOI 185.01), so that it would not be counted as new.

---

## 📡 Data Sources & Provenance

All data are downloaded from the Mikulski Archive for Space Telescopes (MAST), the Gaia archive, and VizieR as they are needed, and the search runs through them in the order set in `config.json`, each step from the newest to the oldest sector.

| Step | Data | Stars |
|---|---|---|
| 1 | TESS full-frame images, MIT Quick-Look Pipeline (QLP) | TESS magnitude ≤ 10 |
| 2 | K2, K2SFF light curves | every K2 target |
| 3 | TESS 2-minute cadence, Science Processing Operations Center (SPOC) | every 2-minute target |
| 4 | TESS QLP | 10 < T ≤ 11 |
| 5 | Kepler long cadence | every Kepler target |
| 6 | TESS QLP | 11 < T ≤ 12 |
| 7 | TESS QLP | 12 < T ≤ 13.5 |

Further sources are opt-in and are not in the default plan. TESS-SPOC full-frame-image light curves (Caldwell et al. 2020; up to about 160,000 stars per sector, many more than the 2-minute SPOC set) are searched when `{"mission": "TESS-SPOC"}` is added to `plan` in `config.json`; they are stored as units 4001 and up and read with the same PDCSAP flux and quality flag as the 2-minute SPOC data.

Light curves of the Zwicky Transient Facility (ZTF; Masci et al. 2019), which are sampled about once a night over several years, are searched for one position at a time with `occulta ztf RA DEC`; they carry the unit codes 5001 (zg), 5002 (zr) and 5003 (zi). Only deep events on bright stars can be found in such sparse data, periods near one day are flagged as aliases, and the vetting of ZTF candidates is not yet part of occulta.

JWST and Ariel are not sources: JWST observes chosen targets in short pointed visits, so there is no survey time series to search blindly, and Ariel launches in 2029 and has no data yet. PLATO and CHEOPS are possible later sources.

The catalogues used for vetting are VSX, GCVS, Gaia DR3 (including its non-single-star orbits), the TESS eclipsing-binary catalogues of Prša et al. (2022) and Kostov et al. (2025), the TESS Input Catalog, and SIMBAD; METHODS.md gives the references, collected in docs/refs.bib.

---

## 🚀 How to Install and Run

### Install with pip (Python 3.10+)

occulta is not on PyPI; pip installs it straight from GitHub:

```
pip install "git+https://github.com/jovian-explorer/occulta.git"
```

To work on the code, or to run the tests, install it from a copy of the repository (or from the source archive on the [Releases](https://github.com/jovian-explorer/occulta/releases) page):

```
git clone https://github.com/jovian-explorer/occulta.git
cd occulta
pip install ".[test]"
```

### Command-line options

```
occulta init                                      # create the working directory and its configuration
occulta survey run --minutes 55 --mode both       # one resumable work session of the blind search
occulta survey status                             # progress and candidates
occulta survey map                                # map of the plan and of the stars searched, on the TESS sector count if tess-point is installed (pip install tess-point)
occulta survey map --refresh-totals               # the same, with the MAST light-curve totals counted now (cached in the home)
occulta reproduce 181149914 --mode both           # repeat the whole analysis of one star from public data
occulta planets vet 402710280                     # vet one star as a planet host
occulta binaries --help                           # the tools for eclipsing binaries
occulta rogues search TIC                         # brightenings with the shape of a microlensing event (free-floating planets)
occulta moons search                              # transit timing, duration and orbital-sampling-effect tests for moons of vetted planets
occulta singles scan                              # single transits and aperiodic dippers among the dips already measured
```

The option `--mode` decides what is hunted. With `planets`, the survey searches the light curves and vets each candidate as a planet, and stellar signals are set aside. With `binaries`, every signal that the planet vetting rejects for a stellar reason, and that no catalogue of variable stars contains, is analysed as an eclipsing binary; with `both`, the default, the two hunts run together.

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

### Where occulta keeps its files

Every command works inside one working directory, which holds the database of the search and the cache of light curves (`data/`), the logs (`logs/`), the binary analyses (`binaries/`), and the findings (`findings/`). It is chosen with `--home`, with the environment variable `OCCULTA_HOME`, or, by default, as the current directory when that folder was created with `occulta init`.

---

## Licence

occulta is free software under the MIT licence (`LICENSE`).

---

## 👨‍🔬 Creator & Lead Architect

* **Keshav Aggarwal**
  * **Current Affiliation:** Student Visitor, Space Physics Laboratory (SPL), Vikram Sarabhai Space Centre (VSSC), Indian Space Research Organisation (ISRO), Thiruvananthapuram.
  * **Doctoral Education:** Ph.D. in Astronomy and Space Engineering, Department of Astronomy, Astrophysics and Space Engineering (DAASE), Indian Institute of Technology (IIT) Indore (Thesis Submitted).
  * **Fellowship:** Prime Minister's Research Fellow (PMRF).
  * **Research Focus:** Planetary/exoplanetary atmospheres, ionospheres, planetary radio occultation, solar-planetary space weather interactions across MOM, Chandrayaan-2, Chandrayaan-3, Akatsuki, and Aditya-L1.

---

## 🧪 Automated Test Suite

```
pip install ".[test]"
python -m pytest                          # offline checks
OCCULTA_ONLINE=1 python -m pytest         # also recovers two known objects from public data (about 15 minutes); on Windows: set OCCULTA_ONLINE=1 first
```

The recovery tests require that the pipeline find, with no prior information, the planet candidate TIC 402710280.01 at 36.340 d and the eccentric eclipsing binary CPD−37 4813 at 10.7598 d. The workflow `.github/workflows/ci.yml` runs the offline checks on Windows, macOS and Linux for every push to `main` and every pull request, on a copy of the tests against the installed package.
