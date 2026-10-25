# Model audit: every method, model and threshold of occulta

This file lists every place where occulta decides, with a method or a number, whether a signal is real, what it is, and which star it belongs to. For each item it gives where the code is, the value used, the literature source when there is one, and whether an offline test checks it. Line numbers refer to the commit that adds this file. Every reference with a DOI was checked against the registry record with Argus (Crossref, DataCite or OpenAlex) and is listed in `docs/refs.bib`.

Status words: **sourced** (the method or value follows the cited work), **choice** (an empirical setting of this pipeline, documented in METHODS.md but taken from no publication), **fixed** (wrong before this audit and corrected, with a test), and **untested** (no offline test exercises it; only the online recovery tests in `tests/test_recovery.py` run it end to end).

## Fixed in this audit

### F1. `occulta planets vet TIC` could never vet a star by hand
* **Where**: `src/occulta/planets.py:196` (`vet`).
* **Problem**: the dips of each light curve were kept only when their depth lay within a factor of 2 of the candidate depth. A star vetted by hand (`occulta planets vet`, the README example) has no candidate depth, so that every dip was dropped and the star was always rejected with "fewer than two consistent dips".
* **Fix**: without a candidate depth every dip with a signal-to-noise ratio of at least 7 enters, and the pairing of dips that follows already requires similar depths (0.65 to 1.55).
* **Test**: `tests/test_offline_science.py::test_vet_by_hand_finds_the_period` recovers a 7.13 d transit from two synthetic sectors.

### F2. Centroid shifts of Kepler and K2 data were converted with the TESS pixel
* **Where**: `src/occulta/binaries/analysis.py:33-34, 316`.
* **Problem**: the in-eclipse centroid shift was multiplied by 21 arcsec per pixel for every unit. Since a TIC star in the Kepler field also receives Kepler light curves, a Kepler shift was overstated by a factor of 5.3 (Kepler pixels are 3.98 arcsec; Koch et al. 2010). Its significance in units of the error was not affected.
* **Fix**: 3.98 arcsec for units of 2000 and above (K2 and Kepler), 21 arcsec for TESS (Ricker et al. 2015).
* **Test**: `test_binary_centroid_shift_uses_the_pixel_of_each_mission`, which also checks the period and the secondary phase of a synthetic eclipsing binary.

### F3. `occulta binaries lighttime` and `lighttime-scan` stopped with a NameError
* **Where**: `src/occulta/binaries/lighttime.py`, `lighttime_scan.py`.
* **Problem**: both read `paths.BIN` without importing `paths`, so that the two commands failed before reading any eclipse time.
* **Test**: `test_lighttime_scan_finds_an_injected_outer_orbit` recovers a 400 d light-travel-time term from synthetic eclipse times.

### F4. METHODS.md disagreed with the code in three places
* Section 6, item 8, said that the point-spread-function model uses the Gaia sources "brighter than the target by up to 4 magnitudes"; the code (`pixels.py:104`) uses every source brighter than the target and every source fainter by less than 4 magnitudes. The pixels are 21 arcsec, not 20.
* Section 6, item 9, said that the light-travel-time fit at the Gaia epoch is compared with fits at trial epochs spread over the outer period; the code (`lighttime.py:45-49`) moves the outer period by up to 3 d and the epoch of periastron by up to 25 d around the Gaia values.
* Section 9 said `pip install occulta`, but occulta is not on PyPI; it now installs from GitHub.

## Search of each light curve (`lightcurves.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| QLP flux column | `read_lc`, l. 35 | `SYS_RM_FLUX`, else `SAP_FLUX` (not the high-passed `KSPSAP_FLUX`) | Huang et al. 2020a, b; Kunimoto et al. 2021 | sourced |
| SPOC and Kepler flux | l. 47 | PDCSAP, quality flag 0 | Smith et al. 2012; Stumpe et al. 2012; Jenkins et al. 2016 | sourced |
| K2 flux | l. 41 | K2SFF `FCOR`, thruster firings (`MOVING`) dropped | Vanderburg & Johnson 2014 | sourced |
| Kepler/K2 time to BTJD | l. 24 | BKJD + 2454833 − 2457000 | mission definitions of BKJD and BTJD | sourced |
| Detrending | `detrend`, l. 65-80 | running median of 2 d per segment, segments split at gaps > 0.5 d, upward outliers > 4 MAD removed | choice | untested |
| Box durations | l. 10 | 1 to 13 h, nine widths | choice | untested |
| Single-event detection | `single_events` | min(white, red) S/N ≥ 7, at most 5 dips per light curve | choice | untested |
| False-alarm level | `process`, l. 231 | strongest dip of the inverted light curve, away from edges | choice (inverted-light-curve test as in common transit-search practice; no single source) | untested |
| V-shape flag | l. 166 | depth of outer quarters / centre < 0.35 | choice | untested |
| Background and centroid flags | l. 178 | excursion > 6σ against the flanks | choice | untested |
| Common-mode systematics | `survey.mark_systematics` | ≥ 5 stars on one CCD within ±0.1 d | choice | untested |
| BLS | `run_bls` | Astropy BLS (Kovács et al. 2002; Astropy Collaboration 2022), 20-min bins, 0.5 d to 15 d (30 d for baselines > 40 d), durations 1-6 h, objective S/N | sourced | untested |
| BLS SDE | l. 205 | (peak − median) / standard deviation of the S/N periodogram | Kovács et al. 2002 define the SDE with the mean of the signal residue; the median used here is slightly larger than the mean for a skewed periodogram, so that the SDE is marginally lower (conservative) | choice, noted |
| BLS acceptance | `survey.promote`, l. 389; `default_config.json` | S/N ≥ 10, SDE ≥ 9, ≥ 3 transits, odd/even < 3σ | choice | untested |

## Candidates across units (`survey.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Dip S/N for stitching | `promote`; config `event_snr` | ≥ 9 and ≥ 1.2 × false-alarm level | choice | untested |
| Depth and duration agreement | l. 369 | depth ratio 0.6-1.6, duration ratio 0.5-2 | choice | untested |
| Known-object match | `known_match`, l. 333 | ±0.3 d of a TOI, CTOI, KOI, K2 or confirmed-planet ephemeris; any EB-catalogue entry sets the star aside | Prša et al. 2022; Kirk et al. 2016 | sourced |

## Vetting as a planet (`planets.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Period from aliases | l. 202-272 | every P = Δt/n of a pair of dips; rejected when a covered epoch shows < 2σ and falls > 4σ short of the expected depth; shortest period within 90 % of the best stacked S/N | choice | tested (F1 test) |
| Trapezoid fit | `trapezoid`, l. 91 | depth, T14, ingress time | choice (standard trapezoid approximation) | tested (F1 test) |
| Radius | l. 337 | Rp = √depth × R★ (TIC), no dilution correction | sourced (geometric) | untested; the TIC contamination ratio is reported but not applied |
| Stellar-size limit | l. 350 | Rp > 24 R⊕ (about 2.1 RJ) | choice: about the radius of the largest known inflated hot Jupiters | untested |
| Depth limit | l. 351 | depth > 25 % | choice | untested |
| V shape | l. 353 | ingress fraction 2 T12 / T14 > 0.8 | choice | untested |
| Odd/even | l. 297-312, 354 | ≥ 2 fully covered transits of each parity, difference > 3σ | choice | untested |
| Secondary eclipse | l. 330, 355 | phase 0.5 only, > 3σ | choice; eccentric secondaries away from phase 0.5 are not tested here (they are in the binary analysis) | untested |
| Gaia RUWE | l. 357 | > 1.4 | Lindegren et al. 2021 | sourced, untested |
| Gaia eclipsing binary | `gaia_checks` | target in `vari_eclipsing_binary` | Mowlavi et al. 2023 | sourced, untested |
| Nearest Gaia source as target | `gaia_checks`, l. 139 | nearest source within 60 arcsec, not matched by the TIC Gaia identifier | choice | untested; a high-proper-motion star could be matched to a neighbour |
| Blend depth | l. 149-152 | neighbour listed when depth / flux fraction < 1 | sourced (dilution geometry) | untested |
| Folded S/N | l. 358 | < 10 rejects | choice | tested (F1 test) |
| Steps and ramps | l. 358-360 | baseline step > 0.5 depth and > 3σ in ≥ half the epochs | choice | untested |
| Individual epochs | l. 362-366 | ≥ 2 epochs at ≥ 5σ with depth 0.4-2.5 × mean | choice | tested (F1 test) |
| Depth changes | l. 367-373 | χ² probability < 10⁻⁶ and ratio > 2.5 | choice | untested |
| Stellar density | l. 339-342, 374-383 | ρ★ = 3π (a/R★)³ / (G P²) with a/R★ = P (1 + k) / (π T14) for a central transit; rejected below 0.2 × the TIC density | Seager & Mallén-Ornelas 2003 | sourced, untested |
| Equilibrium temperature | l. 343 | Teff √(R★ / 2a) (1 − A)^¼ with A = 0.3 | standard; the albedo 0.3 is a choice | reported only |
| Pulsators | l. 384-394 | Lomb-Scargle above 1.5 / T14, semi-amplitude > 0.5 × depth | Lomb 1976; Scargle 1982 | sourced, untested |
| Sector-synchronised systematics | l. 376-379 | another candidate within 0.3 % in period and 0.3 d in epoch | choice | untested |

## Analysis as an eclipsing binary (`binaries/analysis.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Normalisation | `normalise`, l. 68-89 | 2-d running median with dips > 5σ masked, mask grown by 0.75 widths, three iterations | choice | tested (F2 test) |
| Eclipse width | `eclipse_width` | full width at 10 % of the depth (FWHM at 50 %) | choice | tested (F2 test) |
| Eclipse times | `eclipse_times` | empirical template centred on its axis of mirror symmetry, shift-only χ², parabola at the minimum, error inflated by the reduced χ² | choice | tested (F2 test) |
| Linear ephemeris | `ephemeris` | weighted least squares, errors × √χ²ᵣ | standard | tested (F2 test) |
| Period doubling | l. 227-256 | median odd/even > 5σ, relative difference > 5 %, ≥ 75 % of sectors agree in sign | choice | tested, not doubled (F2 test) |
| Secondary search | l. 273-294 | box of half the primary width at every phase, local minimum, adopted when > 8σ | choice | tested (F2 test) |
| Centroid shift | l. 306-316 | median in eclipse minus out of eclipse, 21″ (TESS) or 3.98″ (Kepler/K2) per pixel | Ricker et al. 2015; Koch et al. 2010 | fixed (F2) |
| QLP Sectors 99-104 | l. 215 | excluded from all eclipse times | MAST QLP documentation (timestamps awaiting redelivery) | sourced |
| Individual eclipses | `eclipses.py` | generalised Gaussian f = a + b(t − tc) − D exp(−\|(t − tc)/s\|^g) | choice | untested |
| e cos ω and e sin ω | `report.py` | Sterne's phase displacement and the width ratio (1 + e sin ω)/(1 − e sin ω) | Sterne 1939; Kopal 1959 (book, no DOI) | sourced, untested |
| Roche lobe in triage | `triage.py:25` | R★ / (0.38 a) for equal masses | Eggleton 1983 (r_L/a = 0.379 at q = 1) | sourced, untested |

## Pixel-level source assignment (`binaries/pixels.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Cutouts | l. 100 | 15 × 15 TESScut pixels | TESScut (MAST) | untested |
| In and out of eclipse | l. 120-121 | in: \|Δt\| < 0.25 D; out: 0.75-1.5 D on both sides | choice | untested |
| Sources modelled | l. 104-105 | Gaia DR3 sources with G < G_target + 4, blends within 8″ merged | choice | untested |
| PSF model | `fit_out` | circular Gaussian, common width 0.4-1.6 px and offset ≤ 1 px, flat background; sector skipped when σ > 1.55 px or offset > 0.95 px | choice (approximation; the TESS PRF is not Gaussian) | untested |
| Decision | l. 136-159 | single source with the lowest χ²; χ² ratio to the runner-up | choice | untested |

## Gaia comparison (`binaries/gaia.py`, `lighttime.py`, `lighttime_scan.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Gaia match | `gaia.py:28` | brightest source within 3″ | choice | untested |
| NSS tables | `gaia.py:10` | two-body orbit, acceleration, non-linear spectroscopic, VIM | Gaia Collaboration, Arenou et al. 2023; Halbwachs et al. 2023 | sourced |
| ω from Thiele-Innes elements | `lighttime.py:19` | ω = ½[atan2(B − F, A + G) + atan2(−B − F, A − G)], ambiguous by 180° (absorbed by the sign of K) | Halbwachs et al. 2023 (Thiele-Innes inversion) | sourced, untested |
| Gaia epoch | l. 18 | t_periastron counted from JD 2457389.0 (J2016.0) | Gaia DR3 documentation | sourced |
| Light-travel time | `shape`, l. 22-28 | (1 − e²)/(1 + e cos ν) sin(ν + ω) + e sin ω, amplitude free | Irwin 1952 | sourced, untested |
| Outer-period scan | `lighttime_scan.py` | sinusoid, 20-3000 d, 20 000 log-spaced periods | choice | tested (F3 test) |

## Catalogue cross-match (`binaries/catalogues.py`, `crossmatch.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Catalogues | `VIZ`, `CATS` | VSX, GCVS, Gaia DR3 EB and classes, Prša et al. 2022, IJspeert et al. 2024, Kostov et al. 2025, SIMBAD, ExoFOP | Watson et al. 2006 (no DOI); Samus et al. 2017; Mowlavi et al. 2023; Rimoldini et al. 2023; Wenger et al. 2000 | sourced |
| Radii | `catalogues.py:26, 35`; `crossmatch.py:72, 75` | 10″ for the star, 63″ (three TESS pixels) for neighbours, 180″ for VSX variables ≥ 2 mag brighter | choice | untested (needs the live services) |

## References checked with Argus

All DOIs in `docs/refs.bib` resolve to the cited works. Argus marks nine as "suspect" only for benign reasons: the Gaia DR3 papers (Gaia Collaboration 2023a, b; Halbwachs; Mowlavi; Rimoldini) and Lindegren et al. 2021 carry their subtitle in the registry as a separate field, SIMBAD (Wenger et al. 2000) likewise, and Stumpe et al. 2012 is registered with its first page only. Ricker et al. is registered by Crossref with the online year 2014; METHODS.md cites it as 2015, the year of the journal volume. Watson et al. 2006 (VSX) and Kopal 1959 have no DOI. The Zenodo record of TIC 402710280.01 (Vora 2026) was checked through OpenAlex.

## Open items

* Most rejection tests of the planet vetting run only in the online recovery tests; each could get a synthetic case like F1.
* The planet radius ignores the TIC contamination ratio, which is reported but not applied; a diluted star gives a radius that is too small.
* `gaia_checks` takes the nearest Gaia source within 60″ instead of the TIC's Gaia identifier.
* The pixel-level assignment uses a Gaussian PSF, which is an approximation of the TESS pixel response function.
