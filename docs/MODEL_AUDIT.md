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

### F5. `gaia_checks` could take a neighbour for the target
* **Where**: `src/occulta/planets.py` (`gaia_checks`).
* **Problem**: the target was the Gaia source nearest to the TIC position, and the Gaia identifier that the TIC lists was not used. For a high-proper-motion star, or in a crowded field, the nearest source can be a faint neighbour, so that the RUWE and the eclipsing-binary flag of the wrong star were reported and the real star was listed as a blend.
* **Fix**: the source with the TIC's Gaia identifier is the target when it lies within the 60 arcsec; the nearest source remains the fallback when the TIC has no identifier, and the target is skipped in the list of blends by its identifier.
* **Test**: `test_gaia_checks_use_the_tic_gaia_identifier_not_the_nearest_source`. No threshold changed, so that no new citation was needed.

### F6. The planet radius ignored the dilution of the transit; the planet report was not UTF-8
* **Where**: `src/occulta/planets.py` (`vet`).
* **Problem**: Rp = sqrt(depth) R★ used the observed depth, although the TIC contamination ratio c (flux of the neighbours divided by the flux of the target; Stassun et al. 2019) dilutes the transit, so that a diluted star gave a radius that was too small. Separately, `report.md` was written with the system code page, so that on Windows the em dash of its title was not valid UTF-8.
* **Fix**: depth_true = depth_obs (1 + c); a missing or negative c counts as 0. `rp_earth` is the corrected radius (used for the 24 R⊕ limit), `rp_earth_uncorrected` the old value; both appear in `summary.json` and in `report.md`. The report is written as UTF-8.
* **Test**: `test_planet_radius_is_corrected_for_the_tic_contamination_ratio` (c = 0, 0.5, 3, missing, NaN, negative).

## Search of each light curve (`lightcurves.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| QLP flux column | `read_lc`, l. 35 | `SYS_RM_FLUX`, else `SAP_FLUX` (not the high-passed `KSPSAP_FLUX`) | Huang et al. 2020a, b; Kunimoto et al. 2021 | sourced |
| SPOC and Kepler flux | l. 47 | PDCSAP, quality flag 0 | Smith et al. 2012; Stumpe et al. 2012; Jenkins et al. 2016 | sourced |
| TESS-SPOC FFI flux | `read_lc`, same branch as SPOC | PDCSAP, quality flag 0; units 4000 + sector; opt-in plan step | Caldwell et al. 2020; Jenkins et al. 2016 | sourced; tests/test_tess_spoc.py |
| K2 flux | l. 41 | K2SFF `FCOR`, thruster firings (`MOVING`) dropped | Vanderburg & Johnson 2014 | sourced |
| Kepler/K2 time to BTJD | l. 24 | BKJD + 2454833 − 2457000 | mission definitions of BKJD and BTJD | sourced |
| Detrending | `detrend`, l. 65-80 | running median of 2 d per segment, segments split at gaps > 0.5 d, upward outliers > 4 MAD removed | choice | tested: `test_detrend_removes_a_slow_trend_and_keeps_the_dip`, `test_detrend_works_per_segment_across_a_gap_and_clips_only_upward_outliers` |
| Box durations | l. 10 | 1 to 13 h, nine widths | choice | tested: `test_box_search_uses_nine_widths_from_1_to_13_hours_and_matches_each_dip_width` |
| Single-event detection | `single_events` | min(white, red) S/N ≥ 7, at most 5 dips per light curve | choice | tested: `test_single_event_search_finds_a_dip_above_the_threshold`, `test_single_event_search_ignores_noise_and_dips_below_the_threshold`, `test_single_event_search_reports_each_of_several_dips_once` |
| False-alarm level | `process`, l. 231 | strongest dip of the inverted light curve, away from edges | choice (inverted-light-curve test as in common transit-search practice; no single source) | tested: `test_false_alarm_level_is_the_strongest_inverted_dip_away_from_the_edges` (noise stays below 7; a brightening is found; the same bump at a segment edge is ignored; a real dip raises the level because boxes beside it see it in their flank, which is conservative) |
| V-shape flag | l. 166 | depth of outer quarters / centre < 0.35 | choice | tested: `test_characterise_flags_a_v_shaped_dip_and_a_dip_at_a_gap_edge` |
| Background and centroid flags | l. 178 | excursion > 6σ against the flanks | choice | tested: `test_characterise_flags_a_background_excursion` (background only) |
| Common-mode systematics | `survey.mark_systematics` | ≥ 5 stars on one CCD within ±0.1 d | choice | tested: `test_common_mode_flag_needs_five_events_within_a_fifth_of_a_day_on_one_ccd` |
| BLS | `run_bls` | Astropy BLS (Kovács et al. 2002; Astropy Collaboration 2022), 20-min bins, 0.5 d to 15 d (30 d for baselines > 40 d), durations 1-6 h, objective S/N | sourced | untested |
| BLS SDE | l. 205 | (peak − median) / standard deviation of the S/N periodogram | Kovács et al. 2002 define the SDE with the mean of the signal residue; the median used here is slightly larger than the mean for a skewed periodogram, so that the SDE is marginally lower (conservative) | choice, noted |
| BLS acceptance | `survey.promote`, l. 389; `default_config.json` | S/N ≥ 10, SDE ≥ 9, ≥ 3 transits, odd/even < 3σ | choice | tested at each limit: `test_bls_candidates_need_snr_sde_three_transits_and_a_consistent_odd_even` |

## Candidates across units (`survey.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Dip S/N for stitching | `promote`; config `event_snr` | ≥ 9 and ≥ 1.2 × false-alarm level | choice | tested: `test_stitching_needs_strong_similar_dips_above_the_false_alarm_level` |
| Depth and duration agreement | l. 369 | depth ratio 0.6-1.6, duration ratio 0.5-2 | choice | tested (same test, both sides of each limit) |
| Known-object match | `known_match`, l. 333 | ±0.3 d of a TOI, CTOI, KOI, K2 or confirmed-planet ephemeris; any EB-catalogue entry sets the star aside | Prša et al. 2022; Kirk et al. 2016 | sourced, tested: `test_stitching_skips_dips_of_a_known_ephemeris_and_catalogued_binaries` and the BLS test |

## Vetting as a planet (`planets.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Period from aliases | l. 202-272 | every P = Δt/n of a pair of dips; rejected when a covered epoch shows < 2σ and falls > 4σ short of the expected depth; shortest period within 90 % of the best stacked S/N | choice | tested (F1 test) |
| Trapezoid fit | `trapezoid`, l. 91 | depth, T14, ingress time | choice (standard trapezoid approximation) | tested (F1 test) |
| Radius | l. 358 | Rp = √(depth (1 + c)) × R★ (TIC), c = TIC contamination ratio; both radii reported | sourced (geometric; Stassun et al. 2019 for c) | fixed (F6), tested: `test_planet_radius_follows_the_depth_and_the_tic_radius`, `test_planet_radius_is_corrected_for_the_tic_contamination_ratio` |
| Stellar-size limit | l. 350 | Rp > 24 R⊕ (about 2.1 RJ) | choice: about the radius of the largest known inflated hot Jupiters | tested on both sides of the limit: `test_vetting_rejects_a_stellar_size_radius_only_above_24_earth_radii`, and a deep dip on a star without radius: `test_vetting_rejects_a_deep_dip_on_a_star_without_radius` |
| Depth limit | l. 351 | depth > 25 % | choice | tested: `test_vetting_rejects_a_stellar_depth` |
| V shape | l. 353 | ingress fraction 2 T12 / T14 > 0.8 | choice | tested: `test_vetting_rejects_a_v_shaped_transit` (the clean planet of the F1 test stays below the limit) |
| BLS period refit | `planets.refine_period` (`REFINE_SPAN`, `REFINE_TRIALS`) | before the epoch test, P and P/2 of a BLS candidate are refitted on a grid of 97 periods spanning +-12 steps of P dur / B (one duration moved over the baseline B), sampled at a quarter step, times 9 epoch shifts of up to 0.75 durations; the fold with the deepest in-transit mean (0.7 duration wide) wins. Design choice, no external threshold | tested (`test_thresholds.py`, input periods 1e-4 and 2e-4 d off over 2,000 d; fails without it) |
| P/2 trial | `planets.vet` (`HALF_ODD_EVEN_SIGMA`) | for a BLS candidate, P/2 (if ≥ `P_HALF_MIN` = 0.2 d, about the shortest known planet periods and below the 0.5 d start of the BLS grid, otherwise LHS 3844 b at 0.4629 d could never be reached) competes with P; it needs no absent covered epoch and a depth difference between its odd and even epochs ≤ 3σ (the same limit as the odd/even test) | choice; reason: BLS locked on 2 × P for LHS 3844 b (0.926 d instead of 0.4629 d), and a true P/2 stacks √2 more signal than P, while an eclipsing binary at P/2 fails the parity test | tested: `test_vetting_halves_a_bls_period_that_is_twice_the_true_one`, `test_vetting_does_not_halve_alternating_depths`, `test_vetting_keeps_the_true_bls_period` |
| Odd/even | l. 297-312, 354 | ≥ 2 fully covered transits of each parity, difference > 3σ | choice | tested: `test_vetting_rejects_alternating_depths` |
| Secondary eclipse | l. 330, 355 | phase 0.5 only, > 3σ | choice; eccentric secondaries away from phase 0.5 are not tested here (they are in the binary analysis) | tested: `test_vetting_rejects_a_secondary_eclipse` |
| Gaia RUWE | l. 357 | > 1.4 | Lindegren et al. 2021 | sourced, tested at 1.40 and 1.41 with a mocked archive |
| Gaia eclipsing binary | `gaia_checks` | target in `vari_eclipsing_binary` | Mowlavi et al. 2023 | sourced, tested with a mocked archive |
| Target source | `gaia_checks`, l. 139 | the source with the TIC's Gaia identifier among those within 60 arcsec, else the nearest | choice | fixed (F5), tested |
| Blend depth | l. 149-152 | neighbour listed when depth / flux fraction < 1 | sourced (dilution geometry) | tested: a neighbour 2 mag fainter explains a 5 % dip, not a 20 % one |
| Folded S/N | l. 358 | < 10 rejects | choice | tested (F1 test) |
| Steps and ramps | l. 358-360 | baseline step > 0.5 depth and > 3σ in ≥ half the epochs | choice | tested: `test_vetting_rejects_a_step_instead_of_a_transit` |
| Individual epochs | l. 362-366 | ≥ 2 epochs at ≥ 5σ with depth 0.4-2.5 × mean | choice | tested (F1 test) |
| Depth changes | l. 367-373 | χ² probability < 10⁻⁶ and ratio > 2.5 | choice | tested: `test_vetting_rejects_depths_that_change_between_sectors` (4000 against 1300 ppm) |
| Stellar density | l. 339-342, 374-383 | ρ★ = 3π (a/R★)³ / (G P²) with a/R★ = P (1 + k) / (π T14) for a central transit; rejected below 0.2 × the TIC density | Seager & Mallén-Ornelas 2003 | sourced, tested: the formula against the summary value, and both density branches (ratio 0.10 and 0.013) |
| Equilibrium temperature | l. 343 | Teff √(R★ / 2a) (1 − A)^¼ with A = 0.3 | standard; the albedo 0.3 is a choice | reported only |
| Pulsators | l. 384-394 | Lomb-Scargle above 1.5 / T14, semi-amplitude > 0.5 × depth | Lomb 1976; Scargle 1982 | sourced, tested at 0.7 and 0.2 × the depth |
| Sector-synchronised systematics | l. 376-379 | another candidate within 0.3 % in period and 0.3 d in epoch | choice | tested: `test_vetting_rejects_a_signal_shared_with_another_candidate_within_the_tolerance` |

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
| Individual eclipses | `eclipses.py` | generalised Gaussian f = a + b(t − tc) − D exp(−\|(t − tc)/s\|^g) | choice | tested: `test_generalised_gaussian_fit_recovers_time_depth_and_width` (time, depth, width, noise level; an eclipse cut at its centre is refused) |
| e cos ω and e sin ω | `report.py` | Sterne's phase displacement and the width ratio (1 + e sin ω)/(1 − e sin ω) | Sterne 1939; Kopal 1959 (book, no DOI) | sourced, tested: `test_orbit_solution_recovers_eccentricity_and_argument_of_periastron` (circular case, first-order limit, four (e, ω) pairs recovered, impossible phase refused) |
| Roche lobe in triage | `triage.py:25` | R★ / (0.38 a) for equal masses | Eggleton 1983 (r_L/a = 0.379 at q = 1) | sourced, tested: `test_roche_ratio_of_a_twin_mass_orbit` (Kepler's third law in SI, Eggleton's 0.379, missing mass) |

## Pixel-level source assignment (`binaries/pixels.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Cutouts | l. 100 | 15 × 15 TESScut pixels | TESScut (MAST) | untested (needs the live service; the analysis of a cutout file is tested) |
| In and out of eclipse | l. 120-121 | in: \|Δt\| < 0.25 D; out: 0.75-1.5 D on both sides | choice | tested: `test_pixel_model_finds_which_star_eclipses` |
| Sources modelled | l. 104-105 | Gaia DR3 sources with G < G_target + 4, blends within 8″ merged | choice | tested: `test_sources_closer_than_8_arcsec_are_merged_into_one_blend`; the G cut is applied in `main` and not tested separately |
| PSF model | `fit_out` | circular Gaussian, common width 0.4-1.6 px and offset ≤ 1 px, flat background; sector skipped when σ > 1.55 px or offset > 0.95 px | choice (approximation; the TESS PRF is not Gaussian) | tested on a synthetic Gaussian field (the limits for a skipped sector are not) |
| Decision | l. 136-159 | single source with the lowest χ²; χ² ratio to the runner-up | choice | tested: an eclipse on the target and on a neighbour 3 px away give the right source, and the lost fraction (25 %) is recovered |

## Gaia comparison (`binaries/gaia.py`, `lighttime.py`, `lighttime_scan.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Gaia match | `gaia.py:28` | brightest source within 3″ | choice | tested with a mocked archive: `test_gaia_match_asks_for_the_brightest_source_within_3_arcsec` (query radius, position, brightness order, only graded A and B stars) |
| NSS tables | `gaia.py:10` | two-body orbit, acceleration, non-linear spectroscopic, VIM | Gaia Collaboration, Arenou et al. 2023; Halbwachs et al. 2023 | sourced |
| ω from Thiele-Innes elements | `lighttime.py:19` | ω = ½[atan2(B − F, A + G) + atan2(−B − F, A − G)], ambiguous by 180° (absorbed by the sign of K) | Halbwachs et al. 2023 (Thiele-Innes inversion) | sourced, tested: ω = 40° recovered from A, B, F, G of an orbit with Ω = 70°, i = 60° |
| Gaia epoch | l. 18 | t_periastron counted from JD 2457389.0 (J2016.0) | Gaia DR3 documentation | sourced |
| Light-travel time | `shape`, l. 22-28 | (1 − e²)/(1 + e cos ν) sin(ν + ω) + e sin ω, amplitude free | Irwin 1952 | sourced, tested: a 3 min term of e = 0.4 is fitted with a χ² below 10⁻³ of the linear fit and K = 3.0 min (`test_lighttime_fit_recovers_a_gaia_shaped_term_and_omega_from_thiele_innes`) |
| Outer-period scan | `lighttime_scan.py` | sinusoid, 20-3000 d, 20 000 log-spaced periods | choice | tested (F3 test) |

## Catalogue cross-match (`binaries/catalogues.py`, `crossmatch.py`)

| Item | Where | Value | Source | Status |
|---|---|---|---|---|
| Catalogues | `VIZ`, `CATS` | VSX, GCVS, Gaia DR3 EB and classes, Prša et al. 2022, IJspeert et al. 2024, Kostov et al. 2025, SIMBAD, ExoFOP | Watson et al. 2006 (no DOI); Samus et al. 2017; Mowlavi et al. 2023; Rimoldini et al. 2023; Wenger et al. 2000 | sourced |
| Radii | `catalogues.py:26, 35`; `crossmatch.py:72, 75` | 10″ for the star, 63″ (three TESS pixels) for neighbours, 180″ for VSX variables ≥ 2 mag brighter | choice | untested (needs the live services) |

## References checked with Argus

All DOIs in `docs/refs.bib` resolve to the cited works. Argus marks ten as "suspect" only for benign reasons: the Gaia DR3 papers (Gaia Collaboration 2023a, b; Halbwachs; Mowlavi; Rimoldini) and Lindegren et al. 2021 carry their subtitle in the registry as a separate field, SIMBAD (Wenger et al. 2000) likewise, and Stumpe et al. 2012 is registered with its first page only. Ricker et al. (2015) is registered by Crossref with the online date of 24 October 2014 and the page range 1-10; text and `refs.bib` both give 2015, the year of the journal volume, and the article number 014003, so Argus marks its first page. Watson et al. 2006 (VSX) and Kopal 1959 have no DOI. The Zenodo record of TIC 402710280.01 (Vora 2026) was checked through OpenAlex.

## Open items

* Every threshold of the planet vetting, the stitching, the Gaia checks, the pixel model, the light-travel-time fit, the box durations, the false-alarm level, the eclipse fit, e cos ω and e sin ω, the Roche ratio and the Gaia match radius now has an offline test (batches 3 and 4). Still untested offline: the TESScut download and the catalogue radii (both need the live services).
* The pixel-level assignment uses a Gaussian PSF, which is an approximation of the TESS pixel response function.
