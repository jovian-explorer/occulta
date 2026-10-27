# Methods of occulta: a blind search of TESS, K2 and Kepler photometry for uncatalogued transiting and eclipsing objects

Keshav Aggarwal

## 1. Overview

Space-based photometric surveys have recorded the brightness of several million stars with a precision of tens to hundreds of parts per million, providing an archive in which transiting planets, eclipsing binaries, and other periodic or transient phenomena are recorded whether or not they have been identified. The Transiting Exoplanet Survey Satellite (TESS) alone has produced full-frame images (FFIs) of nearly the entire sky since 2018, while the Kepler and K2 missions observed about half a million stars at higher precision. Although the mission teams and the community have searched these data extensively, the catalogues of planets and eclipsing binaries remain incomplete, in particular for signals with long periods, long durations, or eccentric orbits, which fall outside the assumptions of the standard search pipelines. In this work, we search the public light curves blindly, vet every repeating signal, and report each object that no catalogue lists. The software described here performs every step from the download of the light curves to the final report, and each detection can be reproduced from public data with a single command (Section 9).

## 2. Data

The search uses four sets of light curves, all obtained from the Mikulski Archive for Space Telescopes (MAST).

1. TESS FFI light curves from the MIT Quick-Look Pipeline (QLP; Huang et al. 2020a, 2020b; Kunimoto et al. 2021). We use the systematics-corrected aperture flux (`SYS_RM_FLUX`, or `SAP_FLUX` where the former is absent), since the QLP detrended flux (`KSPSAP_FLUX`) is high-pass filtered and erodes transits longer than a few hours. QLP Sectors 99-104 have not yet been redelivered with corrected timestamps (MAST QLP documentation), and mid-times derived from them are treated as provisional; on stable eclipsing binaries we measure their eclipse times to agree with earlier ephemerides within 2 min in Sectors 99-100 and within about 7 min in Sectors 101-104.
2. TESS 2-minute and 200-second light curves from the Science Processing Operations Center (SPOC; Jenkins et al. 2016), using the Presearch Data Conditioning flux (PDCSAP; Smith et al. 2012; Stumpe et al. 2012).
3. K2 light curves corrected for the spacecraft roll by the K2SFF method (Vanderburg & Johnson 2014), using the best-aperture corrected flux.
4. Kepler long-cadence light curves, quarters 0-17, using PDCSAP flux.

Stellar parameters (TESS magnitude, effective temperature, radius, mass, distance, and the contamination ratio) are taken from the TESS Input Catalog version 8.2 (TIC; Stassun et al. 2019), the K2 Ecliptic Plane Input Catalog, and the Kepler DR25 stellar table (Mathur et al. 2017). Astrometric quality and neighbouring sources are taken from Gaia DR3 (Gaia Collaboration et al. 2023).

The survey is processed one unit at a time (a TESS sector, a K2 campaign, or a Kepler quarter), from the newest to the oldest unit, in the order given by `config.json`: QLP stars with T ≤ 10, every K2 target, every TESS 2-minute target, QLP stars with 10 < T ≤ 11, every Kepler target, and then fainter QLP tiers. Each light curve is searched as soon as it arrives and is deleted afterwards, so that the raw data never exceed 2 GB on disk and the project as a whole never exceeds 30 GB.

## 3. Search of each light curve

Each light curve is normalised by its median, and slow variability is removed with a running median of 2 d width, applied separately to every continuous segment (gaps longer than 0.5 d split segments). Upward outliers beyond 4 times the median absolute deviation (MAD) are removed, while downward points are kept, since they carry the signal. The detrended flux is placed on a uniform time grid at the native cadence.

Individual dips are found with a box search over nine durations (1, 1.5, 2, 3, 4, 5, 6.5, 9, and 13 h). For every box position the mean flux inside the box is compared with the mean flux of two flanking windows, each as long as the box and separated from it by half its width, so that each dip is measured against its local baseline. Two signal-to-noise ratios are computed: a white-noise ratio from the point-to-point scatter and the number of points in the box, and a red-noise ratio from the scatter of all box means of the same duration. The smaller of the two is adopted, which prevents correlated stellar or instrumental noise from producing spurious detections. Up to five dips with a signal-to-noise ratio of at least 7 are recorded per light curve, each masked before the next is sought. The same search is run on the inverted light curve, restricted to boxes away from gaps and segment edges, and the strongest inverted "dip" sets the false-alarm level of that light curve: a real dip must exceed it by a factor of 1.2.

Each dip is characterised by its mid-time, duration, depth, and signal-to-noise ratio, and is flagged when (i) the outer quarters of the box are much shallower than its centre (V shape, ratio below 0.35), (ii) it lies within a few hours of a gap or segment edge, (iii) the background flux or either centroid coordinate departs from its flanking value by more than 6 times its scatter, or (iv) five or more stars on the same camera and CCD show a dip within 0.1 d of the same time, which identifies spacecraft momentum dumps, scattered light, and other common-mode systematics.

Periodic signals with short periods are found with the box least-squares method (BLS; Kovács et al. 2002) as implemented in Astropy (Astropy Collaboration et al. 2022), applied to the light curve binned to 20 min, for periods from 0.5 d to 15 d in a single TESS sector (30 d in a K2 campaign), and durations of 1-6 h. A BLS signal is kept when its signal-to-noise ratio is at least 10, its signal detection efficiency (SDE) is at least 9, it shows at least three transits, and its odd and even depths agree within 3σ. Kepler quarters are searched for individual dips only, since the Kepler pipeline exhausted the short-period signals, while long-period objects that transit once or twice per quarter remain the most likely to have been missed.

## 4. Candidates across units

A star becomes a candidate when (a) two or more dips with a signal-to-noise ratio of at least 9, free of all flags other than proximity to an edge, at least one of them free of every flag, and each above 1.2 times the false-alarm level of its light curve, agree in depth within a factor of 1.6 and in duration within a factor of 2, or (b) a BLS signal passes the criteria of Section 3. Criterion (a) recovers planets and binaries whose period is longer than a sector, which transit once per sector and are therefore invisible to single-sector periodic searches. Signals that coincide in time with a known planet or candidate (TOI, CTOI, Kepler KOI, K2 candidate, or confirmed planet) or that belong to a catalogued eclipsing binary (Prša et al. 2022; Kirk et al. 2016) are set aside at this stage. The lists of known objects are refreshed every week.

## 5. Vetting as a planet

Every light curve of a candidate star is downloaded again and detrended in the same way. The period is resolved by considering every alias P = Δt/n allowed by each pair of dips and testing it against all covered epochs: a period is rejected when it predicts a transit at an epoch that is covered by data and shows no dip, that is, where the measured depth is below 2σ and falls short of the expected depth by more than 4σ, and among the surviving periods the one that stacks the highest combined signal-to-noise ratio is adopted. When the period comes from the box least-squares search, P/2 is tested as well (if P/2 is at least 0.2 d, below the 0.5 d start of the search grid, since a planet with a shorter period can appear there only doubled), because that search often locks on twice the period of a short-period, faint planet; P/2 is kept only if the in-between epochs show a dip and their depth agrees with that of the others within 3σ, so that an eclipsing binary with alternating depths is not halved. A trapezoid model is fitted to the folded light curve to obtain the depth, the total duration T14, and the ingress fraction.

The planet radius follows from the depth and the TIC stellar radius, Rp = R★ √(δ_obs (1 + c)), where c is the TIC contamination ratio (the flux of the neighbouring stars divided by the flux of the target; Stassun et al. 2019); a missing or negative ratio is taken as no dilution, and the radius without the correction is reported next to the corrected one. The candidate is rejected as a planet when any of the following holds: the implied radius exceeds 24 Earth radii; the depth exceeds 25%; the ingress fraction exceeds 0.8 (V shape); the odd and even depths, each measured from at least two fully covered transits, differ by more than 3σ; a secondary eclipse at phase 0.5 exceeds 3σ; Gaia DR3 lists the star as an eclipsing binary or reports a renormalised unit weight error (RUWE) above 1.4 (Lindegren et al. 2021); the folded signal-to-noise ratio is below 10; the baseline differs before and after the dip in at least half of the epochs (steps and ramps); fewer than two fully covered epochs show the transit individually at 5σ with a depth within a factor of 2.5 of the mean; the depth varies between epochs beyond what the errors allow (χ² probability below 10⁻⁶ and a ratio above 2.5, which indicates a diluted blend); the stellar density implied by the transit duration, for a central transit on a circular orbit (Seager & Mallén-Ornelas 2003), is below 0.2 times the catalogue density; another star shows the same period and epoch (a sector-synchronised systematic); fewer than two transits are covered with baseline on both sides; or the star oscillates coherently, at frequencies faster than the transit, with a semi-amplitude above half the transit depth (pulsators such as δ Scuti stars), which is measured with a Lomb-Scargle periodogram (Lomb 1976; Scargle 1982) of the out-of-transit flux. Gaia DR3 is also searched within 60 arcsec for neighbours bright enough to produce the observed depth if they were eclipsed, and these are listed in the report.

## 6. Analysis as an eclipsing binary

Candidates rejected as planets because their signal is stellar (stellar-size radius, a depth above 25%, a V shape, odd/even differences, or a secondary eclipse) are cross-matched with the variable-star catalogues (Section 7). Those that no catalogue lists, and that have no catalogued variable within 63 arcsec (three TESS pixels), are analysed as eclipsing binaries.

1. All sectors are downloaded again. Each sector is normalised by a 2-day running median from which every dip deeper than 5 times the scatter is masked; the mask is grown by 0.75 eclipse widths and the procedure is iterated three times, so that eclipses at any orbital phase, including those of eccentric orbits, are protected from the trend.
2. The width of the primary eclipse is measured from the folded light curve as the full width over which the binned profile is deeper than 10% of its minimum.
3. The mid-time of every well-covered primary eclipse is measured by fitting an empirical template, built from the median folded profile and centred on its axis of mirror symmetry (so that the times refer to mid-eclipse even when the trial epoch is offset), with a shift only; the minimum of the χ² curve is located with a parabola, and the error follows from its curvature, inflated by the reduced χ² when that exceeds unity. A weighted linear ephemeris is fitted to these times, its parameter errors are inflated by the square root of the reduced χ², and the procedure is iterated twice. The residuals (O-C) are kept, together with the coefficient of a quadratic fit when at least six times are available, since eclipse timing variations reveal a third body, apsidal motion, or a changing period.
4. Odd and even primary eclipses are compared inside each sector, which removes the effect of different dilution in different sectors. Only eclipses whose core is fully covered enter the comparison, since an eclipse that falls partly into a data gap would appear shallow, and a sector counts only when it contains at least two eclipses of each parity; the error of the difference is taken from the photometric noise or from the scatter between the depths of eclipses of the same parity, whichever is larger. The period is doubled when the median significance of the difference exceeds 5σ, the median relative difference exceeds 5%, and at least three quarters of the sectors agree in sign, so that a single noisy sector cannot veto the doubling; the deeper eclipse is then assigned to phase zero.
5. The secondary eclipse is searched over all phases outside the primary, with a box of half the primary width, and the most significant dip is adopted. When it exceeds 8σ, its mid-times are measured and fitted in the same way, giving its phase with an error. A secondary phase different from 0.5 measures e cos ω of the orbit (Sterne 1939).
6. The in-eclipse shift of the flux-weighted centroid is measured in each sector, in arcsec (21 arcsec per TESS pixel and 3.98 arcsec per Kepler pixel; Ricker et al. 2015; Koch et al. 2010) and in units of its error. The one-line summary of each binary (`occulta binaries triage`) also compares the catalogue radius of the star with its Roche lobe in an orbit with a companion of equal mass (Eggleton 1983), which rejects giant hosts whose eclipses would require a companion inside the star.

Each binary is then examined by eye in a figure that shows every sample: the longest sector, the full phase curve, the primary and secondary eclipses, and the O-C diagram.

For the systems reported individually, three further analyses are made.

7. Every eclipse covered on both sides is fitted individually with a symmetric generalised-Gaussian profile on a local linear baseline (`occulta binaries eclipses`), which gives its mid-time, depth, and width independently of the template. The secondary phase is measured within each sector relative to the primary eclipses of the same sector, so that it is independent of the long-term ephemeris and of any timing offset common to a sector, and the sector values are averaged with inverse-variance weights (`occulta binaries catalogue`). Eclipse timing variations are accepted only when the QLP and SPOC light curves give the same residuals in every sector they share (`occulta binaries timing`).
8. The source of the eclipses is assigned at the pixel level (`occulta binaries pixels`). For every sector, a 15×15-pixel cutout of the full-frame images is obtained from TESScut. Around each eclipse, the median image within a quarter of the eclipse duration of mid-eclipse is subtracted from the median image between 0.75 and 1.5 durations on both sides, which leaves an image of the light lost in the eclipse. The out-of-eclipse image is modelled as a sum of Gaussian point-spread functions at the positions of all Gaia DR3 sources that are brighter than the target or fainter by less than 4 magnitudes, projected with the world coordinate system of the cutout, with a common width, a common small pointing offset, and a flat background; sources closer than 8 arcsec are merged into one blend, since they cannot be separated in 21-arcsec pixels (Ricker et al. 2015). The difference image is then fitted with the point-spread function of each source in turn, and the source that alone reproduces it with the lowest χ² is taken as the eclipsing star; the ratio of the χ² of the runner-up to that of the best source measures how decisive the assignment is, and sectors in which the out-of-eclipse model fails are discarded.
9. The Gaia DR3 non-single-star catalogue (Gaia Collaboration, Arenou et al. 2023) is searched for an orbital, acceleration, or spectroscopic solution of the brightest Gaia source within 3 arcsec (`occulta binaries gaia`). When an astrometric orbit exists, its Thiele-Innes elements give the semi-major axis and inclination of the photocentre orbit (Halbwachs et al. 2023), and hence a lower limit on the light-travel-time amplitude of an eclipsing pair that moves on it; the eclipse times are then fitted with a light-travel-time term whose period, eccentricity, epoch of periastron, and argument of periastron are fixed at the Gaia values, so that only its amplitude is free (Irwin 1952; `occulta binaries lighttime`), and the fit is repeated with the outer period and the epoch of periastron moved by up to 3 d and 25 d from the Gaia values, which shows how strongly the eclipse times constrain them. Independently, a sinusoidal light-travel-time term is scanned over outer periods from 20 to 3000 d (`occulta binaries lighttime-scan`).

## 7. Catalogue cross-match

A detection is considered new only when the star is absent from the live versions (queried through VizieR, not a cached cross-match service) of the International Variable Star Index (VSX; Watson et al. 2006), the General Catalogue of Variable Stars (GCVS; Samus et al. 2017), the Gaia DR3 eclipsing-binary and variability classifications (Mowlavi et al. 2023; Rimoldini et al. 2023), the TESS eclipsing-binary catalogues of Prša et al. (2022), IJspeert et al. (2024), and Kostov et al. (2025), the variable-star classes of SIMBAD (Wenger et al. 2000), and the ExoFOP lists of TESS Objects of Interest and community TOIs. A radius of 10 arcsec is used for the star itself; catalogued variables between 10 and 63 arcsec are listed as possible contaminants of the TESS aperture.

## 8. Validation

The pipeline recovers known objects without prior information on their periods. From the dips found by the search in TESS data, the vetting stage recovers TIC 402710280.01 at 36.340 d (Vora 2026); from K2 and Kepler data, it recovers K2-3 b at 10.0523 d (published 10.0545 d) and Kepler-10 c at 45.2896 d (published 45.2949 d). Pulsators are rejected, as shown by TIC 192930466, a δ Scuti star whose pulsation beating produced an apparent 5000 ppm dip that passed the earlier tests. The eclipsing-binary analysis was checked on CPD-37 4813 (TIC 181149914), which `occulta reproduce` rediscovers from public data with no prior information; both recoveries form the test suite of the package (`tests/test_recovery.py`).

## 9. Reproducing a detection

```
pip install "git+https://github.com/jovian-explorer/occulta.git"      # or, from a copy of the repository: pip install .
occulta reproduce <TIC> --mode both
```

`occulta reproduce` downloads every light curve of the star from MAST, runs the search of Section 3 on each, stitches the dips into a candidate exactly as the survey does (Section 4) but with no catalogue of known objects loaded, vets it as a planet (Section 5), analyses it as an eclipsing binary when the signal is stellar (Section 6; the option `--mode planets` stops after the planet vetting, and `--mode binaries` analyses the signal as a binary whatever the vetting decides), queries the catalogues (Section 7), and writes all results, figures, the list of data files, and SHA-256 checksums of the code to `reproduce_out/<TIC>/`. The full survey is run with `occulta survey run --minutes 55 --mode both`, which resumes where it stopped, and `occulta survey status` and `occulta survey map` report its progress; with `--mode planets`, the survey searches and vets planet candidates only, and with `--mode binaries` or `--mode both`, it also analyses every uncatalogued stellar signal as an eclipsing binary.

## 10. Limitations

The search relies on the aperture photometry of QLP, SPOC, K2SFF, and Kepler, so that a signal can be attributed to the target only to within the aperture, about one arcminute for TESS; neighbours that could produce each signal are listed, and for the systems that are reported, the source is assigned with difference images of the TESS pixels (Section 6), although stars closer than about 8 arcsec can only be separated by higher-resolution photometry. Eclipses longer than about a day are partly absorbed by the 2-day trend, so that the depths of very long eclipses are lower limits. The period of a binary whose two eclipses are identical cannot be separated from half its value by photometry alone, and such cases are reported at the period that shows two eclipses. QLP Sectors 99-104 await a timestamp correction, and timings from them are provisional; in these sectors, the depths of some eclipses in the QLP light curves also disagree with uniform photometry of the same images, so that depths are taken from the 2-minute light curves or from that photometry.

## References

Astropy Collaboration et al. 2022, ApJ, 935, 167
Eggleton, P. P. 1983, ApJ, 268, 368
Gaia Collaboration, Arenou, F., et al. 2023, A&A, 674, A34
Gaia Collaboration, Vallenari, A., et al. 2023, A&A, 674, A1
Halbwachs, J.-L., et al. 2023, A&A, 674, A9
Huang, C. X., et al. 2020a, RNAAS, 4, 204
Huang, C. X., et al. 2020b, RNAAS, 4, 206
IJspeert, L. W., et al. 2024, A&A, 691, A242
Irwin, J. B. 1952, ApJ, 116, 211
Jenkins, J. M., et al. 2016, Proc. SPIE, 9913, 99133E
Kirk, B., et al. 2016, AJ, 151, 68
Koch, D. G., et al. 2010, ApJL, 713, L79
Kostov, V. B., et al. 2025, ApJS, 279, 50
Kovács, G., Zucker, S., & Mazeh, T. 2002, A&A, 391, 369
Kunimoto, M., et al. 2021, RNAAS, 5, 234
Lindegren, L., et al. 2021, A&A, 649, A2
Lomb, N. R. 1976, Ap&SS, 39, 447
Mathur, S., et al. 2017, ApJS, 229, 30
Mowlavi, N., et al. 2023, A&A, 674, A16
Prša, A., et al. 2022, ApJS, 258, 16
Ricker, G. R., et al. 2015, JATIS, 1, 014003
Rimoldini, L., et al. 2023, A&A, 674, A14
Samus, N. N., et al. 2017, Astronomy Reports, 61, 80
Scargle, J. D. 1982, ApJ, 263, 835
Seager, S., & Mallén-Ornelas, G. 2003, ApJ, 585, 1038
Smith, J. C., et al. 2012, PASP, 124, 1000
Stassun, K. G., et al. 2019, AJ, 158, 138
Sterne, T. E. 1939, MNRAS, 99, 451
Stumpe, M. C., et al. 2012, PASP, 124, 985
Vanderburg, A., & Johnson, J. A. 2014, PASP, 126, 948
Vora, R. 2026, A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280, Zenodo, doi:10.5281/zenodo.23390317
Watson, C. L., Henden, A. A., & Price, A. 2006, SASS, 25, 47
Wenger, M., et al. 2000, A&AS, 143, 9
