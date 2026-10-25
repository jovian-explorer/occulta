# Changelog

## Unreleased

- `occulta planets vet TIC` now vets a star that the survey has not promoted; before, every dip was dropped for want of a candidate depth and the star was always rejected.
- Centroid shifts of Kepler and K2 data are converted with the Kepler pixel (3.98 arcsec) instead of the TESS pixel (21 arcsec), which had overstated them by a factor of 5.3.
- `occulta binaries lighttime` and `occulta binaries lighttime-scan` no longer stop with a NameError.
- The README and METHODS.md install occulta from GitHub (it is not on PyPI); METHODS.md now describes the Gaia sources of the pixel model and the light-travel-time search as the code performs them, and cites the sources of the RUWE limit, the stellar density, the Lomb-Scargle periodogram, the light-travel time, the secondary phase, and the Roche lobe.
- docs/MODEL_AUDIT.md lists every method and threshold with its source and its test; offline tests on synthetic light curves cover the planet vetting, the eclipsing-binary analysis, and the light-travel-time scan.

## 0.0.1 (2026-10-25)

First public pre-release.

- Blind search of TESS QLP and SPOC, K2SFF and Kepler light curves: individual dips and periodic signals in every light curve, dips from different sectors stitched into candidates, and a plan of seven steps worked through from the newest to the oldest sector (`occulta survey run`).
- Vetting of every candidate as a planet: period from all aliases allowed by the dips, trapezoid fit, and rejection tests for stellar radii, V shapes, odd and even depths, secondary eclipses, Gaia astrometry, systematics, depth changes between epochs, stellar density, and coherent pulsations.
- Analysis of stellar signals as eclipsing binaries: ephemerides from a template centred on its axis of symmetry, eclipses fitted one by one, secondary phases measured sector by sector, odd and even eclipses compared with full coverage and empirical scatter, eccentricities from the phase and the duration ratio, and eclipse timing compared between QLP and SPOC (`occulta binaries ...`).
- Assignment of the eclipses to their star with difference images of the TESS pixels and point-spread functions placed at the Gaia DR3 positions (`occulta binaries pixels`).
- Gaia DR3 non-single-star solutions, and fits of the eclipse times with the light-travel time of a Gaia astrometric orbit whose shape and phase are fixed (`occulta binaries gaia`, `occulta binaries lighttime`).
- Live cross-match with VSX, GCVS, Gaia DR3, the TESS eclipsing-binary catalogues, SIMBAD and the TOI lists, through the CfA mirror of VizieR when the CDS server is unreachable.
- `occulta reproduce TIC --mode planets|binaries|both` repeats the analysis of one star from public data and records checksums of the code; the recovery of a published planet candidate and of an eccentric eclipsing binary forms the test suite.
