# ZTF light curves in occulta: design note

This note describes `src/occulta/ztf.py`, which adds a first, unwired path for searching Zwicky Transient Facility (ZTF) photometry for transits and eclipses. The module is deliberately small and is called from `occulta ztf` (see Wiring below).

## Data access

ZTF light curves are downloaded from the public Infrared Science Archive (IRSA) light-curve interface, which accepts either a ZTF object identifier or a sky position and returns a table in CSV, IPAC or VOTable form ([API description](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html)). The public collection needs no login, whereas proprietary data require a password, and the same page limits the search radius of a position query to 0.1667 degrees, which `ztf.py` enforces before any request is sent. The same page asks users who need more than about 1000 light curves to read the HATS catalogue instead of calling this interface, so a survey-wide ZTF search will need a different access path than the one written here. The ZTF survey itself and its cadence are described by Bellm et al. ([2019](https://doi.org/10.1088/1538-3873/aaecbe)), and the data processing and the light-curve products by Masci et al. ([2019](https://doi.org/10.1088/1538-3873/aae8ac)).

The module requests CSV, parses it with the standard library, and uses the columns `oid`, `hjd`, `mag`, `magerr`, `catflags` and `filtercode`. The meaning of these columns, and of the bits of `catflags`, is given in the [ZTF Science Data System Explanatory Supplement](https://irsa.ipac.caltech.edu/data/ZTF/docs/ztf_explanatory_supplement.pdf); the column names were not re-checked against a live response in this session, because the tests run offline, so the first online run should confirm them.

## Filtering and flux conversion

A point is kept only if `catflags` equals zero, its magnitude error is finite and not larger than 0.1 mag, and its time and magnitude are finite. The IRSA page describes `BAD_CATFLAGS_MASK` as a way to drop points with chosen flag bits set, and the stricter condition `catflags == 0` is applied here on the client, after download, so that the choice stays visible in the code.

The three ZTF filters (`zg`, `zr`, `zi`) have different zero points and colour terms, and an eclipse of a star with a hot companion can have a different depth in each. Each filter is therefore normalised to its own median flux before the filters are joined, and the combined series is searched; a single filter can be searched alone with `band=`. For each filter only the object identifier with the most good points is used, since the same star observed in two ZTF fields is stored as two objects with different zero points.

Magnitudes are converted with flux = 10^(-0.4 (m - median m)), and the flux error is 0.4 ln(10) times the flux times the magnitude error. Times are heliocentric Julian dates (`hjd`) minus 2457000, which is the BTJD time base that occulta uses for TESS, K2 and Kepler. The difference between heliocentric and barycentric time is up to a few minutes, which is negligible against the one-hour or longer durations searched here.

## Search method

The search uses `astropy.timeseries.BoxLeastSquares`, since astropy is already a dependency of occulta. The period grid runs from 0.2 d to 30 d (and never beyond the baseline divided by three) and is spaced uniformly in the logarithm of the frequency, with a step of d(ln f) = duration / (2 x baseline), so that the phase slip over the whole baseline stays below half a transit duration. For a baseline of about three years and durations from 1 h to 11 h this gives between about 2 x 10^4 (longest duration) and 3 x 10^5 (shortest duration) trial periods. For each of seven trial durations the thirty best peaks in signal-to-noise are kept, and a candidate must have at least five points inside the box and at least three distinct transit epochs. These two requirements are what prevents a single bad night from being reported as a periodic signal. The reported quantities are period, epoch, duration, depth, signal-to-noise ratio, `n_in_transit` and `n_transits`, and only dimming events are returned.

## Limits of sparse cadence

ZTF observes a given field at most a few times per night and, in the public survey, typically once every few nights (Bellm et al. [2019](https://doi.org/10.1088/1538-3873/aaecbe)). The consequences for a transit search are the following.

- Sampling at about one point per night makes periods near one day, and near simple fractions and multiples of a day, ambiguous, because the phase coverage is nearly identical for all of them; `alias_flags` marks such periods and they must not be taken as the true period without follow-up.
- A transit of a few hours is sampled only if an exposure happens to fall inside it, so that with 1000 nights only a few points per transit epoch, often one, constrain the depth; hence only deep events are detectable, in practice eclipses of binaries and giant planets around bright stars whose photometric error is well below the depth.
- The per-point precision of ZTF photometry is limited to the percent level for all but the brightest stars (see the photometric performance in Masci et al. [2019](https://doi.org/10.1088/1538-3873/aae8ac); the exact figure per magnitude should be taken from that paper), so that Earth-size and sub-Neptune planets cannot be found, and the depth sensitivity must be measured with injections before any non-detection is interpreted.
- Systematic trends and flagged nights are not removed by a running median here, because with nightly sampling a filter cannot distinguish a long eclipse from a trend; red noise therefore remains a source of false positives, and no candidate from this module should be treated as a detection before vetting.

## Wiring (batch 10)

- `sources.py`: ZTF unit codes are 5000 + filter id, that is 5001 (zg), 5002 (zr) and 5003 (zi), one unit per filter because the zero points differ; `unit_name`, `fmt_of_unit` ('ztf') and `units('ZTF')` know them, and `unit_targets` returns None because IRSA offers no per-unit target list.
- `cli.py`: `occulta ztf RA DEC [--radius-arcsec R] [--pmin D] [--pmax D] [--band zg|zr|zi]` fetches and searches one position. `{"mission": "ZTF"}` in `plan` is accepted and the job is marked missing with a log line pointing to this command.
- `tests/test_ztf.py` has offline tests of the units, the command and the plan, and a live test (`OCCULTA_ONLINE=1`, about one minute) that checks the IRSA column names. A query of 2 arcsec at the catalogue position of WASP-12 returned no rows on 28 Oct 2026 (the ZTF object lies 8 arcsec away), so use a radius that matches the field.

## Still to be done

- A star identifier namespace for ZTF objects and the cross-match to Gaia, to attach a ZTF light curve to a catalogued star.
- The vetting steps (odd-even, secondary eclipse, catalogue lookup) of the rest of occulta for ZTF candidates, and an injection test of the detection limit; the false-positive rate on pure noise has been tested only on one synthetic series.
