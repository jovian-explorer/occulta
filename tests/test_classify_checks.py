"""ZTF robustness, contact-binary duty cycle, blends and note handling in occulta.classify."""
import numpy as np

from occulta import classify as c

P, DUR_H = 41.367, 3.0


def _ztf(depth=0.0, n=2800, outlier=False, seed=1):
    """Flat ZTF-like light curve over 7.5 yr, with a box eclipse of the given depth at phase 0, or with one bad epoch
    inside an eclipse window only (outlier)."""
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(2458200.0, 2460970.0, n))
    f = 1 + rng.normal(0, 0.01, n)
    ph = (t - 2458238.7 + 0.5 * P) % P - 0.5 * P
    inn = np.abs(ph) < DUR_H / 48
    f[inn] -= depth
    if outlier:
        t = np.append(t, 2458238.7 + 3 * P)
        f = np.append(f, 0.40)
    return t, f, np.full(len(t), 0.01)


def test_one_bad_ztf_epoch_is_not_an_eclipse():
    # TIC 427344083: a single ZTF point at 0.40 gave a 41 per cent "undiluted depth" and a false eclipsing binary
    eb, note, zd = c.ztf_test(_ztf(outlier=True), 'TESS', P, 0.0, DUR_H, 12043)
    assert not eb and zd == 0 and 'not recovered' in note


def test_a_real_deep_ztf_eclipse_is_kept():
    eb, note, zd = c.ztf_test(_ztf(depth=0.10, n=6000), 'TESS', P, 0.0, DUR_H, 12043)
    assert 'not recovered' not in note and abs(zd - 1e5) < 2e4


def test_contact_duty_cycle():
    assert c.contact_test(4.0, 0.33)[0]                  # dip over half of a 0.33 d period
    assert not c.contact_test(3.0, 3.5)[0]               # hot Jupiter: 4 per cent of the period
    assert not c.contact_test(None, 1.0)[0]


def test_blend_needs_an_unseparated_neighbour():
    near = [{'sep_arcsec': 6.6, 'dG': -0.15, 'eclipse_needed': 0.002}]
    far = [{'sep_arcsec': 30.0, 'dG': 1.0, 'eclipse_needed': 0.01}]
    faint = [{'sep_arcsec': 5.0, 'dG': 7.0, 'eclipse_needed': 0.8}]
    assert c.blend_test(near)[0] and 'at 6.6 arcsec' in c.blend_test(near)[1]
    assert not c.blend_test(far)[0] and c.blend_test(far, sep_max=42)[0]
    assert not c.blend_test(faint)[0] and not c.blend_test([])[0]


def test_old_class_fragments_are_dropped_and_vetting_notes_kept():
    notes = ('periods that are not multiples of the chosen one appear; TESS/Kepler depth 12043 ppm); '
             'ZTF depth 408379 ppm (TESS/Kepler 12043), odd/even 0.0 sigma, secondary -1.6 sigma); '
             'class: planet (Rp 1.12 R_J from R_star 1.05 R_sun; duration 2.99 h)')
    assert c.strip_class(notes) == ['periods that are not multiples of the chosen one appear']
    assert c.strip_class(None) == []


def test_radius_source_is_named():
    t, why = c.label(10000, 1.0, np.nan, np.nan, src='Gaia DR3')
    assert t == 'planet' and '(Gaia DR3)' in why
    assert '(' not in c.label(10000, 1.0, np.nan, np.nan)[1]
