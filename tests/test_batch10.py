"""Batch 10: a transit of a few cadences passes the V-shape cut, a well-sampled V is still rejected, and a secondary
eclipse must be present in both halves of the data to reject a candidate."""
import numpy as np

from test_cuts import _vet, _local, P, T0, DEP, DUR
from test_offline_science import _transits

SHORT = dict(P=2.3, t0=1001.1, depth=1500e-6, dur=0.85 / 24, cand=(2.3, 1001.1, 0.85, 1500.0))


def test_short_transit_is_not_rejected_for_an_unresolved_ingress(tmp_path, monkeypatch):
    """0.85 h at 10 min cadence is five cadences: a flat-bottomed transit, integrated over each exposure, looks like a V
    and the fit returns an ingress fraction near 1; the V must not be preferred over the flat bottom, so no V-shape cut."""
    from occulta.planets import trapezoid_ss
    def u(t):
        return trapezoid_ss(t - 2.3 * np.round((t - 1001.1) / 2.3), 1001.1, 3000e-6, 0.85 / 24, 0.065 / 24, 10 / 1440)
    status, period, notes = _vet(tmp_path, monkeypatch, **{**SHORT, 'depth': 0.0, 'extra': u, 'cand': (2.3, 1001.1, 0.85, 3000.0)})
    assert 'V-shaped' not in notes, notes
    assert status == 'ready', (status, notes)


def test_well_sampled_v_is_still_rejected(tmp_path, monkeypatch):
    def tri(t):
        x, _ = _local(t)
        return np.clip(1 - np.abs(x) / (2 / 24), 0, None) * -3000e-6
    status, _, notes = _vet(tmp_path, monkeypatch, depth=0.0, extra=tri, cand=(P, T0, 4.0, 3000.0))
    assert status == 'rejected' and 'V-shaped' in notes, notes


def test_secondary_in_both_halves_is_rejected(tmp_path, monkeypatch):
    sec = lambda t: _transits(t, P, T0 + P / 2, 500e-6, 3 / 24)
    status, _, notes = _vet(tmp_path, monkeypatch, extra=sec)
    assert status == 'rejected' and 'secondary eclipse' in notes, notes


def test_secondary_in_one_half_only_does_not_reject(tmp_path, monkeypatch):
    """A secondary that exists only in the later sector (like provisional timestamps) is reported, not used."""
    sec = lambda t: np.where(t > 1026, _transits(t, P, T0 + P / 2, 900e-6, DUR), 0.0)
    status, _, notes = _vet(tmp_path, monkeypatch, extra=sec)
    assert 'secondary eclipse' not in notes.split('not confirmed')[0], notes
    assert status == 'ready', (status, notes)
