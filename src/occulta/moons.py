"""Exomoons: transit-timing (TTV) and transit-duration (TDV) variations of every vetted transiting object, their
periodograms, false-alarm probabilities and phase relation, an orbital-sampling-effect (OSE) stack of the folded light
curve around the transit, and upper limits on the mass and radius of a moon where nothing is found.

    occulta moons search [TIC ...] [--include-rejected]
                                  every planet that the vetting marked ready (or the stars given), from the transit
                                  windows that the vetting keeps in findings/<star>/transit_windows.npz; nothing is
                                  downloaded. --include-rejected also searches the transiting companions that the
                                  vetting rejected (stellar radius, V shape, ...).

Results go to findings/<star>/moons.md, moons.json and moons.png, and to the table moon of the database. The survey runs
the search on every newly vetted planet when "search_moons": true is set in config.json.
"""
import argparse, json, time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as _chi2

from . import lightcurves as lcs
from .planets import trapezoid_ss

G = 6.674e-11
RSUN, MSUN, REARTH, MEARTH, AU = 6.957e8, 1.989e30, 6.371e6, 5.972e24, 1.496e11
MIN_TRANSITS = 5        # fewer timed transits give no periodogram
MOON_FAP = 1e-3         # a TTV, TDV or OSE signal counts when its false-alarm probability is below this
HILL_STABLE = 0.4895    # widest stable prograde moon orbit, in Hill radii (Domingos et al. 2006)
N_FLIPS = 4000          # sign flips of whole transit windows for the false-alarm probability of the OSE stack
LIM_CL = 2.71           # one-sided 95 % limit from the profile likelihood (delta -2 ln L)
PHASE_SIGMA = 2.0       # TTV and TDV in quadrature (|dphi| = 90 deg) within this many sigma
PHASE_P = 1e-2          # the second series must show the same frequency with a noise probability below this
OSE_SNR = 7.0           # an OSE dip also needs this S/N, since sign flips of a few windows cannot reach small FAPs
OSE_FAP = 1e-2          # ... and a sign-flip false-alarm probability below this


# ---------------------------------------------------------------- planet mass and orbit
def planet_mass(rp_earth):
    """Planet mass (Earth masses) from its radius with the probabilistic mass-radius relation of Chen & Kipping (2017),
    mean branch: R = 1.008 M^0.279 below 2.04 Earth masses, R proportional to M^0.589 up to 0.414 Jupiter masses, and
    one Jupiter mass for larger radii, where the radius no longer constrains the mass."""
    if not rp_earth or not np.isfinite(rp_earth):
        return np.nan
    r_break1 = 1.008 * 2.04 ** 0.279
    if rp_earth <= r_break1:
        return float((rp_earth / 1.008) ** (1 / 0.279))
    m = 2.04 * (rp_earth / r_break1) ** (1 / 0.589)
    return float(m) if m <= 0.414 * 317.8 else 317.8


def moon_ttv_rms(ms, mp, a_s, a_p, P):
    """rms TTV (s) of a planet of mass mp with a moon of mass ms on a circular, coplanar orbit of radius a_s (m), the
    planet orbiting at a_p (m) in P (s): the barycentre offset a_s ms/mp divided by the orbital speed 2 pi a_p / P,
    times 1/sqrt(2) for moon phases sampled at random (Sartoretti & Schneider 1999; Kipping 2009a)."""
    return a_s * ms / mp * P / (2 * np.pi * a_p) / np.sqrt(2)


def moon_tdv_rms(ms, mp, a_s, a_p, P, T):
    """rms TDV (s) of a transit of duration T (s): the planet's speed about the barycentre, ms sqrt(G / (a_s mp)),
    relative to its orbital speed, times T, times 1/sqrt(2) for random moon phases (Kipping 2009a)."""
    return T * ms * np.sqrt(G / (a_s * mp)) / (2 * np.pi * a_p / P) / np.sqrt(2)


# ---------------------------------------------------------------- timing
def transit_times(t, y, P, t0, t14, dep, t12):
    """Mid-time and duration of every transit that is covered at 80 % and has baseline on both sides, fitted with the
    folded trapezoid of fixed depth and ingress fraction on a local straight line: [(epoch, tc, e_tc, t14, e_t14)] (d)."""
    out = []
    frac = t12 / t14
    for n in np.unique(np.round((t - t0) / P)):
        tc0 = t0 + n * P
        m = np.abs(t - tc0) < t14 / 2 + max(t14, 0.1)
        tw, yw = t[m], y[m]
        if tw.size < 12:
            continue
        cad = float(np.median(np.diff(tw)))
        x = tw - tc0
        if np.sum(np.abs(x) < t14 / 2) * cad < 0.8 * t14 or np.sum(x < -0.6 * t14) < 5 or np.sum(x > 0.6 * t14) < 5:
            continue
        res = lambda p: yw - (p[2] + p[3] * x + trapezoid_ss(x, p[0], dep, p[1], frac * p[1], cad))
        lo, hi = [-0.5 * t14, 0.5 * t14, -1, -1], [0.5 * t14, 1.5 * t14, 1, 1]
        f = least_squares(res, [0, t14, 0, 0], bounds=(lo, hi), x_scale=[t14 / 20, t14 / 20, dep, dep])
        if abs(f.x[0]) > 0.45 * t14 or not 0.55 * t14 < f.x[1] < 1.45 * t14:
            continue                      # at a bound: no transit to time
        s2 = np.sum(f.fun ** 2) / max(1, yw.size - 4)
        cov = np.linalg.pinv(f.jac.T @ f.jac) * s2
        if not (cov[0, 0] > 0 and cov[1, 1] > 0):
            continue
        out.append((int(n), tc0 + f.x[0], float(np.sqrt(cov[0, 0])), float(f.x[1]), float(np.sqrt(cov[1, 1]))))
    return out


def jitter_limit(r, e, cl=LIM_CL):
    """(best, 95 % upper limit) of an extra scatter s of r beyond the errors e: r_i ~ N(0, e_i^2 + s^2), profile
    likelihood, one-sided."""
    top = 5 * max(np.std(r), np.median(e)) + 1e-12
    s = np.linspace(0, top, 4001)
    v = e[None, :] ** 2 + s[:, None] ** 2
    nll = np.sum(np.log(v) + r[None, :] ** 2 / v, axis=1)
    i = int(np.argmin(nll))
    above = np.where((nll - nll[i] > cl) & (s > s[i]))[0]
    return float(s[i]), float(s[above[0]] if above.size else s[-1])


def periodogram(n, r, e):
    """Generalised Lomb-Scargle periodogram of a series sampled once per transit (frequencies in cycles per epoch up to
    the Nyquist frequency 0.5; a moon's period is aliased into this band) and the false-alarm probability of its highest
    peak (Baluev 2008). Returns (frequency, power, FAP)."""
    from astropy.timeseries import LombScargle
    span = max(n.max() - n.min(), 2)
    ls = LombScargle(n, r, e)
    fq, pw = ls.autopower(minimum_frequency=1 / span, maximum_frequency=0.5, samples_per_peak=10)
    i = int(np.argmax(pw))
    fap = float(ls.false_alarm_probability(pw[i], minimum_frequency=1 / span, maximum_frequency=0.5, method='baluev'))
    return float(fq[i]), float(pw[i]), min(1.0, fap)


def sine_phase(n, r, e, f):
    """Amplitude, phase (rad), phase error and the probability that a sinusoid this strong arises from noise
    (chi-square with two degrees of freedom) of r = A cos(2 pi f n - phi) + c, by weighted least squares; the covariance
    is scaled by the reduced chi-square of the fit when that exceeds one."""
    X = np.vstack([np.cos(2 * np.pi * f * n), np.sin(2 * np.pi * f * n), np.ones_like(n, dtype=float)]).T
    W = 1 / e ** 2
    C = np.linalg.pinv(X.T @ (X * W[:, None]))
    a, b, c = C @ (X.T @ (W * r))
    c2 = float(np.sum(W * (r - X @ [a, b, c]) ** 2))
    c2_flat = float(np.sum(W * (r - np.sum(W * r) / np.sum(W)) ** 2))
    C = C * max(1.0, c2 / max(1, r.size - 3))                                 # errors from the scatter
    A = np.hypot(a, b)
    ephi = np.sqrt(max(b * b * C[0, 0] + a * a * C[1, 1] - 2 * a * b * C[0, 1], 0)) / max(A * A, 1e-300)
    return float(A), float(np.arctan2(b, a)), float(ephi), float(_chi2.sf(max(c2_flat - c2, 0), 2))


# ---------------------------------------------------------------- orbital sampling effect
def ose_template(tau, tau_s, h):
    """Fraction of the phases of a moon on a circular orbit, seen edge-on at a projected distance tau_s cos(phi) (in time)
    from the planet, for which the moon is in transit at time tau from mid-transit; h is half the transit duration
    (Heller 2014)."""
    hi = np.arccos(np.clip((tau - h) / tau_s, -1, 1))
    lo = np.arccos(np.clip((tau + h) / tau_s, -1, 1))
    return (hi - lo) / np.pi


def ose_stack(tau, r, win, h, taus, sig, rng, n_flips=N_FLIPS):
    """Fit the OSE template of every trial orbit tau_s to the residuals r outside the planetary transit; returns
    (depth, error, S/N per tau_s, best index, false-alarm probability of the best S/N from sign flips of whole windows
    of the residuals left by the best template)."""
    use = np.abs(tau) > h * 1.05
    tau, r, win = tau[use], r[use], win[use]
    M = np.array([ose_template(tau, ts, h) for ts in taus])          # (k, points)
    B = np.sum(M ** 2, axis=1) / sig ** 2
    ok = B > 0
    uw, inv = np.unique(win, return_inverse=True)
    Aw = np.zeros((len(taus), uw.size))
    for k in range(len(taus)):
        Aw[k] = np.bincount(inv, M[k] * r, minlength=uw.size) / sig ** 2
    A = Aw.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        dep = np.where(ok, -A / B, np.nan)
        err = np.where(ok, 1 / np.sqrt(B), np.nan)
    snr = dep / err
    best = int(np.nanargmax(snr)) if np.any(np.isfinite(snr)) else 0
    # null: the residuals after the best template is removed, with whole windows flipped in sign, so that a real dip
    # does not widen the noise distribution against which it is judged
    d0 = np.nan_to_num(dep[best]) if snr[best] > 0 else 0.0
    Ar = np.array([Aw[k] + d0 * np.bincount(inv, M[k] * M[best], minlength=uw.size) / sig ** 2 for k in range(len(taus))])
    s = rng.choice([-1.0, 1.0], size=(n_flips, uw.size))
    with np.errstate(invalid='ignore', divide='ignore'):
        null = np.nanmax(-(s @ Ar.T) / np.sqrt(np.where(ok, B, np.nan)), axis=1)
    fap = float((1 + np.sum(null >= snr[best])) / (n_flips + 1))
    return dep, err, snr, best, fap


# ---------------------------------------------------------------- one object
def analyse(t, y, summary, rng=None):
    """The moon search of one vetted transiting object: t, y are the transit windows (BTJD, detrended relative flux),
    summary the vetting summary (period_d, t0_bjd, t14_h, depth_ppm, ingress_frac, rp_earth, star). Returns a dict."""
    rng = rng or np.random.default_rng(0)
    P, t0 = summary['period_d'], summary['t0_bjd'] - 2457000.0
    t14, dep = summary['t14_h'] / 24, summary['depth_ppm'] / 1e6
    t12 = min(max(summary.get('ingress_frac') or 0.2, 0.02), 1.0) * t14 / 2
    st = summary.get('star') or {}
    mstar = st.get('mass') if st.get('mass') and np.isfinite(st.get('mass')) else 1.0
    rstar = st.get('rad') if st.get('rad') and np.isfinite(st.get('rad')) else 1.0
    rp = summary.get('rp_earth')
    rp = rp if rp and np.isfinite(rp) else np.sqrt(dep) * rstar * RSUN / REARTH
    mp = planet_mass(rp)
    Ps = P * 86400
    a_p = (G * mstar * MSUN * Ps ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    v_p = 2 * np.pi * a_p / Ps
    r_hill = a_p * (mp * MEARTH / (3 * mstar * MSUN)) ** (1 / 3)
    out = {'tic': summary.get('tic'), 'period_d': P, 'mp_earth_assumed': mp, 'rp_earth': rp, 'a_p_au': a_p / AU,
           'hill_radius_rp': r_hill / (rp * REARTH), 'mstar': mstar, 'rstar': rstar, 'signals': [], 'notes': []}
    if summary.get('star', {}).get('mass') is None:
        out['notes'].append('no stellar mass in the catalogue: 1 solar mass and radius assumed')
    if rp > 14.3:
        out['notes'].append('Jupiter-size radius: the planet mass is set to one Jupiter mass, the radius does not constrain it')

    # ---- timing
    tt = transit_times(t, y, P, t0, t14, dep, t12)
    out['n_transits'] = len(tt)
    out['transits'] = [{'epoch': n, 'tc_btjd': tc, 'e_tc_min': e * 1440, 't14_h': d * 24, 'e_t14_min': ed * 1440}
                       for n, tc, e, d, ed in tt]
    if len(tt) >= MIN_TRANSITS:
        n = np.array([x[0] for x in tt], float)
        tc, etc = np.array([x[1] for x in tt]), np.array([x[2] for x in tt])
        d, ed = np.array([x[3] for x in tt]), np.array([x[4] for x in tt])
        Wt = 1 / etc ** 2
        X = np.vstack([np.ones_like(n), n]).T
        T0, Pf = np.linalg.solve(X.T @ (X * Wt[:, None]), X.T @ (Wt * tc))
        oc = tc - (T0 + Pf * n)
        tdv = d - np.sum(d / ed ** 2) / np.sum(1 / ed ** 2)
        for name, r, e, dof in (('ttv', oc, etc, n.size - 2), ('tdv', tdv, ed, n.size - 1)):
            c2 = float(np.sum((r / e) ** 2))
            f, pw, fap = periodogram(n, r, e)
            best, lim = jitter_limit(r, e)
            out[name] = {'chi2': c2, 'dof': int(dof), 'p_chi2': float(_chi2.sf(c2, dof)), 'rms_s': float(np.std(r) * 86400),
                         'jitter_s': best * 86400, 'limit95_s': lim * 86400, 'freq_per_epoch': f, 'power': pw, 'fap': fap}
        f_ref = out['ttv']['freq_per_epoch'] if out['ttv']['fap'] <= out['tdv']['fap'] else out['tdv']['freq_per_epoch']
        A1, p1, e1, q1 = sine_phase(n, oc, etc, f_ref)
        A2, p2, e2, q2 = sine_phase(n, tdv, ed, f_ref)
        dphi = (np.degrees(p2 - p1) + 180) % 360 - 180
        edphi = float(np.degrees(np.hypot(e1, e2)))
        out['phase'] = {'freq_per_epoch': f_ref, 'ttv_amp_s': A1 * 86400, 'tdv_amp_s': A2 * 86400,
                        'p_ttv_at_freq': q1, 'p_tdv_at_freq': q2, 'dphi_deg': float(dphi), 'e_dphi_deg': edphi,
                        'quadrature': bool(abs(abs(dphi) - 90) <= PHASE_SIGMA * edphi + 1e-9)}
        sig_ttv, sig_tdv = out['ttv']['fap'] < MOON_FAP, out['tdv']['fap'] < MOON_FAP
        # a moon: one series significant over the whole periodogram, the other at that same frequency, in quadrature
        both = (sig_ttv and q2 < PHASE_P) or (sig_tdv and q1 < PHASE_P)
        if both and out['phase']['quadrature']:
            out['signals'].append('moon_timing')
        elif sig_ttv or sig_tdv:
            out['signals'].append('timing_variation')
            out['notes'].append('TTV or TDV significant, but not both in quadrature as a moon requires: another planet, '
                                'spots, or systematics are more likely')
        # ---- mass limits: largest moon mass allowed at each orbit between 2 Rp and the stability limit
        a_s = np.geomspace(2 * rp * REARTH, max(HILL_STABLE * r_hill, 2.0001 * rp * REARTH), 12)
        ms_ttv = out['ttv']['limit95_s'] / moon_ttv_rms(MEARTH, mp * MEARTH, a_s, a_p, Ps)    # signals of a moon of 1 Me
        ms_tdv = out['tdv']['limit95_s'] / moon_tdv_rms(MEARTH, mp * MEARTH, a_s, a_p, Ps, t14 * 86400)
        ms = np.minimum(ms_ttv, ms_tdv)                 # Earth masses
        out['mass_limit'] = {'a_s_rp': (a_s / (rp * REARTH)).tolist(), 'ms_earth': ms.tolist(),
                             'ms_earth_any_orbit': float(np.max(ms)), 'ms_earth_widest_orbit': float(ms[-1])}
        if HILL_STABLE * r_hill < 2 * rp * REARTH:
            out['notes'].append('the Hill sphere is too small for a stable moon')
        out['_series'] = (n, oc, etc, tdv, ed)

    # ---- orbital-sampling-effect stack, relative to each transit's own mid-time where it was measured
    tcs = {x[0]: x[1] for x in tt}
    ep = np.round((t - t0) / P).astype(int)
    tau = t - np.array([tcs.get(k, t0 + k * P) for k in ep])
    cad = np.median(np.diff(np.sort(t)))
    resid = y - trapezoid_ss(tau, 0.0, dep, t14, t12, cad)
    h = t14 / 2 + cad
    sig = lcs.mad_std(resid[np.abs(tau) > h])
    tau_max = min(np.max(np.abs(tau)) - h, HILL_STABLE * r_hill / v_p / 86400)
    tau_min = max(2 * rp * REARTH / v_p / 86400, 0.25 * h)
    if np.unique(ep).size >= MIN_TRANSITS and tau_max > tau_min and np.isfinite(sig) and sig > 0:
        taus = np.geomspace(tau_min, tau_max, 16)
        dd, ee, snr, k, fap = ose_stack(tau, resid, ep, h, taus, sig, rng)
        d95 = np.maximum(np.nan_to_num(dd), 0) + 1.645 * ee
        rs = rstar * RSUN * np.sqrt(d95) / REARTH
        out['ose'] = {'tau_s_h': (taus * 24).tolist(), 'depth_ppm': (dd * 1e6).tolist(), 'snr': snr.tolist(),
                      'best_tau_h': float(taus[k] * 24), 'best_depth_ppm': float(dd[k] * 1e6), 'best_snr': float(snr[k]),
                      'fap': fap, 'rs_limit_earth': rs.tolist(), 'rs_limit_earth_any_orbit': float(np.nanmax(rs))}
        if fap < OSE_FAP and snr[k] >= OSE_SNR:
            out['signals'].append('ose_dip')
        out['_ose'] = (tau, resid, taus[k], h, dd[k])
    elif tau_max <= tau_min:
        out['notes'].append('no OSE stack: every stable moon orbit is crossed within the planetary transit itself')
    if out['signals']:
        out['status'] = 'moon_signal' if 'moon_timing' in out['signals'] or 'ose_dip' in out['signals'] else 'timing_variation'
    else:
        out['status'] = 'no_moon' if out['n_transits'] >= MIN_TRANSITS else 'too_few_transits'
    return out


# ---------------------------------------------------------------- reports
def report_md(o, name):
    L = [f"# {name} — moons: {o['status'].upper()}", '',
         f"Period {o['period_d']:.5f} d, {o['n_transits']} transits timed. Planet radius {o['rp_earth']:.2f} Re, mass "
         f"{o['mp_earth_assumed']:.1f} Me from the mass-radius relation, semi-major axis {o['a_p_au']:.4f} AU, Hill radius "
         f"{o['hill_radius_rp']:.1f} planet radii.", '']
    if 'ttv' in o:
        tv, dv, ph = o['ttv'], o['tdv'], o['phase']
        L += [f"TTV: rms {tv['rms_s'] / 60:.2f} min, chi2 {tv['chi2']:.1f} for {tv['dof']} degrees of freedom, highest "
              f"periodogram peak at {tv['freq_per_epoch']:.3f} cycles per transit with false-alarm probability {tv['fap']:.2g}; "
              f"95 % upper limit on an extra timing scatter {tv['limit95_s'] / 60:.2f} min.", '',
              f"TDV: rms {dv['rms_s'] / 60:.2f} min, chi2 {dv['chi2']:.1f} for {dv['dof']} degrees of freedom, false-alarm "
              f"probability {dv['fap']:.2g}; 95 % upper limit {dv['limit95_s'] / 60:.2f} min.", '',
              f"Phase of the TDV relative to the TTV at {ph['freq_per_epoch']:.3f} cycles per transit: {ph['dphi_deg']:.0f} "
              f"+- {ph['e_dphi_deg']:.0f} deg (a moon gives +-90 deg).", '']
        if o.get('mass_limit'):
            L += [f"Moon mass, 95 % upper limit: {o['mass_limit']['ms_earth_any_orbit']:.3g} Me for any orbit between 2 "
                  f"planet radii and {HILL_STABLE} Hill radii ({o['mass_limit']['ms_earth_widest_orbit']:.3g} Me on the widest).", '']
    if 'ose' in o:
        s = o['ose']
        L += [f"Orbital-sampling-effect stack: strongest extra dip {s['best_depth_ppm']:.0f} ppm at S/N {s['best_snr']:.1f} "
              f"for a moon orbit crossed in {s['best_tau_h']:.1f} h, false-alarm probability {s['fap']:.2g}; moon radius, "
              f"95 % upper limit {s['rs_limit_earth_any_orbit']:.2f} Re.", '']
    if o['notes']:
        L += ['Notes: ' + '; '.join(o['notes']), '']
    L += ['![moons](moons.png)', '']
    return '\n'.join(L)


def plot(o, path, name):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 1, figsize=(7, 8))
    if '_series' in o:
        n, oc, etc, tdv, ed = o['_series']
        axs[0].errorbar(n, oc * 1440, etc * 1440, fmt='o', ms=3, color='0.15')
        axs[1].errorbar(n, tdv * 1440, ed * 1440, fmt='o', ms=3, color='0.15')
    axs[0].set_ylabel('O-C (min)'); axs[1].set_ylabel('duration - mean (min)'); axs[1].set_xlabel('epoch')
    if '_ose' in o:
        tau, r, ts, h, d = o['_ose']
        b = np.linspace(tau.min(), tau.max(), 120)
        i = np.digitize(tau, b)
        m = np.array([np.mean(r[i == k]) if np.any(i == k) else np.nan for k in range(1, b.size)])
        axs[2].plot(tau * 24, r * 1e6, ',', color='0.6')
        axs[2].plot((b[:-1] + b[1:]) / 2 * 24, m * 1e6, 'o', ms=3, color='0.15')
        x = np.linspace(tau.min(), tau.max(), 600)
        axs[2].plot(x * 24, -np.nan_to_num(d) * ose_template(x, ts, h) * 1e6, color='C3')
        axs[2].axvspan(-h * 24, h * 24, color='C0', alpha=0.12)
    axs[2].set_xlabel('hours from mid-transit'); axs[2].set_ylabel('residual (ppm)')
    axs[0].set_title(f"{name}  P = {o['period_d']:.5f} d  {o['status']}", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def run_one(con, d, rng=None):
    """Moon search of the object whose vetting is in folder d; writes moons.md, moons.json, moons.png; returns the status."""
    s = json.loads((d / 'summary.json').read_text())
    z = np.load(d / 'transit_windows.npz')
    o = analyse(np.asarray(z['t'], float), np.asarray(z['f'], float), s, rng)
    name = d.name.replace('TIC', 'TIC ') if d.name.startswith('TIC') else d.name
    plot(o, d / 'moons.png', name)
    (d / 'moons.md').write_text(report_md(o, name), encoding='utf-8')
    (d / 'moons.json').write_text(json.dumps({k: v for k, v in o.items() if not k.startswith('_')}, indent=1, default=float))
    if con is not None and s.get('tic') is not None:
        con.execute('INSERT OR REPLACE INTO moon VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                    (s['tic'], o['status'], o['n_transits'], o.get('ttv', {}).get('rms_s'), o.get('ttv', {}).get('fap'),
                     o.get('tdv', {}).get('fap'), o.get('phase', {}).get('dphi_deg'), o.get('ose', {}).get('best_snr'),
                     o.get('ose', {}).get('fap'), o.get('mass_limit', {}).get('ms_earth_any_orbit'),
                     o.get('ose', {}).get('rs_limit_earth_any_orbit'), time.strftime('%Y-%m-%d %H:%M')))
        con.commit()
    return o['status']


def folders(findings, tics=None, include_rejected=False):
    """Vetting folders with a summary and transit windows, of planets that are ready (or of every transiting object)."""
    out = []
    for d in sorted(Path(findings).glob('*/')):
        if not (d / 'summary.json').exists() or not (d / 'transit_windows.npz').exists():
            continue
        s = json.loads((d / 'summary.json').read_text())
        if tics and s.get('tic') not in tics:
            continue
        if s.get('status') == 'ready' or include_rejected or tics:
            out.append(d)
    return out


def run_survey(con, findings):
    """Called by the survey: every ready planet without a moon search yet."""
    done = {r[0] for r in con.execute('SELECT tic FROM moon')}
    for d in folders(findings):
        tic = json.loads((d / 'summary.json').read_text()).get('tic')
        if tic not in done:
            run_one(con, d)


def main(argv=None):
    from . import paths, survey
    ap = argparse.ArgumentParser(prog='occulta moons', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='sub', required=True)
    s = sp.add_parser('search', help='TTV, TDV and OSE search of vetted transiting objects')
    s.add_argument('tic', nargs='*', type=int)
    s.add_argument('--include-rejected', action='store_true', help='also the transiting companions that the vetting rejected')
    a = ap.parse_args(argv)
    paths.init()
    con = survey.db()
    ds = folders(survey.FIND, set(a.tic) or None, a.include_rejected)
    counts = {}
    for d in ds:
        try:
            st = run_one(con, d)
        except Exception as e:
            st = 'error'
            survey.log('moons failed', d.name, repr(e)[:200])
        counts[st] = counts.get(st, 0) + 1
        print(d.name, st)
    print(f'{len(ds)} objects searched:', ', '.join(f'{k} {v}' for k, v in counts.items()) or 'none')
