"""Progress map of the search: (top) the plan as a grid of TESS sectors x plan steps, coloured by state; (bottom) the sky
positions of every star searched so far, coloured by the state of its light curve, on a sky shaded by the number of TESS
sectors per position when tess-point is installed (pip install tess-point). Usage: occulta survey map"""
import json, sqlite3
from pathlib import Path

from . import paths

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ROOT = paths.HOME
NS = 104
# light curves in MAST on 25 Oct 2026, millions
MAST_TOTALS = {'QLP': 107.6, 'SPOC 2-min': 1.84, 'K2SFF': 0.45, 'Kepler': 0.21}
K2 = [(98.30, 21.59), (173.94, 1.42), (246.12, -22.44), (336.67, -11.10), (59.08, 18.66), (130.16, 16.83), (204.86, -11.29),
      (287.83, -23.36), (16.31, 5.25), (270.35, -21.78), (186.78, -4.03), (260.32, -23.98), (351.66, -5.10), (72.79, 20.79),
      (160.68, 6.85), (233.61, -20.08), (133.65, 18.53), (202.50, -7.72), (130.69, 16.83), (347.21, -4.21)]


def footprint(n=20000):
    """Number of TESS sectors (1-NS) at n random sky positions, as (ra, dec, nsec) arrays, or None without tess-point. The
    first call takes about 7 min and keeps the result in <home>/data/footprint.npz."""
    fp = ROOT / 'data' / 'footprint.npz'
    if not fp.exists():
        try:
            import tess_stars2px as t2p
        except ImportError:
            return None
        print(f'building the TESS footprint at {n} positions (about 7 min, once) ...', flush=True)
        rng = np.random.default_rng(1)
        ra = rng.uniform(0, 360, n); dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
        res = t2p.tess_stars2px_function_entry(np.arange(n), ra, dec)
        ok = res[3] <= NS
        fp.parent.mkdir(parents=True, exist_ok=True)
        np.savez(fp, ra=ra, dec=dec, nsec=np.bincount(res[0][ok], minlength=n))
    f = np.load(fp)
    return f['ra'], f['dec'], f['nsec']


def circle(ax, ra0, dec0, r, **kw):
    """Outline of a field of radius r (deg) on the Mollweide axes, split where it crosses RA = 180 deg."""
    a = np.linspace(0, 2 * np.pi, 120)
    r_, d0 = np.radians(r), np.radians(dec0)
    d = np.degrees(np.arcsin(np.sin(d0) * np.cos(r_) + np.cos(d0) * np.sin(r_) * np.cos(a)))
    rr = ra0 + np.degrees(np.arctan2(np.sin(a) * np.sin(r_) * np.cos(d0), np.cos(r_) - np.sin(d0) * np.sin(np.radians(d))))
    x = -np.radians(((rr + 180) % 360) - 180); y = np.radians(d)
    cut = np.where(np.abs(np.diff(x)) > np.pi)[0] + 1
    for k, (xs, ys) in enumerate(zip(np.split(x, cut), np.split(y, cut))):
        ax.plot(xs, ys, **(kw if k == 0 else {**kw, 'label': None}))


def main(argv=None):
    from . import survey
    if not (ROOT / 'data' / 'hunt.db').exists():      # a new working directory: draw the map of an empty search
        (ROOT / 'data').mkdir(parents=True, exist_ok=True)
        survey.db(ROOT / 'data' / 'hunt.db').close()
    con = sqlite3.connect(f"file:{ROOT / 'data' / 'hunt.db'}?mode=ro", uri=True, timeout=60)
    plan = paths.config()['plan']
    steps = [f"{p['mission']}" + (f" T<={p['tmag']:g}" if 'tmag' in p else (' all' if p['mission'] == 'QLP' else '')) for p in plan]
    # 0 not applicable, 1 planned, 2 queued (pending), 3 running, 4 done
    grid = np.zeros((len(steps), NS), int)
    for i, p in enumerate(plan):
        if p['mission'] in ('QLP', 'SPOC'):
            grid[i, :] = 1
        elif p['mission'] == 'K2':
            grid[i, :] = 0
        elif p['mission'] == 'Kepler':
            grid[i, :] = 0
    counts = {int(s): (n, d, q) for s, n, d, q in con.execute(
        'SELECT sector, COUNT(*), SUM(status=1), SUM(status=0) FROM target GROUP BY sector') if str(s).isdigit()}
    for s, (n, d, q) in counts.items():
        if 1 <= s <= NS:
            grid[0, s - 1] = 4 if q == 0 else 3
    done_lc = sum(d for n, d, q in counts.values()); total_q = sum(n for n, d, q in counts.values())

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_axes([0.11, 0.62, 0.86, 0.33])
    cmap = ListedColormap(['#ffffff', '#d9d9d9', '#f4c542', '#e8743b', '#2e7d32'])
    ax.imshow(grid, aspect='auto', cmap=cmap, vmin=0, vmax=4, interpolation='nearest')
    ax.set_yticks(range(len(steps))); ax.set_yticklabels([f'{i + 1}. {s}' for i, s in enumerate(steps)])
    ax.set_xticks(np.arange(0, NS, 4)); ax.set_xticklabels(np.arange(1, NS + 1, 4)); ax.set_xlabel('TESS sector')
    ax.set_xticks(np.arange(-0.5, NS, 1), minor=True); ax.set_yticks(np.arange(-0.5, len(steps), 1), minor=True)
    ax.grid(which='minor', color='0.6', lw=0.3); ax.tick_params(which='minor', length=0)
    for i, p in enumerate(plan):
        if p['mission'] in ('K2', 'Kepler'):
            ax.text(NS / 2, i, f"{p['mission']} (not sector-based): planned", ha='center', va='center', fontsize=9, color='0.3')
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cmap.colors[1:], ['planned', 'queued', 'running', 'done'])],
              loc='upper left', bbox_to_anchor=(0, 1.18), ncol=4, frameon=False)
    ax.set_title(f'Search plan by sector (newest first). Light curves searched: {done_lc:,} of {total_q:,} queued', loc='right', fontsize=10)

    ax2 = fig.add_axes([0.08, 0.09, 0.84, 0.46], projection='mollweide')
    fp = footprint()
    if fp is None:
        print('tess-point is not installed: map without the TESS sector shading (pip install tess-point)')
    else:
        ra_f = np.radians(((fp[0] + 180) % 360) - 180)
        sc = ax2.scatter(-ra_f, np.radians(fp[1]), s=7, c=np.clip(fp[2], 0, 20), cmap='Greys', vmin=-4, vmax=20,
                         rasterized=True, zorder=1)
        cb = fig.colorbar(sc, ax=ax2, fraction=0.025, pad=0.02, ticks=[0, 5, 10, 15, 20])
        cb.set_label(f'TESS sectors available (1-{NS}); 20 = 20 or more', fontsize=8)
        print('TESS covers %.1f %% of the sky, mean %.1f sectors per position' % (100 * np.mean(fp[2] > 0), np.mean(fp[2])))
    for i, (r0, d0) in enumerate(K2):
        circle(ax2, r0, d0, 6.5, color='#1f5fa8', lw=0.8, zorder=3, label='K2 campaigns' if i == 0 else None)
    circle(ax2, 290.67, 44.5, 7.5, color='#8e24aa', lw=1.2, zorder=3, label='Kepler field')
    print('light curves in MAST (25 Oct 2026): ' + ', '.join(f'{k} {v:g} M' for k, v in MAST_TOTALS.items()))
    rows = np.array(con.execute('SELECT t.ra, t.dec, MIN(g.status) FROM target g JOIN tic t ON t.tic = g.tic '
                                'WHERE t.ra IS NOT NULL GROUP BY g.tic').fetchall(), float)
    if rows.size:
        ra = np.radians(((rows[:, 0] + 180) % 360) - 180); dec = np.radians(rows[:, 1])
        for st, c, l in ((1, '#2e7d32', 'searched'), (0, '#f4c542', 'queued'), (2, '#e8743b', 'failed / retry'), (3, '#9e9e9e', 'skipped')):
            m = rows[:, 2] == st
            ax2.scatter(-ra[m], dec[m], s=0.3, c=c, label=f'{l} ({m.sum():,})', rasterized=True, zorder=2)
    ax2.grid(True, lw=0.3); ax2.set_xticklabels([])
    ax2.legend(loc='lower left', bbox_to_anchor=(0, -0.12), ncol=3, markerscale=15, frameon=False)
    ax2.set_title('Stars searched so far (equatorial, Mollweide); MAST light curves, 25 Oct 2026: ' +
                  ', '.join(f'{k} {v:g} M' for k, v in MAST_TOTALS.items()), fontsize=9)
    out = ROOT / 'findings' / 'progress' / 'progress_map.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110); print(out, done_lc, total_q)



if __name__ == '__main__':
    main()
