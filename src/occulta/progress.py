"""Progress map of the search: (top) the plan as a grid of TESS sectors x plan steps, coloured by state; (bottom) the sky
positions of every star searched so far, coloured by the state of its light curve. Usage: occulta survey map"""
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
def main(argv=None):
    con = sqlite3.connect(f"file:{ROOT / 'data' / 'hunt.db'}?mode=ro", uri=True, timeout=60)
    plan = paths.config()['plan']
    steps = [f"{p['mission']}" + (f" T<={p['tmag']:g}" if 'tmag' in p else (' all' if p['mission'] == 'QLP' else '')) for p in plan]
    NS = 104
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

    ax2 = fig.add_axes([0.08, 0.05, 0.84, 0.5], projection='mollweide')
    rows = np.array(con.execute('SELECT t.ra, t.dec, MIN(g.status) FROM target g JOIN tic t ON t.tic = g.tic '
                                'WHERE t.ra IS NOT NULL GROUP BY g.tic').fetchall(), float)
    if rows.size:
        ra = np.radians(((rows[:, 0] + 180) % 360) - 180); dec = np.radians(rows[:, 1])
        for st, c, l in ((1, '#2e7d32', 'searched'), (0, '#f4c542', 'queued'), (2, '#e8743b', 'failed / retry'), (3, '#9e9e9e', 'skipped')):
            m = rows[:, 2] == st
            ax2.scatter(-ra[m], dec[m], s=0.3, c=c, label=f'{l} ({m.sum():,})', rasterized=True)
    ax2.grid(True, lw=0.3); ax2.set_xticklabels([])
    ax2.legend(loc='lower left', bbox_to_anchor=(0, -0.12), ncol=4, markerscale=15, frameon=False)
    ax2.set_title('Stars searched so far (equatorial, Mollweide)', fontsize=10)
    out = ROOT / 'findings' / 'progress' / 'progress_map.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110); print(out, done_lc, total_q)



if __name__ == '__main__':
    main()
