"""Writes findings/status.html: progress of the search plan, candidates and recent problems."""
import html, sqlite3, time
from pathlib import Path

from . import paths

from . import sources as sources

ROOT = paths.HOME


def write(db_path=ROOT / 'data' / 'hunt.db', out=ROOT / 'findings' / 'status.html'):
    con = sqlite3.connect(db_path, timeout=60)
    q = lambda sql, *a: con.execute(sql, a).fetchall()
    jobs = q('SELECT job, status, n_targets, started, finished FROM job ORDER BY started')
    rows = []
    for jid, st, n, started, fin in jobs:
        u = int(jid.split(':')[-1])
        done = q('SELECT COUNT(*) FROM target WHERE sector=? AND status>0', u)[0][0]
        pct = 100 * done / n if n else 0
        rows.append(f'<tr><td>{html.escape(sources.unit_name(u))}</td><td>{html.escape(jid.split(":")[0])}'
                    f'{"" if jid.split(":")[1] == "None" else " T&le;" + jid.split(":")[1]}</td>'
                    f'<td class=num>{done:,} / {n:,}</td><td><div class=bar><span style="width:{pct:.1f}%"></span></div></td>'
                    f'<td>{html.escape(st)}</td><td>{fin or ""}</td></tr>')
    searched = q('SELECT COUNT(*) FROM target WHERE status=1')[0][0]
    events = q('SELECT COUNT(*) FROM event')[0][0]
    cstat = dict(q('SELECT status, COUNT(*) FROM candidate GROUP BY status'))
    good = q("SELECT tic, status, period, depth_ppm, dur_h, snr, notes FROM candidate "
             "WHERE status IN ('ready','review_user') ORDER BY snr DESC")
    crows = ''.join(
        f'<tr><td><a href="{sources.name_of(t).replace(" ", "")}/report.md">{sources.name_of(t)}</a></td><td>{s}</td>'
        f'<td class=num>{p:.4f}</td><td class=num>{d:.0f}</td><td class=num>{h:.2f}</td><td class=num>{x:.1f}</td>'
        f'<td>{html.escape((n or "")[:160])}</td></tr>' for t, s, p, d, h, x, n in good) or \
        '<tr><td colspan=7>No candidate has passed vetting yet.</td></tr>'
    log = (ROOT / 'logs' / 'hunt.log')
    tail = log.read_text(encoding='utf-8', errors='replace').splitlines()[-4000:] if log.exists() else []
    errs = [l for l in tail if 'vet failed' in l or 'store failed' in l][-8:]
    last = next((l for l in reversed(tail) if 'processed this session' in l), '')
    disk = sum(f.stat().st_size for f in ROOT.rglob('*') if f.is_file()) / 1e9
    page = f"""<!doctype html><html><head><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>Transit hunt status</title><style>
:root{{--bg:#fbfbf9;--fg:#1d1d1b;--mut:#6b6b66;--line:#e2e1dc;--acc:#2f6f8f}}
@media (prefers-color-scheme:dark){{:root{{--bg:#151614;--fg:#e8e7e2;--mut:#9a9993;--line:#2e2f2c;--acc:#7fb6d1}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif;margin:0 auto;max-width:1100px;padding:24px 16px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:16px;margin:28px 0 8px}} .mut{{color:var(--mut)}}
.kpis{{display:flex;flex-wrap:wrap;gap:12px;margin:16px 0}} .kpi{{border:1px solid var(--line);border-radius:8px;padding:10px 14px;min-width:150px}}
.kpi b{{display:block;font-size:20px}} table{{border-collapse:collapse;width:100%;font-size:14px}}
td,th{{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}} .num{{text-align:right;font-variant-numeric:tabular-nums}}
.bar{{background:var(--line);border-radius:4px;height:8px;min-width:120px}} .bar span{{display:block;height:8px;border-radius:4px;background:var(--acc)}}
a{{color:var(--acc)}} code{{font-size:13px}} .wrap{{overflow-x:auto}}
</style></head><body>
<h1>Transit hunt</h1><div class=mut>Updated {time.strftime('%Y-%m-%d %H:%M')} · TESS, K2 and Kepler · raw data purged after each search</div>
<div class=kpis><div class=kpi><b>{searched:,}</b>light curves searched</div><div class=kpi><b>{events:,}</b>dips recorded</div>
<div class=kpi><b>{cstat.get('ready', 0) + cstat.get('review_user', 0)}</b>candidates passing vetting</div>
<div class=kpi><b>{cstat.get('rejected', 0) + cstat.get('rejected_review', 0):,}</b>rejected false positives</div>
<div class=kpi><b>{cstat.get('new', 0)}</b>awaiting vetting</div><div class=kpi><b>{disk:.2f} GB</b>disk used of 30 GB</div></div>
<h2>Candidates passing vetting</h2><div class=wrap><table><tr><th>Star</th><th>Status</th><th class=num>Period (d)</th>
<th class=num>Depth (ppm)</th><th class=num>Duration (h)</th><th class=num>SNR</th><th>Notes</th></tr>{crows}</table></div>
<h2>Search plan</h2><div class=wrap><table><tr><th>Unit</th><th>Pass</th><th class=num>Searched</th><th>Progress</th><th>Status</th><th>Finished</th></tr>
{''.join(rows)}</table></div>
<h2>Recent problems</h2><div class=mut>{'<br>'.join(html.escape(e[:200]) for e in errs) or 'None.'}</div>
<p class=mut>Last progress line: <code>{html.escape(last)}</code></p></body></html>"""
    out.write_text(page, encoding='utf-8')
    return out


if __name__ == '__main__':
    print(write())
