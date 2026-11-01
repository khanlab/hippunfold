"""Aggregate eval_template_tuning.py outputs into a tsv and a self-contained html report.

Can be run from snakemake (rule tune_atlas_report) or from the command line:
  python make_tuning_report.py --tsvs a.tsv b.tsv ... --pngs a.png b.png ... \
      --out_tsv all.tsv --out_html report.html [--param_keys run] [--title ...]
"""

import argparse
import base64
import html
import json
from pathlib import Path

import pandas as pd

# metrics shown in the report: column, label, higher is better
METRICS = [
    ("loo_dice_mean", "LOO subfield Dice (mean)", True),
    ("loo_dice_min", "LOO subfield Dice (worst input)", True),
    ("loo_corr_mean", "LOO metric correlation (mean)", True),
    ("jac_logjac_std", "std log det(J)", False),
    ("jac_folding_frac", "folding fraction (max)", False),
    ("jac_disp_mean_mm", "mean displacement (mm)", None),
]


def load(tsvs, pngs):
    df = pd.concat(
        [pd.read_csv(f, sep="\t", dtype=str) for f in tsvs], ignore_index=True
    )
    for col in df.columns:
        converted = pd.to_numeric(df[col], errors="coerce")
        if (
            col.startswith(("loo_", "jac_", "sharpness_", "n_"))
            and converted.notna().sum() == df[col].notna().sum()
        ):
            df[col] = converted
    df["png"] = [
        base64.b64encode(Path(p).read_bytes()).decode() if p else "" for p in pngs
    ]
    return df


def run_name(row, keys):
    if keys == ["run"]:
        return str(row["run"])
    return " ".join(f"{k}-{row[k]}" for k in keys)


def make_html(df, param_keys, base_params, title):
    varying = [k for k in param_keys if k in df and df[k].nunique() > 1]
    # default modality weights just follow the number of metrics
    if (
        "wt" in varying
        and "mets" in df
        and df.groupby("mets")["wt"].nunique().max() == 1
    ):
        varying.remove("wt")
    fixed = {k: df[k].iloc[0] for k in param_keys if k in df and k not in varying}
    name_keys = varying or [k for k in param_keys if k in df]
    df = df.copy()
    df["name"] = [run_name(r, name_keys) for _, r in df.iterrows()]

    sort_col = "loo_dice_mean" if "loo_dice_mean" in df else "loo_corr_mean"
    df = df.sort_values(sort_col, ascending=False).reset_index(drop=True)

    is_base = [
        bool(base_params)
        and all(str(r[k]) == str(v) for k, v in base_params.items() if k in df)
        for _, r in df.iterrows()
    ]
    df["is_base"] = is_base

    metric_cols = [m for m in METRICS if m[0] in df]
    records = []
    for _, r in df.iterrows():
        rec = {"name": r["name"], "is_base": bool(r["is_base"])}
        rec.update({k: str(r[k]) for k in varying})
        rec.update(
            {c: (None if pd.isna(r[c]) else float(r[c])) for c, _, _ in metric_cols}
        )
        records.append(rec)

    best = df.iloc[0]
    tiles = [
        (
            "best " + {c: l for c, l, _ in metric_cols}[sort_col],
            f"{best[sort_col]:.3f}",
            best["name"],
        )
    ]
    if "loo_dice_unwarped_mean" in df:
        tiles.append(
            (
                "unwarped LOO Dice (no registration)",
                f"{df['loo_dice_unwarped_mean'].iloc[0]:.3f}",
                "same for all runs",
            )
        )
    if df["is_base"].any():
        b = df[df["is_base"]].iloc[0]
        rank = int(df.index[df["is_base"]][0]) + 1
        tiles.append(
            (
                "current config (new_atlas_ants)",
                f"{b[sort_col]:.3f}",
                f"rank {rank} of {len(df)}",
            )
        )
    tiles.append(("runs", str(len(df)), f"{len(varying)} varied parameter(s)"))

    tiles_html = "".join(
        f'<div class="tile"><div class="tlabel">{html.escape(a)}</div>'
        f'<div class="tval">{html.escape(b)}</div><div class="tsub">{html.escape(c)}</div></div>'
        for a, b, c in tiles
    )

    fixed_html = (
        ", ".join(
            f"<code>{html.escape(k)}={html.escape(str(v))}</code>"
            for k, v in fixed.items()
        )
        or "none"
    )

    head_cells = "".join(
        f"<th data-k='{k}'>{html.escape(k)}</th>" for k in varying
    ) + "".join(
        f"<th data-k='{c}' class='num'>{html.escape(l)}</th>" for c, l, _ in metric_cols
    )
    table_rows = []
    for i, r in df.iterrows():
        cells = "".join(f"<td>{html.escape(str(r[k]))}</td>" for k in varying)
        for c, _, _ in metric_cols:
            v = r[c]
            fmt = (
                f"{v:.2e}"
                if c == "jac_folding_frac" and v and v < 0.001
                else f"{v:.4f}"
            )
            cells += f"<td class='num'>{'' if pd.isna(v) else fmt}</td>"
        cls = " class='base'" if r["is_base"] else ""
        table_rows.append(f"<tr{cls}><td class='num'>{i + 1}</td>{cells}</tr>")

    thumbs = "".join(
        f"<figure><img alt='{html.escape(r['name'])}' src='data:image/png;base64,{r['png']}'>"
        f"<figcaption>#{i + 1} {html.escape(r['name'])} &middot; {sort_col}={r[sort_col]:.3f}</figcaption></figure>"
        for i, r in df.head(6).iterrows()
        if r["png"]
    )

    # strip plots only for parameters shared by more than one run
    grouped = [k for k in varying if df[k].nunique() < len(df)]
    data_json = json.dumps(
        {
            "records": records,
            "varying": grouped,
            "metrics": [[c, l, h] for c, l, h in metric_cols],
        }
    )

    return TEMPLATE.format(
        title=html.escape(title),
        tiles=tiles_html,
        fixed=fixed_html,
        sort_label=html.escape({c: l for c, l, _ in metric_cols}[sort_col]),
        head=head_cells,
        rows="".join(table_rows),
        thumbs=thumbs,
        data=data_json.replace("</", "<\\/"),
        sort_col=sort_col,
    )


TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
.viz-root {{
  color-scheme: light;
  --surface-0: #f5f5f3; --surface-1: #fcfcfb; --border: #e2e1dc;
  --text-primary: #0b0b0b; --text-secondary: #52514e; --text-muted: #7a7974;
  --grid: #ecebe7; --series-1: #2a78d6; --series-2: #eb6834;
}}
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .viz-root {{
    color-scheme: dark;
    --surface-0: #121211; --surface-1: #1a1a19; --border: #33332f;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #8f8e86;
    --grid: #2a2a27; --series-1: #3987e5; --series-2: #d95926;
  }}
}}
:root[data-theme="dark"] .viz-root {{
  color-scheme: dark;
  --surface-0: #121211; --surface-1: #1a1a19; --border: #33332f;
  --text-primary: #f0efec; --text-secondary: #c3c2b7; --text-muted: #8f8e86;
  --grid: #2a2a27; --series-1: #3987e5; --series-2: #d95926;
}}
body {{ margin: 0; }}
.viz-root {{ background: var(--surface-0); color: var(--text-primary);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; min-height: 100vh; }}
main {{ max-width: 1180px; margin: 0 auto; padding: 28px 20px 60px; }}
h1 {{ font-size: 22px; margin: 0 0 4px; }}
h2 {{ font-size: 16px; margin: 32px 0 6px; }}
p.note, .sub {{ color: var(--text-secondary); margin: 4px 0 12px; max-width: 900px; }}
code {{ font-size: 12px; background: var(--surface-1); border: 1px solid var(--border); border-radius: 4px; padding: 0 4px; }}
.tiles {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 12px; margin: 18px 0; }}
.tile {{ background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 12px 14px; }}
.tlabel {{ color: var(--text-secondary); font-size: 12px; }}
.tval {{ font-size: 26px; font-weight: 600; font-variant-numeric: tabular-nums; }}
.tsub {{ color: var(--text-muted); font-size: 12px; overflow-wrap: anywhere; }}
.card {{ background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 14px; }}
.controls {{ display: flex; gap: 14px; align-items: center; flex-wrap: wrap; margin-bottom: 8px; color: var(--text-secondary); font-size: 13px; }}
select {{ font: inherit; padding: 2px 6px; background: var(--surface-1); color: var(--text-primary); border: 1px solid var(--border); border-radius: 6px; }}
.legend {{ display: flex; gap: 16px; font-size: 12px; color: var(--text-secondary); }}
.sw {{ display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 5px; vertical-align: -1px; }}
svg text {{ fill: var(--text-secondary); font-size: 11px; }}
svg .axis line, svg .axis path {{ stroke: var(--border); }}
svg .gridl {{ stroke: var(--grid); }}
.small {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 12px; }}
.tip {{ position: fixed; pointer-events: none; background: var(--surface-1); color: var(--text-primary);
  border: 1px solid var(--border); border-radius: 8px; padding: 8px 10px; font-size: 12px;
  box-shadow: 0 4px 14px rgba(0,0,0,.15); display: none; max-width: 320px; z-index: 10; }}
.tip b {{ font-weight: 600; }}
.tip .tv {{ font-variant-numeric: tabular-nums; }}
.tablewrap {{ overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; }}
th, td {{ padding: 5px 8px; border-bottom: 1px solid var(--border); text-align: left; white-space: nowrap; }}
th {{ cursor: pointer; color: var(--text-secondary); font-weight: 600; position: sticky; top: 0; background: var(--surface-1); }}
td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
tr.base td {{ box-shadow: inset 3px 0 0 var(--series-2); }}
.thumbs {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 12px; }}
figure {{ margin: 0; background: #fff; border: 1px solid var(--border); border-radius: 10px; padding: 6px; }}
figure img {{ width: 100%; display: block; }}
figcaption {{ font-size: 11px; color: #52514e; padding: 4px 2px 0; overflow-wrap: anywhere; }}
</style></head>
<body><div class="viz-root"><main>
<h1>{title}</h1>
<p class="sub">Each run is one <code>antsMultivariateTemplateConstruction2.sh</code> template build.
Runs are ranked by {sort_label}: each input's subfield labels are warped to the template and compared with the
majority vote of the other inputs. The subfields are not used to drive the registration, so this measures
anatomical correspondence independently. Regularity is std of log det(J) of the input&rarr;template warps (lower is smoother).
Parameters fixed across runs: {fixed}.
<code>mets</code> is the set of metrics used to build the template (curv, gyr = gyrification, lgyr = legacy gyrification, thick).
The metric correlation averages over whichever metrics a run used, so compare it only between runs with the same <code>mets</code>.</p>
<div class="tiles">{tiles}</div>

<h2>Accuracy vs. regularity</h2>
<p class="note">Upper-left is better: higher agreement with smoother warps. A run with any folding (det(J) &le; 0) is drawn hollow.</p>
<div class="card">
  <div class="controls">
    <label>y: <select id="ysel"></select></label>
    <label>x: <select id="xsel"></select></label>
    <span class="legend"><span><span class="sw" style="background:var(--series-1)"></span>tuning run</span>
    <span id="baseleg"><span class="sw" style="background:var(--series-2)"></span>current config</span></span>
  </div>
  <svg id="scatter" width="100%" height="380" role="img" aria-label="scatter of runs"></svg>
</div>

<h2>Effect of each varied parameter</h2>
<p class="note">Every run's score (y metric above), grouped by one parameter's value; the bar marks the mean over the other parameters.</p>
<div class="small" id="small"></div>

<h2>All runs</h2>
<p class="note">Click a column to sort. The current config row is marked on the left.</p>
<div class="card tablewrap"><table id="tbl"><thead><tr><th data-k="rank" class="num">#</th>{head}</tr></thead><tbody>{rows}</tbody></table></div>

<h2>Top runs</h2>
<p class="note">Template images and the majority vote of the warped subfields.</p>
<div class="thumbs">{thumbs}</div>
</main></div>
<div class="tip" id="tip"></div>
<script>
const D = {data};
const SORT = "{sort_col}";
const tip = document.getElementById("tip");
const NS = "http://www.w3.org/2000/svg";
const metrics = D.metrics.filter(m => D.records.some(r => r[m[0]] != null));
const label = k => (metrics.find(m => m[0] === k) || [k, k])[1];
const fmt = v => v == null ? "" : (Math.abs(v) < 0.001 && v !== 0 ? v.toExponential(2) : v.toFixed(4));
function el(tag, attrs, parent) {{ const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }}
function showTip(ev, r, extra) {{
  let h = `<b>${{r.name}}</b>`;
  for (const k of extra) h += `<div>${{label(k)}}: <span class="tv">${{fmt(r[k])}}</span></div>`;
  if (r.is_base) h += "<div>current config</div>";
  tip.innerHTML = h; tip.style.display = "block";
  const x = Math.min(ev.clientX + 14, innerWidth - tip.offsetWidth - 8);
  tip.style.left = x + "px"; tip.style.top = (ev.clientY + 14) + "px"; }}
function hideTip() {{ tip.style.display = "none"; }}
function ticks(lo, hi, n) {{ const s = (hi - lo) / n; const p = Math.pow(10, Math.floor(Math.log10(s)));
  const st = [1, 2, 5, 10].map(m => m * p).find(m => m >= s); const out = [];
  for (let v = Math.ceil(lo / st) * st; v <= hi + 1e-12; v += st) out.push(+v.toFixed(10)); return out; }}
function extent(vals) {{ let lo = Math.min(...vals), hi = Math.max(...vals);
  if (lo === hi) {{ lo -= 0.01 * Math.abs(lo || 1); hi += 0.01 * Math.abs(hi || 1); }}
  const pad = (hi - lo) * 0.08; return [lo - pad, hi + pad]; }}
function axes(svg, W, H, m, xs, ys, xd, yd, xl, yl) {{
  const g = el("g", {{class: "axis"}}, svg);
  for (const t of ticks(yd[0], yd[1], 5)) {{ el("line", {{x1: m.l, x2: W - m.r, y1: ys(t), y2: ys(t), class: "gridl"}}, g);
    el("text", {{x: m.l - 6, y: ys(t) + 4, "text-anchor": "end"}}, g).textContent = fmt(t).replace(/0+$/, "").replace(/\\.$/, ""); }}
  if (xd) for (const t of ticks(xd[0], xd[1], 5)) {{ el("text", {{x: xs(t), y: H - m.b + 16, "text-anchor": "middle"}}, g).textContent = fmt(t).replace(/0+$/, "").replace(/\\.$/, ""); }}
  el("line", {{x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b}}, g);
  if (xl) el("text", {{x: (m.l + W - m.r) / 2, y: H - 6, "text-anchor": "middle"}}, g).textContent = xl;
  const t = el("text", {{x: 14, y: (m.t + H - m.b) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${{(m.t + H - m.b) / 2}})`}}, g); t.textContent = yl; }}
function dot(svg, cx, cy, r, hollow, base) {{
  const c = base ? "var(--series-2)" : "var(--series-1)";
  return el("circle", {{cx, cy, r: base ? 6 : 5, fill: hollow ? "var(--surface-1)" : c, stroke: hollow ? c : "var(--surface-1)", "stroke-width": 2}}, svg); }}
function scatter() {{
  const svg = document.getElementById("scatter"); svg.innerHTML = "";
  const yk = document.getElementById("ysel").value, xk = document.getElementById("xsel").value;
  const W = svg.clientWidth || 900, H = 380, m = {{l: 64, r: 16, t: 12, b: 42}};
  const R = D.records.filter(r => r[xk] != null && r[yk] != null);
  const xd = extent(R.map(r => r[xk])), yd = extent(R.map(r => r[yk]));
  const xs = v => m.l + (v - xd[0]) / (xd[1] - xd[0]) * (W - m.l - m.r);
  const ys = v => H - m.b - (v - yd[0]) / (yd[1] - yd[0]) * (H - m.t - m.b);
  axes(svg, W, H, m, xs, ys, xd, yd, label(xk), label(yk));
  const order = R.slice().sort((a, b) => a.is_base - b.is_base);
  for (const r of order) {{
    const hollow = (r.jac_folding_frac || 0) > 0;
    const d = dot(svg, xs(r[xk]), ys(r[yk]), 5, hollow, r.is_base);
    const hit = el("circle", {{cx: xs(r[xk]), cy: ys(r[yk]), r: 11, fill: "transparent"}}, svg);
    hit.addEventListener("mousemove", ev => {{ d.setAttribute("r", 7); showTip(ev, r, [yk, xk, "jac_folding_frac"]); }});
    hit.addEventListener("mouseleave", () => {{ d.setAttribute("r", r.is_base ? 6 : 5); hideTip(); }});
  }}
  const best = R.slice().sort((a, b) => b[yk] - a[yk])[0];
  if (best) {{ const right = xs(best[xk]) > (m.l + W - m.r) / 2;
    const t = el("text", {{x: xs(best[xk]) + (right ? -9 : 9), y: ys(best[yk]) - 9, "text-anchor": right ? "end" : "start"}}, svg);
    t.textContent = "best: " + best.name; t.style.fill = "var(--text-primary)"; }}
  small(yk);
}}
function small(yk) {{
  const box = document.getElementById("small"); box.innerHTML = "";
  if (!D.varying.length) {{ box.innerHTML = "<p class='note'>No parameter is shared by more than one run.</p>"; return; }}
  const yd = extent(D.records.filter(r => r[yk] != null).map(r => r[yk]));
  for (const k of D.varying) {{
    const card = document.createElement("div"); card.className = "card"; box.appendChild(card);
    if (new Set(D.records.map(r => r[k])).size > 5) card.style.gridColumn = "1 / -1";
    const svg = el("svg", {{width: "100%", height: 220}}, card);
    const W = card.clientWidth - 28 || 280, H = 220, m = {{l: 52, r: 10, t: 22, b: 30}};
    const vals = [...new Set(D.records.map(r => r[k]))].sort((a, b) => (+a - +b) || String(a).localeCompare(String(b)));
    const bw = (W - m.l - m.r) / vals.length;
    const xs = v => m.l + (vals.indexOf(v) + 0.5) * bw;
    const ys = v => H - m.b - (v - yd[0]) / (yd[1] - yd[0]) * (H - m.t - m.b);
    axes(svg, W, H, m, xs, ys, null, yd, null, "");
    el("text", {{x: m.l, y: 13}}, svg).textContent = k; svg.lastChild.style.fill = "var(--text-primary)";
    vals.forEach(v => {{
      el("text", {{x: xs(v), y: H - m.b + 16, "text-anchor": "middle"}}, svg).textContent = v;
      const rs = D.records.filter(r => r[k] === v && r[yk] != null);
      const mean = rs.reduce((s, r) => s + r[yk], 0) / rs.length;
      el("line", {{x1: xs(v) - bw * 0.3, x2: xs(v) + bw * 0.3, y1: ys(mean), y2: ys(mean), stroke: "var(--text-secondary)", "stroke-width": 2, "stroke-linecap": "round"}}, svg);
      rs.forEach((r, i) => {{
        const jx = xs(v) + ((i % 7) - 3) * Math.min(5, bw / 16);
        const d = dot(svg, jx, ys(r[yk]), 4, (r.jac_folding_frac || 0) > 0, r.is_base); d.setAttribute("r", r.is_base ? 5 : 4);
        const hit = el("circle", {{cx: jx, cy: ys(r[yk]), r: 8, fill: "transparent"}}, svg);
        hit.addEventListener("mousemove", ev => showTip(ev, r, [yk, "jac_logjac_std"]));
        hit.addEventListener("mouseleave", hideTip);
      }});
    }});
  }}
}}
const ysel = document.getElementById("ysel"), xsel = document.getElementById("xsel");
for (const [k, l, h] of metrics) {{
  if (h !== false) ysel.add(new Option(l, k, k === SORT, k === SORT));
  xsel.add(new Option(l, k, k === "jac_logjac_std", k === "jac_logjac_std")); }}
if (![...xsel.options].some(o => o.selected)) xsel.selectedIndex = 0;
if (!D.records.some(r => r.is_base)) document.getElementById("baseleg").remove();
ysel.onchange = xsel.onchange = scatter; addEventListener("resize", scatter); scatter();
// sortable table
const tbody = document.querySelector("#tbl tbody");
document.querySelectorAll("#tbl th").forEach((th, ci) => {{ let asc = false;
  th.onclick = () => {{ asc = !asc; const rows = [...tbody.rows];
    rows.sort((a, b) => {{ const x = a.cells[ci].textContent, y = b.cells[ci].textContent;
      const nx = parseFloat(x), ny = parseFloat(y);
      const c = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : x.localeCompare(y); return asc ? c : -c; }});
    rows.forEach(r => tbody.appendChild(r)); }}; }});
</script>
</body></html>
"""


def main(tsvs, pngs, out_tsv, out_html, param_keys, base_params, title):
    df = load(tsvs, pngs)
    df.drop(columns="png").to_csv(out_tsv, sep="\t", index=False)
    Path(out_html).write_text(make_html(df, param_keys, base_params, title))


if "snakemake" in globals():
    main(
        tsvs=snakemake.input.tsvs,
        pngs=snakemake.input.pngs,
        out_tsv=snakemake.output.tsv,
        out_html=snakemake.output.html,
        param_keys=["label", "hemi"] + list(snakemake.params.tune_keys),
        base_params=dict(snakemake.params.base_params),
        title=f"Template construction tuning: tpl-{snakemake.config['new_atlas_name']}",
    )
else:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsvs", nargs="+", required=True)
    parser.add_argument("--pngs", nargs="+", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--out_html", required=True)
    parser.add_argument("--param_keys", nargs="+", default=["run"])
    parser.add_argument("--title", default="Template construction tuning")
    args = parser.parse_args()
    main(
        args.tsvs,
        args.pngs,
        args.out_tsv,
        args.out_html,
        args.param_keys,
        {},
        args.title,
    )
