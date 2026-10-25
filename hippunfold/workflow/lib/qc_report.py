"""
Self-contained HTML QC reports.

A report is a single HTML file with no external assets: figures are embedded
as inline SVG (with colours mapped to CSS variables, so they follow the light
or dark theme) or as base64 PNG. Typical use:

    report = QCReport("Source/sink terminals", entities={"subject": "001"})
    report.add_summary(results)
    for result in results:
        section = report.add_section(result.label, status=result.status)
        section.add_checks(result)
        section.add_figure(fig, caption="...")
    report.save("report.html")
"""

import base64
import datetime
import html
import io
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Colours used when drawing figures (light theme), and the CSS variable each
# maps to in inline SVG, so figures follow the page theme.
INK = {
    "text_primary": "#0b0b0b",
    "text_secondary": "#52514e",
    "muted": "#898781",
    "axis": "#c3c2b7",
    "grid": "#e1e0d9",
    "surface": "#fcfcfb",
}
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # categorical slots 1-4
SERIES_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500"]

# sequential ramps (light -> dark) for continuous values on surfaces
BLUE_RAMP = ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORANGE_RAMP = ["#fde0d2", "#f6b293", "#eb6834", "#c4521f", "#963c14", "#6b2a0c"]
# diverging ramp (blue -> neutral -> orange) for signed values centred on 0
DIVERGING_RAMP = BLUE_RAMP[4:0:-1] + ["#e6e5e0"] + ORANGE_RAMP[1:5]

_SVG_COLOR_VARS = {
    INK["text_primary"]: "var(--text-primary)",
    INK["text_secondary"]: "var(--text-secondary)",
    INK["muted"]: "var(--text-muted)",
    INK["axis"]: "var(--axis)",
    INK["grid"]: "var(--grid)",
    INK["surface"]: "var(--surface-1)",
    **{c: f"var(--series-{i + 1})" for i, c in enumerate(SERIES)},
}

STATUS_LABEL = {"pass": "Pass", "warn": "Warning", "fail": "Fail"}
STATUS_ICON = {"pass": "✓", "warn": "!", "fail": "✕"}

MPL_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9,
    "text.color": INK["text_secondary"],
    "axes.labelcolor": INK["text_secondary"],
    "axes.edgecolor": INK["axis"],
    "axes.facecolor": "none",
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": INK["grid"],
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": INK["muted"],
    "ytick.color": INK["muted"],
    "xtick.labelcolor": INK["text_secondary"],
    "ytick.labelcolor": INK["text_secondary"],
    "lines.linewidth": 1.5,
    "figure.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
    "legend.frameon": False,
}


def style():
    """Context manager applying the report's matplotlib style."""
    return plt.rc_context(MPL_STYLE)


def fig_to_svg(fig):
    """Render a figure as inline SVG, with report colours mapped to CSS vars."""
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    svg = buf.getvalue()
    svg = svg[svg.index("<svg") :]
    # drop fixed size so it scales to its container
    svg = re.sub(r'(<svg[^>]*?)\s+width="[^"]*"\s+height="[^"]*"', r"\1", svg, count=1)
    svg = svg.replace("font-family:'DejaVu Sans'", "font-family:var(--font)")
    svg = re.sub(r"font-family: *'?DejaVu Sans'?", "font-family: var(--font)", svg)
    for hex_color, var in _SVG_COLOR_VARS.items():
        svg = re.sub(re.escape(hex_color), var, svg, flags=re.IGNORECASE)
    return svg


def fig_to_png(fig, dpi=110):
    """Render a figure as a base64 PNG data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def status_chip(status):
    return (
        f'<span class="chip chip-{status}"><span class="chip-icon" aria-hidden="true">'
        f"{STATUS_ICON[status]}</span>{STATUS_LABEL[status]}</span>"
    )


def legend(items):
    """HTML legend from (label, series index) pairs."""
    keys = "".join(
        f'<span class="key"><span class="swatch" style="background:var(--series-{i + 1})">'
        f"</span>{html.escape(label)}</span>"
        for label, i in items
    )
    return f'<div class="legend">{keys}</div>'


def hex_rgb(colors):
    """Hex colours as an (N, 3) array of RGB in [0, 1]."""
    return np.array([[int(c[i : i + 2], 16) / 255 for i in (1, 3, 5)] for c in colors])


def ramp_rgb(values, ramp, vmin=0.0, vmax=1.0):
    """Map values linearly onto a colour ramp, as an (N, 3) RGB array."""
    ramp = hex_rgb(ramp)
    t = np.clip((np.asarray(values, dtype=float) - vmin) / (vmax - vmin), 0, 1)
    pos = np.nan_to_num(t) * (len(ramp) - 1)
    lo = np.floor(pos).astype(int).clip(0, len(ramp) - 2)
    frac = (pos - lo)[:, None]
    return ramp[lo] * (1 - frac) + ramp[lo + 1] * frac


def view_frame(points):
    """Rotate points to their principal axes, viewing along the smallest."""
    centred = points - points.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    return centred @ vt.T


def face_normals_z(points, faces):
    """z component of the unit face normals (for shading a view along z)."""
    tri = points[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1e-12
    return normals[:, 2]


def render_faces(points2d, depth, faces, face_rgb, normals_z, ax=None, figsize=(4, 4)):
    """Draw faces back to front with simple shading, no axes."""
    from matplotlib.collections import PolyCollection

    shade = 0.55 + 0.45 * np.abs(normals_z)
    colors = np.clip(face_rgb * shade[:, None], 0, 1)
    order = np.argsort(depth[faces].mean(axis=1))
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    ax.add_collection(
        PolyCollection(
            points2d[faces[order]],
            facecolors=colors[order],
            edgecolors="none",
            antialiased=False,
        )
    )
    ax.set_xlim(points2d[:, 0].min(), points2d[:, 0].max())
    ax.set_ylim(points2d[:, 1].min(), points2d[:, 1].max())
    ax.set_aspect("equal")
    ax.axis("off")
    return fig


def _format_value(value):
    if isinstance(value, float):
        return f"{value:.3g}"
    if isinstance(value, (list, tuple)):
        return ", ".join(_format_value(v) for v in value)
    return str(value)


class Section:
    def __init__(self, title, status=None, subtitle="", anchor=None):
        self.title = title
        self.status = status
        self.subtitle = subtitle
        self.anchor = anchor or re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
        self.parts = []

    def add_html(self, content):
        self.parts.append(content)

    def add_text(self, text):
        self.parts.append(f"<p>{html.escape(text)}</p>")

    def add_heading(self, text):
        self.parts.append(f"<h3>{html.escape(text)}</h3>")

    def add_checks(self, result):
        rows = []
        order = {"fail": 0, "warn": 1, "pass": 2}
        checks = sorted(
            result.checks,
            key=lambda c: order[
                "pass"
                if c["passed"]
                else ("fail" if c["severity"] == "error" else "warn")
            ],
        )
        for c in checks:
            status = (
                "pass"
                if c["passed"]
                else ("fail" if c["severity"] == "error" else "warn")
            )
            rows.append(
                f"<tr><td>{status_chip(status)}</td>"
                f'<td><span class="check-name">{html.escape(c["name"])}</span>'
                f'<span class="desc">{html.escape(c.get("description", ""))}</span></td>'
                f"<td>{html.escape(c['message'])}</td>"
                f'<td class="muted">{html.escape(c["severity"])}</td></tr>'
            )
        self.parts.append(
            '<div class="table-wrap"><table class="checks"><thead><tr><th>Status</th>'
            "<th>Check</th><th>Result</th><th>Severity</th></tr></thead><tbody>"
            + "".join(rows)
            + "</tbody></table></div>"
        )

    def add_metrics(self, result, names=None):
        names = names or list(result.metrics)
        rows = "".join(
            f"<tr><td>{html.escape(n)}"
            f'<span class="desc">{html.escape(result.metrics[n].get("description", ""))}'
            f'</span></td><td class="num">{html.escape(_format_value(result.metrics[n]["value"]))}'
            "</td></tr>"
            for n in names
        )
        self.parts.append(
            '<details class="metrics"><summary>Metrics</summary><div class="table-wrap">'
            '<table><thead><tr><th>Metric</th><th class="num">Value</th></tr></thead>'
            f"<tbody>{rows}</tbody></table></div></details>"
        )

    def add_figure(self, fig, caption="", fmt="svg", extra_html=""):
        """Add a matplotlib figure (closed after rendering), as svg or png."""
        if fmt == "svg":
            body = f'<div class="svg-fig">{fig_to_svg(fig)}</div>'
        else:
            body = f'<img src="{fig_to_png(fig)}" alt="{html.escape(caption)}">'
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        self.parts.append(f"<figure>{extra_html}{body}{cap}</figure>")

    def render(self):
        chip = status_chip(self.status) if self.status else ""
        sub = (
            f'<p class="muted">{html.escape(self.subtitle)}</p>'
            if self.subtitle
            else ""
        )
        return (
            f'<section class="card" id="{self.anchor}"><header class="card-head">'
            f"<h2>{html.escape(self.title)}</h2>{chip}</header>{sub}"
            + "".join(self.parts)
            + "</section>"
        )


class QCReport:
    def __init__(self, title, entities=None, description=""):
        self.title = title
        self.entities = dict(entities or {})
        self.description = description
        self.summary = ""
        self.sections = []
        self.footer = ""

    def add_section(self, title, status=None, subtitle="", anchor=None):
        section = Section(title, status=status, subtitle=subtitle, anchor=anchor)
        self.sections.append(section)
        return section

    def add_summary(self, results, extra_columns=None):
        """
        Stat tiles and a table of results, with links to sections whose anchors
        are the result labels. extra_columns maps column names to a function of
        the result giving the cell text.
        """
        from lib.qc import overall_status

        extra_columns = extra_columns or {}
        n_checks = sum(len(r.checks) for r in results)
        n_err = sum(len(r.failed("error")) for r in results)
        n_warn = sum(len(r.failed("warning")) for r in results)
        status = overall_status(results)
        tiles = [
            ("Overall", status_chip(status)),
            ("Checks run", f'<span class="tile-value">{n_checks}</span>'),
            ("Errors", f'<span class="tile-value">{n_err}</span>'),
            ("Warnings", f'<span class="tile-value">{n_warn}</span>'),
        ]
        tiles_html = "".join(
            f'<div class="tile"><span class="tile-label">{label}</span>{value}</div>'
            for label, value in tiles
        )
        head = "".join(f"<th>{html.escape(c)}</th>" for c in extra_columns)
        rows = ""
        for r in results:
            anchor = re.sub(r"[^a-z0-9]+", "-", r.label.lower()).strip("-")
            cells = "".join(
                f"<td>{html.escape(str(fn(r)))}</td>" for fn in extra_columns.values()
            )
            rows += (
                f'<tr><td><a href="#{anchor}">{html.escape(r.label)}</a></td>'
                f"<td>{status_chip(r.status)}</td>{cells}"
                f'<td class="num">{len(r.failed("error"))}</td>'
                f'<td class="num">{len(r.failed("warning"))}</td></tr>'
            )
        self.summary = (
            f'<div class="tiles">{tiles_html}</div>'
            '<div class="card"><div class="table-wrap"><table class="summary-table"><thead>'
            "<tr><th>Input</th>"
            f'<th>Status</th>{head}<th class="num">Errors</th><th class="num">Warnings</th>'
            f"</tr></thead><tbody>{rows}</tbody></table></div></div>"
        )

    def render(self):
        ents = " · ".join(
            f"{k}-{html.escape(str(v))}" for k, v in self.entities.items()
        )
        now = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
        desc = (
            f'<p class="lede">{html.escape(self.description)}</p>'
            if self.description
            else ""
        )
        body = (
            f'<header class="page-head"><p class="eyebrow">HippUnfold QC</p>'
            f"<h1>{html.escape(self.title)}</h1>"
            f'<p class="muted">{ents} · generated {now}</p>{desc}</header>'
            + self.summary
            + "".join(s.render() for s in self.sections)
            + (f'<footer class="muted">{self.footer}</footer>' if self.footer else "")
        )
        title = f"{self.title} QC"
        if self.entities:
            title += " — " + " ".join(f"{k}-{v}" for k, v in self.entities.items())
        return _PAGE.format(title=html.escape(title), css=_CSS, body=body)

    def save(self, path):
        with open(path, "w") as f:
            f.write(self.render())


_PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{css}</style>
</head>
<body><main>{body}</main></body>
</html>
"""

_DARK_TOKENS = """
  color-scheme: dark;
  --page: #0d0d0d;
  --surface-1: #1a1a19;
  --surface-2: #232321;
  --text-primary: #f0efec;
  --text-secondary: #c3c2b7;
  --text-muted: #898781;
  --axis: #383835;
  --grid: #2c2c2a;
  --border: rgba(255,255,255,0.10);
  --series-1: {s1}; --series-2: {s2}; --series-3: {s3}; --series-4: {s4};
  --link: #6da7ec;
""".format(**{f"s{i + 1}": c for i, c in enumerate(SERIES_DARK)})

_CSS = (
    """
:root {
  color-scheme: light;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --page: #f9f9f7;
  --surface-1: #fcfcfb;
  --surface-2: #f3f2ef;
  --text-primary: #0b0b0b;
  --text-secondary: #52514e;
  --text-muted: #898781;
  --axis: #c3c2b7;
  --grid: #e1e0d9;
  --border: rgba(11,11,11,0.10);
  --series-1: #2a78d6; --series-2: #eb6834; --series-3: #1baf7a; --series-4: #eda100;
  --good: #0ca30c; --warning: #fab219; --critical: #d03b3b;
  --link: #256abf;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {"""
    + _DARK_TOKENS
    + """  }
}
:root[data-theme="dark"] {"""
    + _DARK_TOKENS
    + """}
* { box-sizing: border-box; }
html, body { margin: 0; }
body { background: var(--page); color: var(--text-primary); font: 14px/1.5 var(--font); }
main { max-width: 1120px; margin: 0 auto; padding: 32px 16px 64px; }
a { color: var(--link); }
h1 { font-size: 26px; line-height: 1.2; margin: 4px 0 6px; font-weight: 650; }
h2 { font-size: 17px; margin: 0; font-weight: 620; }
h3 { font-size: 14px; margin: 22px 0 8px; color: var(--text-secondary); font-weight: 600; }
p { margin: 6px 0; }
.eyebrow { text-transform: uppercase; letter-spacing: .08em; font-size: 11px;
  color: var(--text-muted); margin: 0; font-weight: 600; }
.muted { color: var(--text-muted); }
.lede { color: var(--text-secondary); max-width: 72ch; }
.page-head { margin-bottom: 20px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 12px; margin-bottom: 16px; }
.tile { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px;
  padding: 12px 14px; display: flex; flex-direction: column; gap: 6px; align-items: flex-start; }
.tile-label { font-size: 12px; color: var(--text-muted); }
.tile-value { font-size: 24px; font-weight: 600; line-height: 1.1; }
.card { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px;
  padding: 18px 20px; margin: 16px 0; }
.card-head { display: flex; gap: 12px; align-items: center; justify-content: space-between;
  flex-wrap: wrap; }
.chip { display: inline-flex; align-items: center; gap: 6px; font-size: 12px; font-weight: 600;
  padding: 2px 10px 2px 4px; border-radius: 999px; border: 1px solid var(--border);
  background: var(--surface-2); color: var(--text-primary); white-space: nowrap; }
.chip-icon { display: inline-grid; place-items: center; width: 18px; height: 18px;
  border-radius: 50%; color: #fff; font-size: 11px; }
.chip-pass .chip-icon { background: var(--good); }
.chip-warn .chip-icon { background: var(--warning); color: #0b0b0b; }
.chip-fail .chip-icon { background: var(--critical); }
.table-wrap { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th { text-align: left; font-weight: 600; color: var(--text-muted); font-size: 12px;
  border-bottom: 1px solid var(--axis); padding: 6px 10px 6px 0; }
td { border-bottom: 1px solid var(--grid); padding: 7px 10px 7px 0; vertical-align: top; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
.check-name { display: block; font-weight: 600; }
.desc { display: block; color: var(--text-muted); font-size: 12px; }
details.metrics { margin-top: 12px; }
details.metrics summary { cursor: pointer; color: var(--text-secondary); font-weight: 600; }
figure { margin: 18px 0 0; }
figure img, .svg-fig svg { width: 100%; height: auto; display: block; }
figcaption { color: var(--text-muted); font-size: 12px; margin-top: 6px; max-width: 80ch; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 12px;
  color: var(--text-secondary); margin-bottom: 6px; }
.key { display: inline-flex; align-items: center; gap: 6px; }
.swatch { width: 12px; height: 12px; border-radius: 3px; display: inline-block; }
.panels { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
  gap: 16px; }
.panel { background: var(--surface-2); border-radius: 10px; padding: 12px; }
.panel h4 { margin: 0 0 6px; font-size: 13px; font-weight: 600; }
.panel img { width: 100%; height: auto; display: block; margin-top: 8px; }
.ramp { height: 8px; border-radius: 4px; }
.ramp-ap { background: linear-gradient(90deg, """
    + ", ".join(BLUE_RAMP)
    + """); }
.ramp-pd { background: linear-gradient(90deg, """
    + ", ".join(ORANGE_RAMP)
    + """); }
.ramp-div { background: linear-gradient(90deg, """
    + ", ".join(DIVERGING_RAMP)
    + """); }
.ramp-labels { display: flex; justify-content: space-between; font-size: 11px;
  color: var(--text-muted); }
code { font-size: 12px; background: var(--surface-2); padding: 1px 4px; border-radius: 4px; }
footer { margin-top: 32px; font-size: 12px; }
.summary-table { min-width: 620px; }
@media (max-width: 640px) {
  .card { padding: 14px; }
  table.checks thead { display: none; }
  table.checks tr { display: block; border-bottom: 1px solid var(--grid); padding: 8px 0; }
  table.checks td { display: block; border: 0; padding: 2px 0; }
  table.checks td:last-child { display: none; }
}
footer dt { font-weight: 600; color: var(--text-secondary); margin-top: 8px; }
footer dd { margin: 2px 0 0; }
"""
)
