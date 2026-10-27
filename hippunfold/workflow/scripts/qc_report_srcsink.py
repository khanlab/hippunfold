import re

import nibabel as nib
from lib.qc import QCResult
from lib.qc_report import QCReport, fig_to_png, legend
from lib.qc_srcsink import DEFAULT_THRESHOLDS, plot_loop, plot_surfaces
from lib.surface import TERMINAL_NAMES
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def entities_of(path):
    return dict(re.findall(r"(hemi|label)-([a-zA-Z0-9]+)", path))


def find(paths, entities):
    """Path whose hemi/label entities match"""
    matches = [p for p in paths if entities_of(p) == entities]
    if len(matches) != 1:
        raise ValueError(f"Expected one input for {entities}, found {matches}")
    return matches[0]


results = [
    QCResult.load(path, warn_checks=snakemake.params.warn_checks)
    for path in snakemake.input.metrics
]

subj_entities = {
    k: v for k, v in snakemake.wildcards.items() if k in ("subject", "session")
}
report = QCReport(
    "Source/sink terminals",
    entities=subj_entities,
    description=(
        "The AP and PD source/sink terminals set the boundary conditions for the "
        "Laplace-Beltrami coordinates of each surface. They are fit as four "
        "contiguous arcs of the surface's outer boundary loop, alternating AP and "
        "PD so that each src is opposite its sink. Checks below flag terminals that "
        "are out of order for the hemisphere, fragmented, or far from their "
        "volumetric labels."
    ),
)
report.add_summary(
    results,
    extra_columns={
        "Terminal order": lambda r: r.metrics["terminal_order"]["value"],
        "Loop vertices": lambda r: r.metrics["loop_vertices"]["value"],
        "Relabelled": lambda r: f"{r.metrics['relabelled_fraction']['value']:.1%}",
    },
)

terminal_legend = legend([(name, k) for k, name in enumerate(TERMINAL_NAMES)])

for result in results:
    logger.info(f"Adding section for {result.label}")
    section = report.add_section(result.label, status=result.status)
    section.add_checks(result)
    section.add_metrics(result)

    section.add_heading("Terminals along the boundary loop")
    section.add_figure(
        plot_loop(result),
        caption=(
            "Lines: distance of each loop vertex to each volumetric terminal. "
            "Shading and top strip: fitted terminal arcs. Bottom strip: nearest "
            "terminal before fitting. Differences between the strips are the "
            "vertices the fit changed."
        ),
        extra_html=terminal_legend,
    )

    surf = nib.load(find(snakemake.input.surf, result.entities))
    points, faces = surf.agg_data()
    coords_ap = nib.load(find(snakemake.input.coords_ap, result.entities)).agg_data()
    coords_pd = nib.load(find(snakemake.input.coords_pd, result.entities)).agg_data()
    figs = plot_surfaces(result, points, faces, coords_ap, coords_pd)

    titles = ["Terminals on the boundary", "AP coordinate", "PD coordinate"]
    ramp_key = (
        '<div class="ramp ramp-{}"></div><div class="ramp-labels">'
        "<span>0 (src)</span><span>1 (sink)</span></div>"
    )
    keys = [terminal_legend, ramp_key.format("ap"), ramp_key.format("pd")]
    panels = "".join(
        f'<div class="panel"><h4>{title}</h4>{key}<img src="{fig_to_png(fig)}" '
        f'alt="{title}"></div>'
        for title, key, fig in zip(titles, keys, figs)
    )
    section.add_heading("Surface (midthickness, viewed face-on)")
    section.add_html(
        f'<div class="panels">{panels}</div><p class="muted" style="font-size:12px">'
        "Laplace-Beltrami coordinates solved from these terminals. Smooth gradients "
        "running src to sink are expected; abrupt jumps or plateaus point to "
        "misplaced terminals.</p>"
    )

checks_help = {
    "terminal_order_hemisphere": "error — terminal order along the loop (in face "
    "winding order) is AP src → PD src → AP sink → PD sink for hemi-L, and the "
    "mirror image for hemi-R.",
    "nearest_terminal_fragmentation": "warning — nearest-terminal labels, before "
    "fitting, already form 4 arcs.",
    "relabelled_fraction": "warning — at most "
    f"{DEFAULT_THRESHOLDS['max_relabelled_fraction']:.0%} (default) of the loop "
    "changed by the fit.",
    "ap_terminal_extent": "warning — each AP terminal covers at most "
    f"{DEFAULT_THRESHOLDS['max_ap_fraction']:.0%} (default) of the loop.",
    "terminal_distance": "warning — median distance of each arc to its "
    f"volumetric terminal is at most {DEFAULT_THRESHOLDS['max_terminal_distance']} "
    "mm (default).",
    "band_connected": "warning — every boundary band vertex is connected to the loop.",
}
report.footer = (
    "<dl>"
    + "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in checks_help.items())
    + "</dl><p>Failed error checks stop the run after this report is written; use "
    "<code>--qc_warn_only</code> to continue anyway.</p>"
)

report.save(snakemake.output.html)
logger.info(f"Report saved to {snakemake.output.html}")
