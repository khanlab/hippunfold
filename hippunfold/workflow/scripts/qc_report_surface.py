import re

import nibabel as nib
import numpy as np
from lib.qc import QCResult
from lib.qc_report import QCReport, fig_to_png
from lib.qc_surface import (
    COMPARE_METRICS,
    DEFAULT_THRESHOLDS,
    DEFECT_MAP_CAPTION,
    DEFECT_TYPES,
    STAGE_DESCS,
    STAGE_LABELS,
    defect_panel,
    defect_table,
    face_geometry,
    plot_distortion,
    plot_histograms,
    unfold_distortion,
)
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


def topology_banner(result):
    m = {k: result.metrics[k]["value"] for k in ("V", "E", "F", "chi")}
    c, b, g = (
        result.metrics[k]["value"] for k in ("components", "boundary_loops", "genus")
    )
    disk = c == 1 and b == 1 and g == 0
    verdict = (
        "a disk, as expected"
        if disk
        else f"not a disk ({c} components, {b} boundary loops, genus {g:g})"
    )
    return (
        '<div class="panel" style="margin-top:12px"><p style="font-size:15px;margin:0">'
        f"χ = V − E + F = {m['V']:,} − {m['E']:,} + {m['F']:,} = <b>{m['chi']}</b> "
        f"= 2c − 2g − b with c = {c}, g = {g:g}, b = {b}</p>"
        f'<p class="muted" style="margin:4px 0 0">Topologically {verdict}.</p></div>'
    )


results = [QCResult.load(path) for path in snakemake.input.metrics]
stage_results = {}
for path in snakemake.input.stage_metrics:
    r = QCResult.load(path)
    stage_results[(r.data["stage"], r.label)] = r
stages = list(snakemake.params.stages)


def stages_of(result):
    """(stage, result) for each intermediate stage of a final result, then final"""
    found = [
        (s, stage_results[(s, result.label)])
        for s in stages
        if (s, result.label) in stage_results
    ]
    return found + [("final", result)]


def stage_surface(stage, entities):
    if stage == "final":
        return find(snakemake.input.surf, entities)
    paths = [p for p in snakemake.input.stage_surf if f"desc-{stage}_" in p]
    return find(paths, entities)


subj_entities = {
    k: v for k, v in snakemake.wildcards.items() if k in ("subject", "session")
}
report = QCReport(
    "Surface topology and quality",
    entities=subj_entities,
    description=(
        "Each hippocampal and dentate surface should be a compact, orientable "
        "2-manifold homeomorphic to a closed disk: one connected component, genus 0 "
        "and one boundary loop, so χ = V − E + F = 1. Checks below flag handles, "
        "holes, non-manifold or misoriented elements, self-intersections, poor "
        "triangles, and folded-over or collapsed faces in the unfolded (AP, PD) "
        "mapping. Every defect found is listed, located on the surface, and the "
        "surface is compared before and after processing steps such as repair."
    ),
)


def fmt_topology(r):
    return "{} / {:g} / {}".format(
        r.metrics["chi"]["value"],
        r.metrics["genus"]["value"],
        r.metrics["boundary_loops"]["value"],
    )


report.add_summary(
    results,
    extra_columns={
        "χ / g / b": fmt_topology,
        "Defects": lambda r: " → ".join(
            str(sr.metrics["defects"]["value"]) for _, sr in stages_of(r)
        ),
        "Self-intersecting": lambda r: sum(
            v["value"]
            for k, v in r.metrics.items()
            if k.startswith("self_intersecting_faces")
        ),
        "Flipped (unfold)": lambda r: r.metrics["unfold_flipped_faces"]["value"],
    },
)


for result in results:
    logger.info(f"Adding section for {result.label}")
    section = report.add_section(result.label, status=result.status)
    section.add_html(topology_banner(result))
    section.add_checks(result)
    section.add_metrics(result)

    surf = nib.load(find(snakemake.input.surf, result.entities))
    points, faces = surf.agg_data()
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    coords_ap = nib.load(find(snakemake.input.coords_ap, result.entities)).agg_data()
    coords_pd = nib.load(find(snakemake.input.coords_pd, result.entities)).agg_data()
    defects = result.data["defects"]

    staged = stages_of(result)
    if len(staged) > 1:
        section.add_heading("Before and after surface processing")
        section.add_html(
            '<p class="muted" style="font-size:12px">'
            + "; ".join(
                f"<b>{STAGE_LABELS.get(s, s)}</b>: {STAGE_DESCS[s]}"
                for s, _ in staged[:-1]
            )
            + ". Unfolded-mapping checks need the final surface, so are only run "
            "there.</p>"
        )
        section.add_comparison(
            [(STAGE_LABELS.get(s, s), r) for s, r in staged], names=COMPARE_METRICS
        )

    section.add_heading("Defects")
    panels = ""
    for stage, stage_result in staged:
        if stage == "final":
            p, f = points, faces
        else:
            p, f = nib.load(stage_surface(stage, result.entities)).agg_data()
        n = stage_result.metrics["defects"]["value"]
        panels += defect_panel(
            f"{STAGE_LABELS.get(stage, stage)}: {n} defects", stage_result, p, f
        )
    section.add_html(f'<div class="panels">{panels}</div>{DEFECT_MAP_CAPTION}')
    section.add_html(defect_table(defects, result.metrics["defects"]["value"]))
    for stage, stage_result in staged[:-1]:
        section.add_html(
            f'<details class="metrics"><summary>{STAGE_LABELS.get(stage, stage)} '
            f"defects ({stage_result.metrics['defects']['value']})</summary>"
            f"{defect_table(stage_result.data['defects'], stage_result.metrics['defects']['value'])}</details>"
        )

    _, _, flipped, log_ratio = unfold_distortion(points, faces, coords_ap, coords_pd)
    figs = plot_distortion(points, faces, coords_ap, coords_pd, flipped, log_ratio)
    ramp_key = (
        '<div class="ramp ramp-div"></div><div class="ramp-labels"><span>¼×</span>'
        "<span>1×</span><span>4×</span></div>"
    )
    panels = "".join(
        f'<div class="panel"><h4>{title}</h4>{ramp_key}<img src="{fig_to_png(fig)}" '
        f'alt="{title}"></div>'
        for title, fig in zip(
            [
                "Areal distortion, native (face-on)",
                "Areal distortion, unfolded (PD →, AP ↑)",
            ],
            figs,
        )
    )
    section.add_heading("Unfolded mapping")
    section.add_html(
        f'<div class="panels">{panels}</div><p class="muted" style="font-size:12px">'
        "Colour: each face's share of the unfolded area relative to its share of the "
        "native area. Grey: faces collapsed to zero area. Circles: faces folded over "
        "in the unfolded mapping.</p>"
    )

    min_angles = {}
    for name, paths in (
        ("midthickness", snakemake.input.surf),
        ("inner", snakemake.input.inner),
        ("outer", snakemake.input.outer),
    ):
        p, f = nib.load(find(paths, result.entities)).agg_data()
        min_angles[name] = face_geometry(np.asarray(p, dtype=float), f)[1]
    section.add_heading("Triangle quality and distortion")
    section.add_figure(
        plot_histograms(
            min_angles,
            log_ratio,
            DEFAULT_THRESHOLDS["sliver_angle"],
            DEFAULT_THRESHOLDS["distortion_log2"],
        ),
        caption=(
            "Left: smallest angle of each triangle, per surface; dashed line is the "
            "sliver threshold. Right: log2 areal distortion of the unfolded mapping "
            "(collapsed faces excluded); dashed lines mark the distortion threshold."
        ),
    )

checks_help = {
    "single_component": "error — the surface is one connected piece.",
    "has_boundary": "error — the surface has at least one boundary loop.",
    "orientable": "error if the surface is non-orientable; warning if some faces are "
    "wound opposite to the rest.",
    "manifold": "warning — no edge with more than 2 faces, no vertex joining "
    "separate fans of faces.",
    "genus": "warning — no handles (g = 0), e.g. from grey matter bridges between "
    "adjacent folds.",
    "holes": "warning — a single boundary loop; holes are free boundaries in the "
    "Laplace-Beltrami solve.",
    "degenerate_faces": "warning — no zero-area or duplicate faces on any surface.",
    "sliver_faces": "warning — no triangle with a min angle below "
    f"{DEFAULT_THRESHOLDS['sliver_angle']:g}° (default) on any surface.",
    "self_intersections": "warning — the midthickness, inner and outer surfaces do "
    "not intersect themselves.",
    "unfold_flipped": "warning — no faces folded over in the unfolded (AP, PD) "
    "mapping.",
    "unfold_collapsed": "warning — no faces with zero unfolded area.",
    "unfold_distortion": "warning — at most "
    f"{DEFAULT_THRESHOLDS['max_distorted_fraction']:.0%} (default) of faces change "
    f"area share by more than {2 ** DEFAULT_THRESHOLDS['distortion_log2']:g}×.",
}
types_help = "".join(
    f"<dt>{t}</dt><dd>{desc} ({cat})</dd>" for t, (cat, desc) in DEFECT_TYPES.items()
)
report.footer = (
    "<p><b>Checks</b></p><dl>"
    + "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in checks_help.items())
    + "</dl><p><b>Defect types</b></p><dl>"
    + types_help
    + "</dl><p>Failed error checks stop the run after this report is written; use "
    "<code>--qc_warn_only</code> to continue anyway.</p>"
)

report.save(snakemake.output.html)
logger.info(f"Report saved to {snakemake.output.html}")
