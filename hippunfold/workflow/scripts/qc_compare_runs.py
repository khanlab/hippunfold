#!/usr/bin/env python
"""
Compare the QC metrics of two HippUnfold runs (e.g. before and after a change
to the workflow) as a self-contained HTML report.

    python hippunfold/workflow/scripts/qc_compare_runs.py BEFORE_DIR AFTER_DIR \
        --qc surface --labels main new -o comparison.html

Each *_desc-<qc>_qcmetrics.json in BEFORE_DIR is paired with the file at the
same path in AFTER_DIR, and their checks and metrics are shown side by side.
For the surface QC, defect maps are added where the native midthickness
surfaces were kept (e.g. runs with --notemp).
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import nibabel as nib
from lib.qc import QCResult
from lib.qc_report import QCReport
from lib.qc_surface import (
    COMPARE_METRICS,
    DEFECT_MAP_CAPTION,
    defect_panel,
    defect_table,
)

UNFOLD_METRICS = [
    "unfold_flipped_faces",
    "unfold_collapsed_faces",
    "unfold_distorted_fraction",
    "unfold_log2_distortion_iqr",
]


def surface_path(run_dir, json_path):
    """Native corobl midthickness of a surface QC JSON, if it was kept."""
    name = json_path.name
    ents = dict(re.findall(r"(sub|ses|hemi|label)-([a-zA-Z0-9]+)", name))
    prefix = "_".join(f"{k}-{ents[k]}" for k in ("sub", "ses") if k in ents)
    sub_dir = Path(run_dir) / f"sub-{ents['sub']}"
    if "ses" in ents:
        sub_dir = sub_dir / f"ses-{ents['ses']}"
    path = (
        sub_dir
        / "surf"
        / f"{prefix}_hemi-{ents['hemi']}_space-corobl_den-native_label-"
        f"{ents['label']}_midthickness.surf.gii"
    )
    return path if path.exists() else None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("before", type=Path, help="output dir of the first run")
    parser.add_argument("after", type=Path, help="output dir of the second run")
    parser.add_argument("--qc", default="surface", help="QC name (default: surface)")
    parser.add_argument(
        "--labels", nargs=2, default=["Before", "After"], help="names of the runs"
    )
    parser.add_argument("-o", "--output", type=Path, default=Path("qc_comparison.html"))
    args = parser.parse_args()

    pairs = []
    for before in sorted(
        args.before.glob(f"sub-*/**/qc/*_desc-{args.qc}_qcmetrics.json")
    ):
        after = args.after / before.relative_to(args.before)
        if after.exists():
            pairs.append((before, after))
        else:
            print(f"No match in {args.after} for {before.name}, skipped")
    if not pairs:
        sys.exit(f"No desc-{args.qc} QC metrics found in both runs")

    report = QCReport(
        f"{args.qc} QC comparison",
        entities={"runs": " vs ".join(args.labels)},
        description=(
            f"Checks and metrics of the desc-{args.qc} QC, for {args.labels[0]} "
            f"({args.before}) and {args.labels[1]} ({args.after})."
        ),
    )
    afters = []
    failed_before = {}  # by id of the after result
    for before_path, after_path in pairs:
        before, after = QCResult.load(before_path), QCResult.load(after_path)
        subject = re.search(r"sub-[a-zA-Z0-9]+(_ses-[a-zA-Z0-9]+)?", before_path.name)
        after.entities = {"sub": subject.group(0)[4:], **after.entities}
        afters.append(after)
        failed_before[id(after)] = len(before.failed())
        section = report.add_section(after.label, status=after.status)
        names = COMPARE_METRICS + UNFOLD_METRICS if args.qc == "surface" else None
        section.add_comparison(
            [(args.labels[0], before), (args.labels[1], after)], names
        )
        if args.qc != "surface":
            continue

        panels = ""
        for label, run_dir, json_path, result in (
            (args.labels[0], args.before, before_path, before),
            (args.labels[1], args.after, after_path, after),
        ):
            surf = surface_path(run_dir, json_path)
            if surf is None:
                continue
            points, faces = nib.load(surf).agg_data()
            n = result.metrics["defects"]["value"]
            panels += defect_panel(f"{label}: {n} defects", result, points, faces)
        if panels:
            section.add_heading("Defects")
            section.add_html(f'<div class="panels">{panels}</div>{DEFECT_MAP_CAPTION}')
        for label, result in ((args.labels[0], before), (args.labels[1], after)):
            section.add_html(
                f'<details class="metrics"><summary>{label} defects '
                f"({result.metrics['defects']['value']})</summary>"
                f"{defect_table(result.data['defects'], result.metrics['defects']['value'])}</details>"
            )

    report.add_summary(
        afters,
        extra_columns={
            f"Failed checks ({args.labels[0]} → {args.labels[1]})": lambda r: (
                f"{failed_before[id(r)]} → {len(r.failed())}"
            ),
        },
    )
    report.save(args.output)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
