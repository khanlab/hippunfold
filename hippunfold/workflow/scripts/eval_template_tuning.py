"""Evaluate one antsMultivariateTemplateConstruction2.sh run (gen_atlas_reg_ants_tune).

Metrics, all computed in the average (template) space:
  - loo_dice_*: each subject's subfield labels are warped to the template and compared
    (Dice, averaged over labels) to the majority vote of the other subjects' warped
    labels. Subfields are not used to drive the registration, so this is an
    independent measure of anatomical correspondence. The same measure on the
    unwarped images is reported as loo_dice_unwarped_* for reference.
  - loo_corr_<metric>: Pearson correlation of each subject's warped metric with the
    mean of the other subjects' warped metric (leave-one-out template).
  - sharpness_<metric>: mean gradient magnitude of the template, divided by its std.
  - jac_*: from the subject->template warps: fraction of voxels with det(J) <= 0
    (folding), std of log det(J), and mean / max displacement magnitude (mm).

Can be run from snakemake (script: directive) or from the command line:
  python eval_template_tuning.py --template_dir DIR --metrics_csv CSV --root ROOT \
      --out_tsv OUT.tsv [--out_png OUT.png] [--subfields_suffix subfieldsfromnative]
"""

import argparse
import glob
import json
import os
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


def load2d(path):
    data = np.asarray(nib.load(path).dataobj)
    return np.squeeze(data)


def find_one(pattern):
    matches = sorted(glob.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"expected 1 file matching {pattern}, found {len(matches)}")
    return matches[0]


def strip_nii(path):
    return Path(path).name.removesuffix(".nii.gz")


def majority_vote(label_maps, labels):
    counts = np.stack([(label_maps == lbl).sum(axis=0) for lbl in labels], axis=-1)
    return labels[np.argmax(counts, axis=-1)]


def dice(a, b, labels):
    scores = []
    for lbl in labels:
        if lbl == 0:
            continue
        ina, inb = a == lbl, b == lbl
        denom = ina.sum() + inb.sum()
        if denom > 0:
            scores.append(2.0 * np.logical_and(ina, inb).sum() / denom)
    return float(np.mean(scores))


def loo_dice(label_maps):
    labels = np.unique(label_maps)
    per_subj = []
    for i in range(len(label_maps)):
        others = np.delete(label_maps, i, axis=0)
        per_subj.append(dice(label_maps[i], majority_vote(others, labels), labels))
    return per_subj


def loo_corr(imgs):
    corrs = []
    for i in range(len(imgs)):
        loo_mean = np.delete(imgs, i, axis=0).mean(axis=0)
        x, y = imgs[i].ravel(), loo_mean.ravel()
        ok = np.isfinite(x) & np.isfinite(y)
        corrs.append(float(np.corrcoef(x[ok], y[ok])[0, 1]))
    return corrs


def sharpness(img):
    gx, gy = np.gradient(img.astype(np.float64))
    std = img.std()
    return float(np.hypot(gx, gy).mean() / std) if std > 0 else float("nan")


def warp_labels(lbl_file, warp_file, tmpdir):
    out = os.path.join(tmpdir, "warped_lbl.nii.gz")
    subprocess.run(
        [
            "antsApplyTransforms",
            "-d",
            "2",
            "-i",
            lbl_file,
            "-r",
            lbl_file,
            "-o",
            out,
            "-t",
            warp_file,
            "-n",
            "GenericLabel",
        ],
        check=True,
    )
    return np.rint(load2d(out)).astype(np.int32)


def jacobian_stats(warp_file, tmpdir):
    out = os.path.join(tmpdir, "jac.nii.gz")
    subprocess.run(
        ["CreateJacobianDeterminantImage", "2", warp_file, out, "0", "0"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    jac = load2d(out).astype(np.float64)
    disp = np.squeeze(np.asarray(nib.load(warp_file).dataobj)).astype(np.float64)
    mag = np.linalg.norm(disp, axis=-1)
    pos = jac[jac > 0]
    return {
        "folding_frac": float((jac <= 0).mean()),
        "logjac_std": float(np.log(pos).std()) if pos.size else float("nan"),
        "disp_mean_mm": float(mag.mean()),
        "disp_max_mm": float(mag.max()),
    }


def save_png(out_png, template_imgs, metric_names, voted, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(template_imgs) + (voted is not None)
    fig, axes = plt.subplots(1, n, figsize=(1.6 * n, 3.2))
    axes = np.atleast_1d(axes)
    for ax, img, name in zip(axes, template_imgs, metric_names):
        lo, hi = np.nanpercentile(img, [2, 98])
        ax.imshow(img.T, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        ax.set_title(name, fontsize=8)
    if voted is not None:
        axes[-1].imshow(voted.T, origin="lower", cmap="tab10", vmin=0, vmax=9)
        axes[-1].set_title("subfields (vote)", fontsize=8)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(title, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def evaluate(
    template_dir,
    metrics_csv,
    out_tsv,
    out_png,
    params,
    lbl_files=None,
    root=None,
    subfields_suffix=None,
):
    rows = [
        line.strip().split(",")
        for line in Path(metrics_csv).read_text().splitlines()
        if line.strip()
    ]
    metric_names = [strip_nii(f).rsplit("_", 1)[-1] for f in rows[0]]

    results = dict(params)
    results["n_inputs"] = len(rows)

    warps = [
        find_one(f"{template_dir}/input*-{strip_nii(row[0])}-1Warp.nii.gz")
        for row in rows
    ]

    # metric agreement and template sharpness
    template_imgs = []
    for k, name in enumerate(metric_names):
        warped = np.stack(
            [
                load2d(
                    find_one(
                        f"{template_dir}/input*-modality{k}-{strip_nii(row[k])}-WarpedToTemplate.nii.gz"
                    )
                )
                for row in rows
            ]
        ).astype(np.float64)
        results[f"loo_corr_{name}"] = float(np.mean(loo_corr(warped)))
        tpl = load2d(f"{template_dir}/template{k}.nii.gz").astype(np.float64)
        template_imgs.append(tpl)
        results[f"sharpness_{name}"] = sharpness(tpl)
    results["loo_corr_mean"] = float(
        np.mean([results[f"loo_corr_{name}"] for name in metric_names])
    )

    voted = None
    with tempfile.TemporaryDirectory() as tmpdir:
        jac = [jacobian_stats(w, tmpdir) for w in warps]
        for key in jac[0]:
            vals = [j[key] for j in jac]
            agg = np.max if key in ("folding_frac", "disp_max_mm") else np.mean
            results[f"jac_{key}"] = float(agg(vals))

        if lbl_files is None and subfields_suffix:
            # same naming as the metric images in the csv
            lbl_files = [
                os.path.join(
                    root, row[0].rsplit("_", 1)[0] + f"_{subfields_suffix}.nii.gz"
                )
                for row in rows
            ]
        if lbl_files:
            if len(lbl_files) != len(rows):
                raise ValueError("Number of subfield images and csv rows must match.")
            unwarped = np.stack(
                [np.rint(load2d(f)).astype(np.int32) for f in lbl_files]
            )
            warped = np.stack(
                [warp_labels(f, w, tmpdir) for f, w in zip(lbl_files, warps)]
            )
            d_unwarped = loo_dice(unwarped)
            d_warped = loo_dice(warped)
            results["loo_dice_mean"] = float(np.mean(d_warped))
            results["loo_dice_min"] = float(np.min(d_warped))
            results["loo_dice_unwarped_mean"] = float(np.mean(d_unwarped))
            results["loo_dice_per_input"] = json.dumps([round(d, 4) for d in d_warped])
            voted = majority_vote(warped, np.unique(warped))

    import pandas as pd

    pd.DataFrame([results]).to_csv(out_tsv, sep="\t", index=False)

    if out_png:
        title = " ".join(f"{k}-{v}" for k, v in params.items())
        save_png(out_png, template_imgs, metric_names, voted, title)


if "snakemake" in globals():
    evaluate(
        template_dir=snakemake.input.avgtemplate_dir,
        metrics_csv=snakemake.input.metrics_csv,
        lbl_files=list(snakemake.input.subfields),
        out_tsv=snakemake.output.tsv,
        out_png=snakemake.output.png,
        params=dict(snakemake.params.tune_params),
    )
else:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template_dir", required=True)
    parser.add_argument("--metrics_csv", required=True)
    parser.add_argument(
        "--root", required=True, help="dir that paths in metrics_csv are relative to"
    )
    parser.add_argument("--subfields_suffix", default="subfieldsfromnative")
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--out_png", default=None)
    parser.add_argument(
        "--name", default=None, help="run name to store in the tsv (default: dir name)"
    )
    args = parser.parse_args()
    evaluate(
        template_dir=args.template_dir,
        metrics_csv=args.metrics_csv,
        root=args.root,
        subfields_suffix=args.subfields_suffix or None,
        out_tsv=args.out_tsv,
        out_png=args.out_png,
        params={"run": args.name or Path(args.template_dir).name},
    )
