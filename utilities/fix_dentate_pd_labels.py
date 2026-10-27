#!/usr/bin/env python
"""
Check the dentate PD src (9) / sink (10) convention of tissue segmentations
(as used with --modality dsegtissue) and write corrected copies.

The src should be the exposed edge of the dentate gyrus at the surface of the
hippocampus (margo denticulatus) and the sink its buried end inside the hilus
(see hippunfold/resources/label_lut/manual_dseg.tsv). Segmentations with the
labels reversed are written with 9 and 10 exchanged; others are copied as is.

    fix_dentate_pd_labels.py IN_DSEG [IN_DSEG ...] --out-dir OUT_DIR
    fix_dentate_pd_labels.py IN_DSEG [IN_DSEG ...]          (check only)
"""

import argparse
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hippunfold/workflow"))
from lib.dentate_pd_labels import (
    BACKGROUND,
    check_convention,
    swap_src_sink,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dsegs", nargs="+", type=Path)
    parser.add_argument(
        "--out-dir", type=Path, help="write corrected copies here (same file names)"
    )
    args = parser.parse_args()

    n_swapped = 0
    for path in args.dsegs:
        img = nib.load(path)
        dseg = np.asarray(img.dataobj)
        convention, d = check_convention(dseg)
        bg = (
            f"background next to src {d['src_neighbours'].get(BACKGROUND, 0):.0%}, "
            f"next to sink {d['sink_neighbours'].get(BACKGROUND, 0):.0%}"
            if convention != "missing"
            else "no src/sink labels"
        )
        print(f"{path.name}: {convention} ({bg})")
        if args.out_dir is None:
            continue
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out = args.out_dir / path.name
        if convention == "swapped":
            n_swapped += 1
            fixed = swap_src_sink(dseg).astype(img.get_data_dtype())
            nib.save(nib.Nifti1Image(fixed, img.affine, img.header), out)
            print(f"  -> {out} (9 and 10 exchanged)")
        else:
            shutil.copy2(path, out)
    if args.out_dir is not None:
        print(f"{n_swapped} of {len(args.dsegs)} segmentations swapped")


if __name__ == "__main__":
    main()
