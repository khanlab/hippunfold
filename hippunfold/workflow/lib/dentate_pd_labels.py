"""
Check the convention of the dentate PD src (9) and sink (10) labels of a tissue
segmentation, and swap them if they are reversed.

The templates place the src at the exposed edge of the dentate gyrus at the
surface of the hippocampus (margo denticulatus), so it borders the background
(0), and the sink at its buried end inside the hilus, enclosed by CA4 grey
matter (1) and the SRLM (2). Reversed labels mirror the unfolded dentate.
"""

import numpy as np
from scipy import ndimage as ndi

SRC, SINK, BACKGROUND = 9, 10, 0


def neighbour_fractions(dseg, label, ring=2):
    """fraction of each label among the voxels within `ring` voxels of `label`"""
    mask = dseg == label
    if not mask.any():
        return {}
    shell = ndi.binary_dilation(mask, iterations=ring) & ~mask
    values, counts = np.unique(dseg[shell], return_counts=True)
    return {int(v): float(c / counts.sum()) for v, c in zip(values, counts)}


def check_convention(dseg, ring=2):
    """
    Returns (convention, details), where convention is "template" (src borders
    the background more than the sink does), "swapped", or "missing" (no src
    or sink label).
    """
    src = neighbour_fractions(dseg, SRC, ring)
    sink = neighbour_fractions(dseg, SINK, ring)
    details = {"src_neighbours": src, "sink_neighbours": sink}
    if not src or not sink:
        return "missing", details
    margin = src.get(BACKGROUND, 0.0) - sink.get(BACKGROUND, 0.0)
    details["background_margin"] = margin
    return ("template" if margin > 0 else "swapped"), details


def swap_src_sink(dseg):
    """copy of dseg with the dentate PD src and sink labels exchanged"""
    out = dseg.copy()
    out[dseg == SRC] = SINK
    out[dseg == SINK] = SRC
    return out
