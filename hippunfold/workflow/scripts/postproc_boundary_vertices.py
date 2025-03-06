import pyvista as pv
import numpy as np
import nibabel as nib
import nibabel.gifti as gifti
from collections import Counter
from scipy.signal import argrelextrema
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def write_src_sink_label_gifti(label_data, out_label_gifti, structure_metadata):

    # Create a Label Table (LUT)
    label_table = gifti.GiftiLabelTable()

    # Define Background label (key 0)
    background_label = gifti.GiftiLabel(
        key=0, red=1.0, green=1.0, blue=1.0, alpha=0.0
    )  # Transparent
    background_label.label = "Background"
    label_table.labels.append(background_label)

    # Define src label (key 1)
    src_label = gifti.GiftiLabel(key=1, red=1.0, green=0.0, blue=0.0, alpha=1.0)
    src_label.label = "Source"
    label_table.labels.append(src_label)

    # Define sink label (key 2)
    sink_label = gifti.GiftiLabel(key=2, red=0.0, green=0.0, blue=1.0, alpha=1.0)
    sink_label.label = "Sink"
    label_table.labels.append(sink_label)

    # write the data and the label table to gifti
    gii_img = gifti.GiftiImage(
        darrays=[gifti.GiftiDataArray(label_data, intent="NIFTI_INTENT_LABEL")],
        labeltable=label_table,
    )

    # set structure metadata
    gii_img.meta["AnatomicalStructurePrimary"] = structure_metadata

    nib.save(gii_img, out_label_gifti)


nmin = snakemake.params.min_terminal_vertices
max_iterations = snakemake.params.max_iterations  # Prevent infinite loops
shifting_epsilon = snakemake.params.shifting_epsilon

logger.info("Loading surface from GIFTI...")


edges = nib.load(snakemake.input.edges).agg_data()
ap_src = nib.load(snakemake.input.ap_src).agg_data()
ap_sink = nib.load(snakemake.input.ap_sink).agg_data()
pd_src = nib.load(snakemake.input.pd_src).agg_data()
pd_sink = nib.load(snakemake.input.pd_sink).agg_data()

# get structure metadata from edges
structure_metadata = nib.load(snakemake.input.edges).meta["AnatomicalStructurePrimary"]


logger.info("Assigning labels apsrc, apsink, pdsrc, pdsink")

distances = np.vstack(
    (ap_src[edges == 1], ap_sink[edges == 1], pd_src[edges == 1], pd_sink[edges == 1])
).T
shifting_factors = np.zeros((4))
num_labels = 4  # starts at 0

for _ in range(max_iterations):
    # shift distances
    shifted_distances = distances - shifting_factors

    # Assign labels based on min scaled distance
    labels = np.argmin(shifted_distances, axis=1)

    # Count occurrences per label
    unique, counts = np.unique(labels, return_counts=True)
    label_counts = dict(zip(unique, counts))
    # Ensure all labels in range(num_labels) are present, setting missing ones to 0
    label_counts = {k: label_counts.get(k, 0) for k in range(num_labels)}

    logger.info(f"Iteration {_}, have label_counts={label_counts}")

    # Check if all labels meet nmin
    if all(count >= nmin for count in label_counts.values()):
        break  # Stop if all labels are sufficiently represented

    # Update scaling factors for underrepresented labels
    for k in range(num_labels):
        if label_counts[k] < nmin:
            shifting_factors[
                k
            ] += shifting_epsilon  # Increase competitiveness of the label
            logger.info(
                f"Shifting distances of label {k} to bump up competitiveness, factor = {shifting_factors[k]}"
            )

# Ensure all labels are represented
for k in range(num_labels):
    if label_counts[k] < nmin:
        raise ValueError(
            f"Label {k} has less than minimum number of vertices, {nmin}, label_counts={label_counts}"
        )

logger.info(["Final label counts:", label_counts])


ap_srcsink = np.zeros((len(edges)), dtype=np.int32)
idx_edges = np.where(edges == 1)[0]
ap_srcsink[idx_edges[labels == 0]] = 1
ap_srcsink[idx_edges[labels == 1]] = 2

pd_srcsink = np.zeros((len(edges)), dtype=np.int32)
pd_srcsink[idx_edges[labels == 2]] = 1
pd_srcsink[idx_edges[labels == 3]] = 2


# save to gifti:
write_src_sink_label_gifti(
    ap_srcsink, snakemake.output.ap, structure_metadata=structure_metadata
)
write_src_sink_label_gifti(
    pd_srcsink, snakemake.output.pd, structure_metadata=structure_metadata
)
