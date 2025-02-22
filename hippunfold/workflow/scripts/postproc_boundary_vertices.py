import pyvista as pv
import numpy as np
import nibabel as nib
import nibabel.gifti as gifti
from collections import Counter
from scipy.signal import argrelextrema


def write_src_sink_label_gifti(label_data, out_label_gifti):

    # Create a Label Table (LUT)
    label_table = gifti.GiftiLabelTable()

    # Define Background label (key 0)
    background_label = gifti.GiftiLabel(
        key=0, red=1.0, green=1.0, blue=1.0, alpha=0.0
    )  # Transparent
    background_label.label = "Background"
    label_table.labels.append(background_label)

    # Define src label (key 1)
    src_label = gifti.GiftiLabel(
        key=1, red=1.0, green=0.0, blue=0.0, alpha=1.0
    )  # Red color
    src_label.label = "Source"
    label_table.labels.append(src_label)

    # Define sink label (key 2)
    sink_label = gifti.GiftiLabel(
        key=1, red=1.0, green=0.0, blue=0.0, alpha=1.0
    )  # Red color
    sink_label.label = "Sink"
    label_table.labels.append(sink_label)

    # write the data and the label table to gifti
    gii_img = gifti.GiftiImage(
        darrays=[gifti.GiftiDataArray(label_data, intent="NIFTI_INTENT_LABEL")],
        labeltable=label_table,
    )
    nib.save(gii_img, out_label_gifti)


nmin = snakemake.params.min_terminal_vertices

logger.info("Loading surface from GIFTI...")

edges = nib.load(snakemake.input.edges).agg_data()
ap_src = nib.load(snakemake.input.ap_src).agg_data()
ap_sink = nib.load(snakemake.input.ap_sink).agg_data()
pd_src = nib.load(snakemake.input.pd_src).agg_data()
pd_sink = nib.load(snakemake.input.pd_sink).agg_data()


logger.info("Assigning labels apsrc, apsink, pdsrc, pdsink")

distances = np.vstack(
    (ap_src[edges == 1], ap_sink[edges == 1], pd_src[edges == 1], pd_sink[edges == 1])
).T
scaling_factors = np.ones((4))
num_labels = 4  # starts at 0

max_iterations = 10  # Prevent infinite loops
for _ in range(max_iterations):
    # Scale distances
    scaled_distances = distances * scaling_factors

    # Assign labels based on min scaled distance
    labels = np.argmin(scaled_distances, axis=1)

    # Count occurrences per label
    unique, counts = np.unique(labels, return_counts=True)
    label_counts = dict(zip(unique, counts))
    # Ensure all labels in range(num_labels) are present, setting missing ones to 0
    label_counts = {k: label_counts.get(k, 0) for k in range(num_labels)}

    # Check if all labels meet nmin
    if all(count >= nmin for count in label_counts.values()):
        break  # Stop if all labels are sufficiently represented

    # Update scaling factors for underrepresented labels
    for k in range(num_labels):
        if label_counts[k] < nmin:
            scaling_factors[k] *= 0.9  # Increase competitiveness of the label

# Ensure all labels are represented
logger.info(["Final label counts:", label_counts])

if 0 in label_counts.values():
    raise ValueError("Encountered a src/sink label with zero vertices, {label_counts}")

ap_srcsink = np.zeros((len(edges)), dtype=np.int32)
idx_edges = np.where(edges == 1)[0]
ap_srcsink[idx_edges[labels == 0]] = 1
ap_srcsink[idx_edges[labels == 1]] = 2

pd_srcsink = np.zeros((len(edges)), dtype=np.int32)
pd_srcsink[idx_edges[labels == 2]] = 1
pd_srcsink[idx_edges[labels == 3]] = 2

# save to gifti:
write_src_sink_label_gifti(ap_srcsink, snakemake.output.ap)
write_src_sink_label_gifti(pd_srcsink, snakemake.output.pd)
