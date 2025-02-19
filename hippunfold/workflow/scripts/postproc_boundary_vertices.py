import pyvista as pv
import numpy as np
import nibabel as nib
import nibabel.gifti as gifti
from collections import Counter
from scipy.signal import argrelextrema

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
    label_counts = {
        k: counts[i] if k in unique else 0 for i, k in enumerate(range(num_labels))
    }

    # Check if all labels meet nmin
    if all(count >= nmin for count in label_counts.values()):
        break  # Stop if all labels are sufficiently represented

    # Update scaling factors for underrepresented labels
    for k in range(num_labels):
        if label_counts[k] < nmin:
            scaling_factors[k] *= 0.9  # Increase competitiveness of the label

# Ensure all labels are represented
logger.info(["Final label counts:", label_counts])

ap_srcsink = np.zeros((len(edges)), dtype=np.int32)
idx_edges = np.where(edges == 1)[0]
ap_srcsink[idx_edges[labels == 0]] = 1
ap_srcsink[idx_edges[labels == 1]] = 2
gii_img = gifti.GiftiImage(
    darrays=[gifti.GiftiDataArray(ap_srcsink, intent="NIFTI_INTENT_LABEL")]
)
nib.save(gii_img, snakemake.output.ap)

pd_srcsink = np.zeros((len(edges)), dtype=np.int32)
pd_srcsink[idx_edges[labels == 2]] = 1
pd_srcsink[idx_edges[labels == 3]] = 2
gii_img = gifti.GiftiImage(
    darrays=[gifti.GiftiDataArray(pd_srcsink, intent="NIFTI_INTENT_LABEL")]
)
nib.save(gii_img, snakemake.output.pd)
