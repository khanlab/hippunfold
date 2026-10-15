import nibabel as nib
import numpy as np
from scipy.ndimage import (
    convolve,
    distance_transform_edt,
    generate_binary_structure,
    label,
)

# --- Load NIfTI files ---
binary_img = nib.load(snakemake.input.gm)
label_img = nib.load(snakemake.input.subfields)

gm = binary_img.get_fdata() == 1
labels = np.asarray(label_img.dataobj).astype(np.int32)

# --- Crop to the region of interest for speed ---
bbox = tuple(
    slice(max(idx.min() - 1, 0), idx.max() + 2) for idx in np.nonzero(gm | (labels > 0))
)
gm_crop = gm[bbox]
labels_crop = labels[bbox]
label_values = np.unique(labels_crop[labels_crop > 0])

# --- Keep only the largest connected piece of each subfield ---
# (surface to volume mapping can leave small disconnected pieces where the surface folds
# back on itself; these are unlabelled here and refilled from their neighbours below)
connectivity = generate_binary_structure(3, 3)
for value in label_values:
    pieces, n_pieces = label(labels_crop == value, structure=connectivity)
    if n_pieces > 1:
        sizes = np.bincount(pieces.ravel())
        sizes[0] = 0
        labels_crop[(pieces > 0) & (pieces != sizes.argmax())] = 0

# --- Grow labels into unlabelled GM one voxel layer at a time ---
# (labels only spread through face-connected GM, so they can't jump across a fold or
# the SRLM to whichever voxel is nearest in straight-line distance)
neighbours = generate_binary_structure(3, 1).astype(np.uint8)
neighbours[1, 1, 1] = 0

unlabelled = gm_crop & (labels_crop == 0)
while unlabelled.any() and label_values.size > 0:
    counts = np.stack(
        [
            convolve(
                (labels_crop == value).astype(np.uint8), neighbours, mode="constant"
            )
            for value in label_values
        ]
    )
    reached = unlabelled & (counts.max(axis=0) > 0)
    if not reached.any():
        break
    labels_crop[reached] = label_values[counts.argmax(axis=0)][reached]
    unlabelled &= ~reached

# --- GM not connected to any label: fall back to straight-line nearest label ---
if unlabelled.any() and label_values.size > 0:
    nearest_idx = distance_transform_edt(
        labels_crop == 0, return_distances=False, return_indices=True
    )
    labels_crop[unlabelled] = labels_crop[tuple(nearest_idx)][unlabelled]

# --- Create output image (GM voxels only) ---
output_data = np.zeros(gm.shape, dtype=np.uint8)
output_data[bbox] = np.where(gm_crop, labels_crop, 0)

# --- Save result ---
output_img = nib.Nifti1Image(
    output_data, affine=binary_img.affine, header=binary_img.header
)
nib.save(output_img, snakemake.output.subfields)
