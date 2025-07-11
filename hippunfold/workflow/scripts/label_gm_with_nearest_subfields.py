import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree

# --- Load NIfTI files ---
binary_img = nib.load(snakemake.input.gm)
label_img = nib.load(snakemake.input.subfields)

binary_data = binary_img.get_fdata().astype(np.uint8)
label_data = label_img.get_fdata().astype(np.uint8)

# --- Get voxel coordinates ---
foreground_coords = np.argwhere(binary_data == 1)
label_coords = np.argwhere(label_data > 0)
label_values = label_data[label_data > 0]

# --- Build KD-Tree on label coords ---
tree = cKDTree(label_coords)

# --- Find nearest labels for each foreground voxel ---
_, nearest_idx = tree.query(foreground_coords)
nearest_labels = label_values[nearest_idx]

# --- Create output image ---
output_data = np.zeros_like(binary_data)
for i, coord in enumerate(foreground_coords):
    output_data[tuple(coord)] = nearest_labels[i]

# --- Save result ---
output_img = nib.Nifti1Image(
    output_data, affine=binary_img.affine, header=binary_img.header
)
nib.save(output_img, snakemake.output.subfields)
