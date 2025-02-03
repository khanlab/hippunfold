import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure

# Load NIfTI images from Snakemake inputs
gm_nii = nib.load(snakemake.input.gm_mask)
pd_nii = nib.load(snakemake.input.srcsinkPD)
io_nii = nib.load(snakemake.input.srcsinkIO)
ap_nii = nib.load(snakemake.input.srcsinkAP)

# Extract data arrays
gm_data = gm_nii.get_fdata().astype(int)
pd_data = pd_nii.get_fdata().astype(int)
io_data = io_nii.get_fdata().astype(int)
ap_data = ap_nii.get_fdata().astype(int)

# Create an empty composite label map
composite = np.zeros_like(gm_data)

# Step 1: Preserve GM label but remap to 8
composite[gm_data == 1] = 8  # GM now has intensity 8

# Step 2: Dilate GM to get surrounding region (5x5x5 kernel)
struct = generate_binary_structure(3, 3)  # 3D connectivity
dilated_gm = binary_dilation(
    gm_data == 1, structure=struct, iterations=2
)  # Larger dilation

# Get only the new voxels from dilation (i.e., GM surroundings)
boundary_mask = dilated_gm & (gm_data == 0)

# Step 3: Apply priority labeling for the boundary region (with dilation)


def dilate_labels(label_data, label_values, iterations=1):
    """Dilate specific labels in a binary mask."""
    mask = np.isin(label_data, label_values)
    return binary_dilation(mask, structure=struct, iterations=iterations)


priority_labels = np.zeros_like(composite)

# First label (will get overwritten by other overlapping): AP src/sink
ap_src_mask = dilate_labels(ap_data, [5, 6])
priority_labels[ap_src_mask] = ap_data[ap_src_mask]


# Priority 1: IO src/sink
io_src_mask = dilate_labels(io_data, [3, 4])
priority_labels[io_src_mask] = io_data[io_src_mask]

# Priority 3: PD src/sink
pd_src_mask = dilate_labels(pd_data, [1, 2])
priority_labels[pd_src_mask] = pd_data[pd_src_mask]

# Step 4: Assign boundary voxels with priority
composite[boundary_mask] = priority_labels[boundary_mask]

# Step 5: Assign remaining unlabelled boundary voxels to 100 to mark for later
composite[(boundary_mask) & (composite == 0)] = 100

# Save composite as NIfTI using Snakemake output path
composite_nii = nib.Nifti1Image(composite, affine=gm_nii.affine, header=gm_nii.header)
nib.save(composite_nii, snakemake.output.dseg)
