import numpy as np
import nibabel as nib
from scipy.stats import mode

# Load input NIfTI files
nii_files = snakemake.input
nii_images = [nib.load(f) for f in nii_files]
nii_data = [img.get_fdata() for img in nii_images]

# Ensure all images have the same shape
shape = nii_data[0].shape
if not all(img.shape == shape for img in nii_data):
    raise ValueError("All input NIfTI images must have the same dimensions.")

# Stack along a new axis and compute the majority vote
nii_stack = np.stack(nii_data, axis=-1)
majority_vote, _ = mode(nii_stack, axis=-1, keepdims=False)

# Save the result as a new NIfTI file
out_img = nib.Nifti1Image(
    majority_vote.squeeze(), nii_images[0].affine, nii_images[0].header
)
nib.save(out_img, snakemake.output[0])
