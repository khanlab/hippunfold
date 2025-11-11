#!/usr/bin/env python3
import nibabel as nib
import numpy as np
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

logger.info("Remapping synthesized layer segmentation...")
# --- Input / Output ---
in_nii = snakemake.input[0]
out_nii = snakemake.output[0]
logger.info(f"Input: {in_nii}")
logger.info(f"Output: {out_nii}")

# --- Load ---
img = nib.load(in_nii)
data = img.get_fdata().astype(int)
logger.info(f"Data shape: {data.shape}, unique values: {np.unique(data)}")

# --- Initialize output as zeros ---
out = np.zeros_like(data, dtype=int)

# --- Remap values ---
mapping = {1: 1, 2: 1, 3: 1, 4: 1, 5: 8, 6: 8, 7: 2, 8: 7, 9: 5, 10: 6, 11: 2}

for k, v in mapping.items():
    out[data == k] = v

logger.info(f"Remapped unique values: {np.unique(out)}")
# --- Save ---
nib.save(nib.Nifti1Image(out, img.affine, img.header), out_nii)
print("Saved to:", out_nii)
