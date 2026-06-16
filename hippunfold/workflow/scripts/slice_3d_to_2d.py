import nibabel as nib
import numpy as np

img = nib.load(snakemake.input.img)
img_data = img.get_fdata()[:, :, 0]
nib.Nifti1Image(
    np.clip(img_data, snakemake.params.clip_min, snakemake.params.clip_max), img.affine
).to_filename(snakemake.output.img)
