import nibabel as nib
import numpy as np
from lib.myelin import divide_t1_by_t2

t1_nib = nib.load(snakemake.input.t1)
t2_nib = nib.load(snakemake.input.t2)

if t1_nib.shape != t2_nib.shape or not np.allclose(
    t1_nib.affine, t2_nib.affine, atol=1e-4
):
    raise ValueError(
        f"T1w and T2w are not on the same grid: {t1_nib.shape}, {t2_nib.shape}"
    )

ratio = divide_t1_by_t2(t1_nib.get_fdata(), t2_nib.get_fdata())

header = t1_nib.header.copy()
header.set_data_dtype(np.float32)
out_nib = nib.Nifti1Image(ratio.astype(np.float32), t1_nib.affine, header)
out_nib.header.set_slope_inter(1, 0)
out_nib.to_filename(snakemake.output.t1overt2)
