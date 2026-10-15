import nibabel as nib
import numpy as np
from scipy.ndimage import distance_transform_edt

img = nib.load(snakemake.input.dseg)
data = np.asarray(img.dataobj).astype(np.int32)

fill_mask = np.isin(data, snakemake.params.fill_labels)
if fill_mask.any():
    # for every voxel, index of the nearest voxel that is not being filled
    nearest_idx = distance_transform_edt(
        fill_mask,
        sampling=img.header.get_zooms()[:3],
        return_distances=False,
        return_indices=True,
    )
    data = data[tuple(nearest_idx)]

out_img = nib.Nifti1Image(
    data.astype(img.get_data_dtype()), affine=img.affine, header=img.header
)
nib.save(out_img, snakemake.output.dseg)
