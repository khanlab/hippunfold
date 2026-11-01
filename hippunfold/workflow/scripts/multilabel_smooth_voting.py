"""Multi-label smoothing fusion of subfield labels in the average unfolded space.

For each subject: split the label map into one binary channel per label,
gaussian smooth each channel, then warp it to the average space with linear
interpolation. The warped channels are averaged across subjects and each voxel
is assigned the label with the highest average (similar to c3d -vote after
-split / -smooth / -reslice).
"""

import os
import subprocess
import tempfile

import nibabel as nib
import numpy as np
from scipy.ndimage import gaussian_filter

img_files = snakemake.input.imgs
warp_files = snakemake.input.warps
sigma_mm = float(snakemake.params.sigma_mm)

if len(img_files) != len(warp_files):
    raise ValueError("Number of subfield images and warps must match.")

ref_nii = nib.load(img_files[0])
ref_shape = ref_nii.shape

# union of labels across all subjects (including background)
labels = np.unique(
    np.concatenate(
        [np.unique(np.asarray(nib.load(f).dataobj).astype(np.int32)) for f in img_files]
    )
)
print(f"labels: {labels.tolist()}")

accum = np.zeros(ref_shape + (len(labels),), dtype=np.float64)

with tempfile.TemporaryDirectory() as tmpdir:
    in_chan = os.path.join(tmpdir, "chan.nii.gz")
    out_chan = os.path.join(tmpdir, "chan_warped.nii.gz")

    for img_file, warp_file in zip(img_files, warp_files):
        nii = nib.load(img_file)
        if nii.shape != ref_shape:
            raise ValueError(f"{img_file} has shape {nii.shape}, expected {ref_shape}")
        data = np.asarray(nii.dataobj).astype(np.int32)

        # sigma in voxels, no smoothing along singleton dims
        zooms = nii.header.get_zooms()[: data.ndim]
        sigma_vox = [0.0 if n == 1 else sigma_mm / z for n, z in zip(data.shape, zooms)]

        hdr = nii.header.copy()
        hdr.set_data_dtype(np.float32)

        for i, label in enumerate(labels):
            chan = gaussian_filter(
                (data == label).astype(np.float32), sigma=sigma_vox, mode="reflect"
            )
            nib.Nifti1Image(chan, nii.affine, hdr).to_filename(in_chan)

            subprocess.run(
                [
                    "antsApplyTransforms",
                    "-d",
                    "2",
                    "-i",
                    in_chan,
                    "-o",
                    out_chan,
                    "-t",
                    warp_file,
                    "-r",
                    in_chan,
                    "-n",
                    "Linear",
                ],
                check=True,
            )
            accum[..., i] += np.asarray(nib.load(out_chan).dataobj).reshape(ref_shape)

voted = labels[np.argmax(accum, axis=-1)]

out_hdr = ref_nii.header.copy()
out_hdr.set_data_dtype(np.int16)
out_img = nib.Nifti1Image(voted.astype(np.int16).squeeze(), ref_nii.affine, out_hdr)
nib.save(out_img, snakemake.output[0])
