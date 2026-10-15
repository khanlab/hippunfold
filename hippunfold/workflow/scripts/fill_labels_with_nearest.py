import nibabel as nib
import numpy as np
from scipy.ndimage import binary_closing, distance_transform_edt

img = nib.load(snakemake.input.dseg)
data = np.asarray(img.dataobj).astype(np.int32)
voxel_size = np.array(img.header.get_zooms()[:3], dtype=float)


def fill_with_nearest(data, fill_mask, source_mask):
    """Replace the voxels in fill_mask with the label of the nearest voxel in source_mask"""
    nearest_idx = distance_transform_edt(
        ~source_mask,
        sampling=voxel_size,
        return_distances=False,
        return_indices=True,
    )
    data[fill_mask] = data[tuple(nearest_idx)][fill_mask]


def ellipsoid(radius_mm):
    """Structuring element with the given radius in mm, so its size in voxels scales with
    the resolution of the image (and per axis for anisotropic voxels)"""
    # small tolerance, as voxel sizes are stored as float32 (e.g. 0.6mm / 0.3mm -> 1.9999999)
    radius_vox = radius_mm / voxel_size + 1e-4
    half_width = np.floor(radius_vox).astype(int)
    grid = np.ogrid[tuple(slice(-w, w + 1) for w in half_width)]
    return sum((g / r) ** 2 for g, r in zip(grid, radius_vox)) <= 1


bg_mask = np.isin(data, snakemake.params.bg_labels)
data[bg_mask] = 0

fill_mask = np.isin(data, snakemake.params.fill_labels)
if fill_mask.any():
    fill_with_nearest(
        data, fill_mask, ~fill_mask & ~np.isin(data, snakemake.params.exclude_labels)
    )

# background voxels that fill a concavity of the tissue take the nearest tissue label
if bg_mask.any() and snakemake.params.closing_radius_mm > 0:
    tissue_mask = np.isin(data, snakemake.params.tissue_labels)
    closed_mask = binary_closing(
        tissue_mask, structure=ellipsoid(snakemake.params.closing_radius_mm)
    )
    concavity_mask = bg_mask & closed_mask
    if concavity_mask.any():
        fill_with_nearest(data, concavity_mask, tissue_mask)

out_img = nib.Nifti1Image(
    data.astype(img.get_data_dtype()), affine=img.affine, header=img.header
)
nib.save(out_img, snakemake.output.dseg)
