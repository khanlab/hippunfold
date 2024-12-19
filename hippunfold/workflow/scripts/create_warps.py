import nibabel as nib
import numpy as np
import naturalneighbor
from scipy.stats import zscore
from scipy.ndimage import generic_filter, binary_dilation
from astropy.convolution import convolve

logfile = open(snakemake.log[0], "w")
print(f"Start", file=logfile, flush=True)


def extrapolate_with_convolve(input_array, mask_array):
    """
    Perform convolution to replace nearby NaN values with interpolated values,
    while preserving values at all non-NaN voxels.

    Parameters:
    - input_array: ndarray
        The input array with possible NaN values to be replaced.
    - mask_array: ndarray
        A boolean mask where True indicates regions to preserve in the input array.

    Returns:
    - output_array: ndarray
        The array with NaN values replaced using interpolation from the neighborhood.
    """
    # 1. Create working array: set to NaN outside mask, else input_array
    working_array = np.where(mask_array, input_array, np.nan)

    # 2. Define a convolution kernel (3x3x3 Box Kernel)
    # set up filter (18NN)
    hl = np.ones([5, 5, 5])
    hl = hl / np.sum(hl)

    # 3. Perform convolution on the working array with NaNs
    convolved_array = convolve(
        working_array, hl, nan_treatment="interpolate", preserve_nan=False
    )

    # 4. Restore original values inside the mask
    working_array = np.where(mask_array, input_array, convolved_array)

    return working_array


def convert_warp_to_itk(warp):
    """Convert warp to ITK convention by negating X and Y"""
    warp_itk = warp.copy()
    warp_itk[:, :, :, 0, 0] = -warp_itk[:, :, :, 0, 0]
    warp_itk[:, :, :, 0, 1] = -warp_itk[:, :, :, 0, 1]
    return warp_itk


def summary(name, array):
    """simple function to print stats on an array"""
    print(
        f"{name}: shape={array.shape} mean={np.nanmean(array)} max={np.nanmax(array)} min={np.nanmin(array)}, numNaNs={np.count_nonzero(np.isnan(array))}, type={array.dtype.name}",
        file=logfile,
        flush=True,
    )
    return


# load unfolded coordinate map
# unfold_ref_nib = nib.load(snakemake.input.unfold_ref_nii)
unfold_phys_coords_nib = nib.load(snakemake.input.unfold_phys_coords_nii)
native_ref_coords_nib = nib.load(snakemake.input.native_ref_coords_nii)

# load laplace coords
coord_ap_nib = nib.load(snakemake.input.coords_ap)
coord_pd_nib = nib.load(snakemake.input.coords_pd)
coord_io_nib = nib.load(snakemake.input.coords_io)

coord_ap = coord_ap_nib.get_fdata()
coord_pd = coord_pd_nib.get_fdata()
coord_io = coord_io_nib.get_fdata()

# get mask of coords
lbl_nib = nib.load(snakemake.input.labelmap)
lbl = lbl_nib.get_fdata()
idxgm = np.zeros(lbl.shape)
for i in snakemake.params.gm_labels:
    idxgm[lbl == i] = 1
mask = idxgm == 1

# extrapolate coords to widen domain
coord_ap = extrapolate_with_convolve(coord_ap, mask)
coord_io = extrapolate_with_convolve(coord_io, mask)
coord_pd = extrapolate_with_convolve(coord_pd, mask)


# then, dilate the mask, to get a larger domain for the warp
structuring_element = np.ones((3, 3, 3), dtype=bool)
mask = binary_dilation(mask, structuring_element)
mask = binary_dilation(mask, structuring_element)


num_mask_voxels = np.sum(mask)
print(f"num_mask_voxels {num_mask_voxels}", file=logfile, flush=True)

# get indices of mask voxels
idxgm = np.flatnonzero(mask)  # matlab: idxgm = find(mask ==1);
summary("idxgm", idxgm)
sz = mask.shape

# Part 1: unfold2native warps

# For this, we want to define a mapping of each grid point in unfolded space
# to points in the native space. We have mappings from native to unfold from
# the laplace solution, but these map to a set of unstructured points in the
# unfolded space. So we use these points and interpolate on a regular grid in
# the unfolded space, using scipy's griddata (equivalent to matlab scatteredInterpolant).


coord_flat_ap = coord_ap[mask == True]
coord_flat_pd = coord_pd[mask == True]
coord_flat_io = coord_io[mask == True]

summary("coord_flat_ap", coord_flat_ap)
summary("coord_flat_pd", coord_flat_pd)
summary("coord_flat_io", coord_flat_io)

# unravel indices of mask voxels into subscripts...
(i_L, j_L, k_L) = np.unravel_index(idxgm, sz)
summary("i_L", i_L)

# ... and stack into vectors ...
native_coords_mat = np.vstack((i_L, j_L, k_L, np.ones(i_L.shape)))
summary("native_coords_mat", native_coords_mat)


# ... then,apply native image affine to get world coords ...
aff = coord_ap_nib.affine
print(f"affine: {aff}, affine shape: {aff.shape}", file=logfile, flush=True)
native_coords_phys = aff @ native_coords_mat
native_coords_phys = np.transpose(native_coords_phys[:3, :])
summary("native_coords_phys", native_coords_phys)

# get unfolded grid from file:
unfold_grid_phys = unfold_phys_coords_nib.get_fdata()

# unfolded space dims
unfold_dims = unfold_grid_phys.shape[:3]
summary("unfold_grid_phys", unfold_grid_phys)


# scattered interpolation / griddata:

# we have points defined by coord_flat_{ap,pd,io}, and corresponding value as native_coords_phys[:,i]
# and we want to interpolate on a grid in the unfolded space

# add some noise to avoid perfectly overlapping datapoints!
points = (
    coord_flat_ap * unfold_dims[0]
    + (np.random.rand(coord_flat_ap.shape[0]) - 0.5) * 1e-6,
    coord_flat_pd * unfold_dims[1]
    + (np.random.rand(coord_flat_ap.shape[0]) - 0.5) * 1e-6,
    coord_flat_io * unfold_dims[2]
    + (np.random.rand(coord_flat_ap.shape[0]) - 0.5) * 1e-6,
)

# get unfolded grid (from 0 to 1, not world coords), using meshgrid:
#  note: indexing='ij' to swap the ordering of x and y
epsilon = snakemake.params.epsilon
(unfold_gx, unfold_gy, unfold_gz) = np.meshgrid(
    np.linspace(
        0 + float(epsilon[0]), unfold_dims[0] - float(epsilon[0]), unfold_dims[0]
    ),
    np.linspace(
        0 + float(epsilon[1]), unfold_dims[1] - float(epsilon[1]), unfold_dims[1]
    ),
    np.linspace(
        0 + float(epsilon[2]), unfold_dims[2] - float(epsilon[2]), unfold_dims[2]
    ),
    indexing="ij",
)
summary("points", points)


# perform the interpolation

interp_ap = naturalneighbor.griddata(
    points,
    native_coords_phys[:, 0],
    [[0, unfold_dims[0], 1], [0, unfold_dims[1], 1], [0, unfold_dims[2], 1]],
)
interp_pd = naturalneighbor.griddata(
    points,
    native_coords_phys[:, 1],
    [[0, unfold_dims[0], 1], [0, unfold_dims[1], 1], [0, unfold_dims[2], 1]],
)
interp_io = naturalneighbor.griddata(
    points,
    native_coords_phys[:, 2],
    [[0, unfold_dims[0], 1], [0, unfold_dims[1], 1], [0, unfold_dims[2], 1]],
)


# prepare maps for writing as warp file:

# combine and reshape interpolated map to 5d (4th dim singleton)
mapToNative = np.zeros(unfold_grid_phys.shape)
mapToNative[:, :, :, 0, 0] = interp_ap
mapToNative[:, :, :, 0, 1] = interp_pd
mapToNative[:, :, :, 0, 2] = interp_io
summary("mapToNative", mapToNative)

mapToNative[np.isnan(mapToNative)] = 0


# mapToNative has the absolute coordinates, but we want them relative to the
# unfolded grid, so we subtract it out:
displacementToNative = mapToNative - unfold_grid_phys
summary("dispacementToNative", displacementToNative)


# write to file
dt = unfold_phys_coords_nib.get_fdata()
dt = dt.dtype.name
warp_unfold2native_nib = nib.Nifti1Image(
    displacementToNative.astype(dt),
    unfold_phys_coords_nib.affine,
    unfold_phys_coords_nib.header,
)
warp_unfold2native_nib.to_filename(snakemake.output.warp_unfold2native)

# write itk transform to file
f = convert_warp_to_itk(displacementToNative)
warpitk_native2unfold_nib = nib.Nifti1Image(
    f.astype(dt), unfold_phys_coords_nib.affine, unfold_phys_coords_nib.header
)

warpitk_native2unfold_nib.to_filename(snakemake.output.warpitk_native2unfold)


# Part 2: native2unfold warps

# This part is easier, since we essentially have the mapping from native to unfold
# in the Laplace solution, so we use these to construct displacement maps (by dealing
# with world coordinates properly).

# First we want to convert the Laplace coordinate from 0-1 to world coordinates.
# The image affine from the unfolded grid takes points from 0 to N to world coords, so
# just need to un-normalize, then multiply by affine


# reshape for multiplication (affine * vec)
coord_flat_unnorm_vec = np.stack(
    (
        coord_flat_ap_unnorm,
        coord_flat_pd_unnorm,
        coord_flat_io_unnorm,
        np.ones(coord_flat_ap_unnorm.shape),
    ),
    axis=-1,
)

summary("coord_flat_unnorm_vec", coord_flat_unnorm_vec)

# multiply by affine
uvw1_phys = np.transpose(
    unfold_phys_coords_nib.affine @ np.transpose(coord_flat_unnorm_vec)
)
summary("uvw1_phys", uvw1_phys)

uvw_phys = uvw1_phys[:, :3]
summary("uvw_phys", uvw_phys)

# now we have the absolute unfold world coords for each native grid point
# but we need the displacement from the native grid point world coord

# the native world coords are in native_coords_phys
# so we subtract it
displace_to_unfold_vec = uvw_phys - native_coords_phys

summary("displace_to_unfold_vec", displace_to_unfold_vec)

# create new shape as 5d vector image in native space
native_map_shape = np.ones(
    5,
)
native_map_shape[:3] = mask.shape
native_map_shape[-1] = 3


# create and populate it
native_to_unfold = np.zeros(native_map_shape.astype("int"))
for d in range(3):
    temp_coords_img = np.zeros(mask.shape)
    temp_coords_img[mask == True] = displace_to_unfold_vec[:, d]
    native_to_unfold[:, :, :, 0, d] = temp_coords_img

summary("native_to_unfold", native_to_unfold)

# now, can write it to file
dt = coord_ap_nib.get_fdata()
dt = dt.dtype.name
warp_native2unfold_nib = nib.Nifti1Image(
    native_to_unfold.astype(dt), coord_ap_nib.affine, coord_ap_nib.header
)
warp_native2unfold_nib.to_filename(snakemake.output.warp_native2unfold)

# and save ITK warp too
f = convert_warp_to_itk(native_to_unfold)
warpitk_unfold2native_nib = nib.Nifti1Image(
    f.astype(dt), native_ref_coords_nib.affine, native_ref_coords_nib.header
)
warpitk_unfold2native_nib.to_filename(snakemake.output.warpitk_unfold2native)
