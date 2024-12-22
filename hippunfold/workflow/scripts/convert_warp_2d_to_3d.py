import numpy as np
import nibabel as nib

# Read 2D xfm nifti
xfm2d_nib = nib.load(snakemake.input.warp)
xfm2d_vol = xfm2d_nib.get_fdata()

# Read 3d ref nifti
ref3d_nib = nib.load(snakemake.input.ref)
ref3d_vol = ref3d_nib.get_fdata()


# Define the new shape
Nx, Ny, Nz = ref3d_vol.shape[:3]

if Nx != xfm2d_vol.shape[0] or Ny != xfm2d_vol.shape[1]:
    print(f'ref_vol: {ref3d_vol.shape}, warp_vol: {xfm2d_vol.shape}')
    raise ValueError('Ref nifti and warp nifti must have the same X and Y dimensions')


# Create a new array initialized with zeros
xfm3d_vol = np.zeros((Nx, Ny, Nz, 1, 3), dtype=xfm2d_vol.dtype)

# Insert the original array, which replicates in each zslice. Leaves the z-displacement untouched (as zero)
xfm3d_vol[...,:2] = xfm2d_vol

# Save as nifti
nib.Nifti1Image(xfm3d_vol,affine=ref3d_nib.affine,header=xfm2d_nib.header).to_filename(snakemake.output.warp)
