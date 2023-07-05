import nibabel as nib
import numpy as np

warp2d = nib.load(snakemake.input.warp2d)
unfold_phys_coords_nib = nib.load(snakemake.input.unfold_phys_coords_nii)

IOdim = unfold_phys_coords_nib.get_fdata().shape[2]

field = warp2d.get_fdata()
shp = np.array(field.shape)
field_z = np.concatenate((field, np.zeros((np.append(shp[:4], 1)))), axis=4)
field_3d = np.repeat(field_z, IOdim, 2)

warp3d = nib.Nifti1Image(
    field_3d, unfold_phys_coords_nib.affine, unfold_phys_coords_nib.header
)
warp3d.to_filename(snakemake.output.warp3d)
