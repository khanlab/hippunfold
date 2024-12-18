import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation
from astropy.convolution import convolve

coords_nib = nib.load(snakemake.input.coords)
label_nib = nib.load(snakemake.input.labelmap)

coords = coords_nib.get_fdata()
label = label_nib.get_fdata()

#get the source and sink regions
src_mask = np.zeros(label.shape,dtype=bool)
for lbl in snakemake.params.src_labels:
    src_mask = np.where(label == lbl, True,src_mask)

sink_mask = np.zeros(label.shape,dtype=bool)
for lbl in snakemake.params.sink_labels:
    sink_mask = np.where(label == lbl, True,sink_mask)

gm_mask = np.zeros(label.shape,dtype=bool)
for lbl in snakemake.params.gm_labels:
    gm_mask = np.where(label == lbl, True,gm_mask)

gm_noDG_mask = np.zeros(label.shape,dtype=bool)
for lbl in snakemake.params.gm_noDG_labels:
    gm_noDG_mask = np.where(label == lbl, True,gm_noDG_mask)



#set everywhere outside coords mask to be nan, so they don't contribute to extrapolation
coords_withnan = np.where(gm_mask, coords, np.nan)

#perform convolution with astropy to fill in nans
hl = np.ones([3, 3, 3])
hl = hl / np.sum(hl)
convolved_coords = convolve(coords_withnan, hl, nan_treatment='interpolate', preserve_nan=False)

#dilating the mask into the src and sink
structuring_element = np.ones((3, 3, 3), dtype=bool)
dilated_mask = binary_dilation(gm_noDG_mask,structuring_element)  & (src_mask | sink_mask)
updated_mask = np.where(dilated_mask,1,gm_noDG_mask)

#get the updated coords:
updated_coords = np.where(gm_mask,coords,
                        np.where(dilated_mask, convolved_coords, 0) )


#save the extrapolated coords and the updated mask
nib.Nifti1Image(convolved_coords, coords_nib.affine, coords_nib.header).to_filename(snakemake.output.coords)
nib.Nifti1Image(updated_mask, coords_nib.affine, coords_nib.header).to_filename(snakemake.output.mask)

