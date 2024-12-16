import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation
from astropy.convolution import convolve

coords_nib = nib.load(snakemake.input.coords)
label_nib = nib.load(snakemake.input.labelmap)

coords = coords_nib.get_fdata()
label = label_nib.get_fdata()

gm_mask = np.zeros(label.shape,dtype=bool)
for lbl in snakemake.params.gm_labels:
    gm_mask = np.where(label == lbl, True,gm_mask)



#set everywhere outside coords mask to be nan, so they don't contribute to extrapolation
coords_withnan = np.where(gm_mask, coords, np.nan)

#perform convolution with astropy to fill in nans
hl = np.ones([3, 3, 3])
hl = hl / np.sum(hl)
convolved_coords = convolve(coords_withnan, hl, nan_treatment='interpolate', preserve_nan=False)

#dilating the mask into the src and sink
structuring_element = np.ones((3, 3, 3), dtype=bool)
dilated_mask = binary_dilation(gm_mask,structuring_element) 

#for dentate, we use the convolved coords (hopefully the smoothing is helpful here for making a surface

#save the extrapolated coords and the updated mask
nib.Nifti1Image(convolved_coords, coords_nib.affine, coords_nib.header).to_filename(snakemake.output.coords)
nib.Nifti1Image(dilated_mask, coords_nib.affine, coords_nib.header).to_filename(snakemake.output.mask)

