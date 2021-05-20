import nibabel as nib
import numpy as np
from  scipy.interpolate import griddata

unfold_ref_nib = nib.load(snakemake.input.unfold_ref_nii)
unfold_phys_coords_nib = nib.load(snakemake.input.unfold_phys_coords_nii)

coord_ap_nib = nib.load(snakemake.input.coord_ap)
coord_pd_nib = nib.load(snakemake.input.coord_pd)
coord_io_nib = nib.load(snakemake.input.coord_io)

coord_ap = coord_ap_nib.get_fdata()
coord_pd = coord_pd_nib.get_fdata()
coord_io = coord_io_nib.get_fdata()

mask = (coord_ap > 0) & (coord_pd > 0) & (coord_io > 0) # matlab: mask = (coord_ap>0 & coord_pd>0 & coord_io>0);
idxgm = np.flatnonzero(mask) #matlab: idxgm = find(mask ==1);
print(f'idxgm shape: {idxgm.shape}')
sz = mask.shape


coord_flat_ap = coord_ap[mask==True]  #matlab: Laplace_AP = coord_ap(mask==1);
coord_flat_pd = coord_pd[mask==True]
coord_flat_io = coord_io[mask==True]


(i_L,j_L,k_L) = np.unravel_index(idxgm,sz)  # matlab: [i_L,j_L,k_L]=ind2sub(sz,idxgm);

print(f'i_L: shape: {i_L.shape}')
native_coords_mat = np.vstack((i_L,j_L,k_L,np.ones(i_L.shape))) #matlab: native_coords_mat = [i_L-1, j_L-1, k_L-1,ones(size(i_L))]';

print(f'native_coords_mat: shape: {native_coords_mat.shape}')
#now, apply native image affine to get world coords

print(f'affine: {coord_ap_nib.affine}, affine shape: {coord_ap_nib.affine.shape}')

native_coords_phys = coord_ap_nib.affine @ native_coords_mat
native_coords_phys = np.transpose(native_coords_phys[:3,:])

#griddata
#griddata(points, values, xi, method='linear', fill_value=nan, rescale=False)[source]¶
