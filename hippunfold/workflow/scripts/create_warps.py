import nibabel as nib
import numpy as np
from  scipy.interpolate import griddata


def summary(name, array):
    print(f'{name}: shape={array.shape} mean={array.mean()} max={array.max()} min={array.min()}')
    return

unfold_ref_nib = nib.load(snakemake.input.unfold_ref_nii)
unfold_phys_coords_nib = nib.load(snakemake.input.unfold_phys_coords_nii)

coord_ap_nib = nib.load(snakemake.input.coord_ap)
coord_pd_nib = nib.load(snakemake.input.coord_pd)
coord_io_nib = nib.load(snakemake.input.coord_io)

coord_ap = coord_ap_nib.get_fdata()
coord_pd = coord_pd_nib.get_fdata()
coord_io = coord_io_nib.get_fdata()

mask = (coord_ap > 0) & (coord_pd > 0) & (coord_io > 0) # matlab: mask = (coord_ap>0 & coord_pd>0 & coord_io>0);
num_mask_voxels = np.sum(mask>0)
print(f'num_mask_voxels {num_mask_voxels}')

idxgm = np.flatnonzero(mask) #matlab: idxgm = find(mask ==1);
summary('idxgm',idxgm)

print(f'idxgm shape: {idxgm.shape}')
sz = mask.shape


coord_flat_ap = coord_ap[mask==True]  #matlab: Laplace_AP = coord_ap(mask==1);
coord_flat_pd = coord_pd[mask==True]
coord_flat_io = coord_io[mask==True]

summary('coord_flat_ap',coord_flat_ap)

(i_L,j_L,k_L) = np.unravel_index(idxgm,sz)  # matlab: [i_L,j_L,k_L]=ind2sub(sz,idxgm);

summary('i_L',i_L)
native_coords_mat = np.vstack((i_L,j_L,k_L,np.ones(i_L.shape))) #matlab: native_coords_mat = [i_L-1, j_L-1, k_L-1,ones(size(i_L))]';

summary('native_coords_mat',native_coords_mat)


#now, apply native image affine to get world coords

print(f'affine: {coord_ap_nib.affine}, affine shape: {coord_ap_nib.affine.shape}')


native_coords_phys = coord_ap_nib.affine @ native_coords_mat
native_coords_phys = np.transpose(native_coords_phys[:3,:])

summary('native_coords_phys',native_coords_phys)

#get unfolded grid:
unfold_grid_phys = unfold_phys_coords_nib.get_fdata() 

summary('unfold_grid_phys',unfold_grid_phys)


#griddata
# matlab:  interp_X = scatteredInterpolant(Laplace_AP,Laplace_PD,Laplace_IO,native_coords_phys(:,1),interp,extrap);

# we have points defined by coord_flat_{ap,pd,io}, and corresponding value as native_coords_phys[:,i]
# and we want to interpolate on a grid in the unfolded space

# to get the unfolded grid
points = (coord_flat_ap,coord_flat_pd,coord_flat_io)

(unfold_gx, unfold_gy, unfold_gz) = np.meshgrid(np.linspace(0,1,unfold_grid_phys.shape[0]),
                                        np.linspace(0,1,unfold_grid_phys.shape[1]),
                                        np.linspace(0,1,unfold_grid_phys.shape[2]),indexing='ij')

unfold_xi = (unfold_gx, unfold_gy, unfold_gz)
summary('unfold_gx',unfold_gx)
summary('unfold_gy',unfold_gy)
summary('unfold_gz',unfold_gz)

interp_ap = griddata(points,
                        values=native_coords_phys[:,0],
                        xi=unfold_xi,
                        method='linear',
                        fill_value=0)
interp_pd = griddata(points,
                        values=native_coords_phys[:,1],
                        xi=unfold_xi,
                        method='linear',
                        fill_value=0)
interp_io = griddata(points,
                        values=native_coords_phys[:,2],
                        xi=unfold_xi,
                        method='linear',
                        fill_value=0)

summary('interp_ap',interp_ap)

#write interp_ap
nib.Nifti1Image(interp_ap,unfold_ref_nib.affine, unfold_ref_nib.header).to_filename('test_interp_ap.nii.gz')


#combine and reshape interpolated map to 5d (4th dim singleton)
mapToNative = np.zeros(unfold_grid_phys.shape)
mapToNative[:,:,:,0,0] = interp_ap
mapToNative[:,:,:,0,1] = interp_pd
mapToNative[:,:,:,0,2] = interp_io

summary('mapToNative',mapToNative)

displacementToNative = mapToNative - unfold_grid_phys

summary('dispacementToNative',displacementToNative)

#mapToNative saved as abswarp_unfold2native.nii #not used by snakemake

#displacementToNative saved as Warp_unfold2native.nii, snakemake.output.warp_unfold2native
warp_unfold2native_nib = nib.Nifti1Image(displacementToNative,
                                        unfold_phys_coords_nib.affine,
                                        unfold_phys_coords_nib.header)
warp_unfold2native_nib.to_filename(snakemake.output.warp_unfold2native)
