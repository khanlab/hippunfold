import os
import numpy as np
import nibabel as nib

in_warp = snakemake.input[0]
out_jac = snakemake.output[0]

warp_nib = nib.load(in_warp)
warp_data = warp_nib.get_fdata()
print(warp_data.shape)

#shape is (256, 128, 16, 1, 3)

#remove the singleton dimension to make easier to work with
warp_data = np.squeeze(warp_data)

# for u, v, w compute gradient
gx_u, gx_v, gx_w = np.gradient(warp_data[:,:,:,0])
gy_u, gy_v, gy_w = np.gradient(warp_data[:,:,:,1])
gz_u, gz_v, gz_w = np.gradient(warp_data[:,:,:,2])

#make shape vec for output
shape_jac = np.array(warp_data.shape)
shape_jac[3] = 9

jac_data = np.zeros(shape_jac)
jac_data[:,:,:,0] = gx_u
jac_data[:,:,:,1] = gy_u
jac_data[:,:,:,2] = gz_u
jac_data[:,:,:,3] = gx_v
jac_data[:,:,:,4] = gy_v
jac_data[:,:,:,5] = gz_v
jac_data[:,:,:,6] = gx_w
jac_data[:,:,:,7] = gy_w
jac_data[:,:,:,8] = gz_w

jac_nib = nib.Nifti1Image(jac_data, warp_nib.affine)

nib.save(jac_nib,out_jac)


