import numpy as np
import nibabel as nib

warp = nib.load('reg1InverseWarp.nii.gz').get_fdata()
surf = nib.load('sub-01_hemi-L_space-unfolded_den-unfoldiso_label-hipp_inner.surf.gii')
vertices = surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data

# apply warp (2D)
warp = np.transpose(warp, (1,0,2,3,4))
vertices[:,0] = vertices[:,0] + warp[:,:,:,:,0].flatten()
vertices[:,1] = vertices[:,1] + warp[:,:,:,:,1].flatten()

surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data = vertices
nib.save(surf,'alignedSurf.surf.gii')
