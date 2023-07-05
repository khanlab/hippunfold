import numpy as np
import nibabel as nib

warp = nib.load(snakemake.input.invwarp[0])
field = warp.get_fdata()
surf = nib.load(snakemake.input.gii)
vertices = surf.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data

# put in same space
# v = np.concatenate((vertices,np.ones((len(vertices),1))),1)
# v = np.matmul(v,warp.affine)
v = vertices / warp.affine[0, 0]
v = v - np.min(v, axis=0)
v = v.astype(int)

# apply warp (2D)
for i in range(len(vertices)):
    vertices[i, 0] = vertices[i, 0] + field[v[i, 0], v[i, 1], 0, 0, 0]
    vertices[i, 1] = vertices[i, 1] + field[v[i, 0], v[i, 1], 0, 0, 1]

surf.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data = vertices
nib.save(surf, snakemake.output.gii)
