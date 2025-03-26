import nibabel as nib
import numpy as np
from scipy.spatial import Delaunay

# Load metric files
metric_ref = nib.load(snakemake.input.metric_ref).get_fdata()

# Create a meshgrid spanning the metric space
x = np.linspace(0, metric_ref.shape[0] - 1, metric_ref.shape[0])
y = np.linspace(0, metric_ref.shape[1] - 1, metric_ref.shape[1])


coords_y, coords_x = np.meshgrid(y, x)
all_points = np.column_stack((coords_x.ravel(), coords_y.ravel()))

# Perform Delaunay triangulation
tri = Delaunay(all_points)

# Add placeholder z-value (will update after transformation)
points = np.column_stack((all_points, np.full(all_points.shape[0], 0)))

# Identify valid vertices (non-zero and finite values across all metrics)
valid_mask = (metric_ref > 0).ravel()
filtered_points = points[valid_mask]


# apply transformation from vox to phys
affine = nib.load(snakemake.input.metric_ref).affine
coords = np.hstack((filtered_points, np.ones((filtered_points.shape[0], 1))))
transformed = affine @ coords.T
filtered_points = transformed.T[:, :3]

# now set z-level based on params
filtered_points[:, 2] = snakemake.params.z_level

# Remap indices for triangulation
new_indices = np.full(points.shape[0], -1, dtype=int)
new_indices[valid_mask] = np.arange(filtered_points.shape[0])
valid_faces_mask = np.all(valid_mask[tri.simplices], axis=1)
filtered_faces = new_indices[tri.simplices[valid_faces_mask]]

# Save as a Gifti file
points_darray = nib.gifti.GiftiDataArray(
    data=filtered_points.astype(np.float32),
    intent="NIFTI_INTENT_POINTSET",
    datatype="NIFTI_TYPE_FLOAT32",
)
tri_darray = nib.gifti.GiftiDataArray(
    data=filtered_faces.astype(np.int32),
    intent="NIFTI_INTENT_TRIANGLE",
    datatype="NIFTI_TYPE_INT32",
)
gifti = nib.GiftiImage()
gifti.add_gifti_data_array(points_darray)
gifti.add_gifti_data_array(tri_darray)
gifti.to_filename(snakemake.output.surf_gii)
