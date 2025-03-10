import nibabel as nib
import numpy as np
from scipy.spatial import Delaunay

# Load metric files
metric_files = [nib.load(metric).get_fdata() for metric in snakemake.input.metric_nii]
all_metrics = np.stack(metric_files, axis=2)

# Create a meshgrid spanning the metric space
x = np.linspace(0, all_metrics.shape[0] - 1, all_metrics.shape[0])
y = np.linspace(0, all_metrics.shape[1] - 1, all_metrics.shape[1])
coords_x, coords_y = np.meshgrid(y, x)
all_points = np.column_stack((coords_x.ravel(), coords_y.ravel()))

# Perform Delaunay triangulation
tri = Delaunay(all_points)

# Add z-values
z_level = snakemake.params.z_level
points = np.column_stack((all_points, np.full(all_points.shape[0], z_level)))

# Identify valid vertices (non-zero and finite values across all metrics)
valid_mask = ~np.all((all_metrics < 1e-6) | (all_metrics > 1e6), axis=2).ravel()
filtered_points = points[valid_mask]

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
gifti.to_filename(snakemake.output.surf)
