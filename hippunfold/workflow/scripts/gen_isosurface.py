import pyvista as pv
import nibabel as nib
import numpy as np
from vtk import vtkNIFTIImageReader


def write_surface_to_gifti(in_surface, out_surf_gii):

    faces = in_surface.faces
    faces = faces.reshape((int(faces.shape[0] / 4), 4))[:, 1:4]
    points = in_surface.points

    points_darray = nib.gifti.GiftiDataArray(
        data=points, intent="NIFTI_INTENT_POINTSET", datatype="NIFTI_TYPE_FLOAT32"
    )

    tri_darray = nib.gifti.GiftiDataArray(
        data=faces, intent="NIFTI_INTENT_TRIANGLE", datatype="NIFTI_TYPE_INT32"
    )

    gifti = nib.GiftiImage()
    gifti.add_gifti_data_array(points_darray)
    gifti.add_gifti_data_array(tri_darray)

    gifti.to_filename(out_surf_gii)


def remove_nan_points_faces(vertices, faces):
    # Step 1: Identify valid vertices (no NaN values)
    nan_mask = np.isnan(vertices).any(axis=1)
    valid_vertices = ~nan_mask  # True for valid rows
    valid_indices = np.where(valid_vertices)[0]

    # Step 2: Create a mapping from old to new indices
    new_indices_map = -np.ones(
        vertices.shape[0], dtype=int
    )  # Default to -1 for invalid vertices
    new_indices_map[valid_indices] = np.arange(len(valid_indices))

    # Step 3: Update the faces array to remove references to invalid vertices
    # Replace old indices with new ones, and remove faces with invalid vertices
    new_faces = []
    for face in faces:
        # Map old indices to new ones
        mapped_face = new_indices_map[face]
        if np.all(mapped_face >= 0):  # Include only faces with all valid vertices
            new_faces.append(mapped_face)

    new_faces = np.array(new_faces)

    # Step 4: Remove invalid vertices from the array
    new_vertices = vertices[valid_vertices]

    return (new_vertices, new_faces)


from scipy.ndimage import binary_dilation


def get_adjacent_voxels(mask_a, mask_b):
    """
    Create a mask for voxels where label A is adjacent to label B.

    Parameters:
    - mask_a (np.ndarray): A 3D binary mask for label A.
    - mask_b (np.ndarray): A 3D binary mask for label B.

    Returns:
    - np.ndarray: A 3D mask where adjacent voxels for label A and label B are marked as True.
    """
    # Dilate each mask to identify neighboring regions
    dilated_a = binary_dilation(mask_a)
    dilated_b = binary_dilation(mask_b)

    # Find adjacency: voxels of A touching B and B touching A
    adjacency_mask = (dilated_a.astype("bool") & mask_b.astype("bool")) | (
        dilated_b.astype("bool") & mask_a.astype("bool")
    )

    return adjacency_mask


# Load the coords image
coords_img = nib.load(snakemake.input.coords)
coords = coords_img.get_fdata()


# Load the nan mask
nan_mask_img = nib.load(snakemake.input.nan_mask)
nan_mask = nan_mask_img.get_fdata()

# Load the sink mask
sink_mask_img = nib.load(snakemake.input.sink_mask)
sink_mask = sink_mask_img.get_fdata()

# Load the src mask
src_mask_img = nib.load(snakemake.input.src_mask)
src_mask = src_mask_img.get_fdata()


affine = coords_img.affine

# Get voxel spacing from the header
voxel_spacing = coords_img.header.get_zooms()[:3]

# Create a PyVista grid
grid = pv.ImageData()
grid.dimensions = np.array(coords.shape) + 1  # Add 1 to include the boundary voxels
grid.spacing = (1, 1, 1)  # Use unit spacing and zero origin since we will apply affine
grid.origin = (0, 0, 0)

# update the coords data to add the nans and sink
coords[nan_mask == 1] = np.nan
coords[sink_mask == 1] = 1.1  # since sink being zero creates a false boundary

# we also need to use a nan mask for the voxels where src and sink meet directly
# (since this is another false boundary)..
src_sink_nan_mask = get_adjacent_voxels(sink_mask, src_mask)
coords[src_sink_nan_mask == 1] = np.nan


# Add the scalar field
grid.cell_data["values"] = coords.flatten(order="F")
grid = grid.cells_to_points("values")

# apply the affine
tfm_grid = grid.transform(
    affine, inplace=False
)  # Apply the rotation part of the affine


# the contour function produces the isosurface

surface = tfm_grid.contour([snakemake.params.threshold], method="contour").decimate(
    0.9
)  # fill_holes(snakemake.params.max_hole_size)
# surface = tfm_grid.contour([snakemake.params.threshold],method='contour').fill_holes(snakemake.params.max_hole_size)

# surface = surface.decimate(float(snakemake.params.decimate_percent) / 100.0)

# faces from pyvista surface are formatted with number of verts each row
# reshape and remove the first col to get Nx3
faces = surface.faces
faces = faces.reshape((int(faces.shape[0] / 4), 4))[:, 1:4]

points = surface.points


# with nans in background we end up with nan vertices, we can remove
# these to end up with an open contour..
new_points, new_faces = remove_nan_points_faces(points, faces)


# Step 1: Prepare the PolyData
# PyVista expects faces in a flat array with the number of points in each face as the first value
faces_flat = np.hstack(
    [[3] + list(face) for face in new_faces]
)  # Add '3' for triangular faces

# Create a new PolyData object
polydata = pv.PolyData(new_points, faces_flat)

#  Ensure the PolyData is clean (optional)
# polydata.clean(inplace=True)  # Removes unused points, degenerate cells, etc.

#  Apply decimation (optional)
# polydata.decimate_pro(snakemake.params.decimate_percent / 100.0, inplace=True)

# write to gifti
write_surface_to_gifti(polydata, snakemake.output.surf_gii)
