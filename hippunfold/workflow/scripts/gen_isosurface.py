import pyvista as pv
import nibabel as nib
import numpy as np
from copy import deepcopy
import pyvista as pv
import pygeodesic.geodesic as geodesic


def write_surface_to_gifti(points, faces, out_surf_gii):

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


import numpy as np


def apply_affine_transform(vertices, affine_matrix, inverse=False):
    """
    Apply an affine transformation to a 3D mesh.

    Parameters:
    - vertices (np.ndarray): (N, 3) array of vertex coordinates.
    - affine_matrix (np.ndarray): (4, 4) affine transformation matrix.
    - inverse (bool): If True, applies the inverse transformation.

    Returns:
    - transformed_vertices (np.ndarray): (N, 3) array of transformed vertex coordinates.
    """
    if inverse:
        affine_matrix = np.linalg.inv(affine_matrix)

    # Convert vertices to homogeneous coordinates
    ones = np.ones((vertices.shape[0], 1))
    homogeneous_vertices = np.hstack([vertices, ones])

    # Apply affine transformation
    transformed_homogeneous = homogeneous_vertices @ affine_matrix.T

    # Convert back to Cartesian coordinates
    transformed_vertices = transformed_homogeneous[:, :3]

    return transformed_vertices


def remove_nan_vertices(vertices, faces):
    """
    Removes vertices containing NaNs and updates faces accordingly.

    Parameters:
    - vertices (np.ndarray): (N, 3) array of vertex positions.
    - faces (np.ndarray): (M, 3) array of triangular face indices.

    Returns:
    - new_vertices (np.ndarray): Filtered (N', 3) array of valid vertex positions.
    - new_faces (np.ndarray): Filtered (M', 3) array of updated face indices.
    """
    # Identify valid (non-NaN) vertices
    valid_mask = ~np.isnan(vertices).any(axis=1)

    # Create a mapping from old indices to new indices
    new_indices = np.full(
        vertices.shape[0], -1, dtype=int
    )  # Default -1 for invalid ones
    new_indices[valid_mask] = np.arange(valid_mask.sum())  # Renumber valid vertices

    # Filter out faces that reference removed vertices
    valid_faces_mask = np.all(
        valid_mask[faces], axis=1
    )  # Keep only faces with valid vertices
    new_faces = new_indices[faces[valid_faces_mask]]  # Remap face indices

    # Filter vertices
    new_vertices = vertices[valid_mask]

    return new_vertices, new_faces


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
coords[src_mask == 1] = -0.1
coords[sink_mask == 1] = 1.1


# Add the scalar field
grid.cell_data["values"] = coords.flatten(order="F")
grid = grid.cells_to_points("values")

# apply the affine
tfm_grid = grid.transform(
    affine, inplace=False
)  # Apply the rotation part of the affine


# the contour function produces the isosurface
surface = tfm_grid.contour([snakemake.params.threshold], method="contour").decimate(0.9)
# faces from pyvista surface are formatted with number of verts each row
# reshape and remove the first col to get Nx3
faces = surface.faces
faces = faces.reshape((int(faces.shape[0] / 4), 4))[:, 1:4]
points = surface.points
points, faces = remove_nan_vertices(points, faces)

## JD clean - instead of trimming surfaces with a nan mask, we
# keep vertices that overlap with good coord values. We then apply
# some surface-based morphological opening and closing to keep
# vertices along holes in the dg

# this is equivalent to wb_command -volume-to-surface-mapping -enclosing
# apply inverse affine to surface to get back to matrix space
V = apply_affine_transform(points, affine, inverse=True)
V = V.astype(int)
# sample coords
coord_at_V = np.zeros((len(V)))
for i in range(len(V)):
    coord_at_V[i] = coords[
        V[i, 0], V[i, 1], V[i, 2]
    ]  # really hope there's no x-y switching fuckery here!

# keep vertices that are in a nice coordinate range
epsilon = snakemake.params.coords_epsilon
good_v = np.where(np.logical_and(coord_at_V < (1 - epsilon), coord_at_V > epsilon))[0]

geoalg = geodesic.PyGeodesicAlgorithmExact(points, faces)
# morphological open
maxdist, _ = geoalg.geodesicDistances(good_v, None)
bad_v = np.where(maxdist > snakemake.params.morph_openclose_dist)[0]
# morphological close
maxdist, _ = geoalg.geodesicDistances(bad_v, None)
bad_v = np.where(maxdist < snakemake.params.morph_openclose_dist)[0]

# toss bad vertices
points[bad_v, :] = np.nan
points, faces = remove_nan_vertices(points, faces)


# apply largest connected component
faces_pv = np.hstack([np.full((faces.shape[0], 1), 3), faces])
mesh = pv.PolyData(points, faces_pv)
mesh_cc = mesh.extract_largest()
points = mesh_cc.points  # This gives you the vertices
faces = mesh_cc.faces
faces = faces.reshape((int(faces.shape[0] / 4), 4))[:, 1:4]

# write to gifti
write_surface_to_gifti(points, faces, snakemake.output.surf_gii)
