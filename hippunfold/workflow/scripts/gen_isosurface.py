import pyvista as pv
import nibabel as nib
import numpy as np
from copy import deepcopy


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


import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra


def compute_geodesic_distances(vertices, faces, source_indices):
    """
    Computes geodesic distances from a set of source vertices to all other vertices on a 3D surface mesh.

    Parameters:
    - vertices (np.ndarray): (N, 3) array of vertex positions.
    - faces (np.ndarray): (M, 3) array of triangular face indices.
    - source_indices (list or np.ndarray): Indices of source vertices.

    Returns:
    - distances (np.ndarray): (N,) array of geodesic distances from the source vertices.
    """
    num_vertices = len(vertices)

    # Create adjacency matrix
    row, col, weight = [], [], []
    for f in faces:
        for i in range(3):
            v1, v2 = f[i], f[(i + 1) % 3]  # Pairwise edges in the triangle
            dist = np.linalg.norm(vertices[v1] - vertices[v2])  # Euclidean edge length
            row.append(v1)
            col.append(v2)
            weight.append(dist)
            row.append(v2)
            col.append(v1)
            weight.append(dist)  # Ensure symmetry

    graph = sp.csr_matrix((weight, (row, col)), shape=(num_vertices, num_vertices))

    # Compute geodesic distances using Dijkstra's algorithm
    distances = dijkstra(csgraph=graph, directed=False, indices=source_indices)

    # If multiple sources, take the minimum distance to any of them
    if isinstance(source_indices, (list, np.ndarray)) and len(source_indices) > 1:
        distances = np.min(distances, axis=0)

    return distances


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
V = deepcopy(points)
V[:, :] = V - affine[:3, 3].T
for xyz in range(3):
    V[:, xyz] = V[:, xyz] * (1 / affine[xyz, xyz])
V = V.astype(int)
# sample coords
coord_at_V = np.zeros((len(V)))
for i in range(len(V)):
    coord_at_V[i] = coords[
        V[i, 0], V[i, 1], V[i, 2]
    ]  # really hope there's no x-y switching fuckery here!

# keep vertices that are in a nice coordinate range
epsilon - snakemake.params.coords_epsilon
good_v = np.where(np.logical_and(coord_at_V < (1 - epsilon), coord_at_V > epsilon))[0]

# morphological open
maxdist = compute_geodesic_distances(points, faces, good_v)
bad_v = np.where(maxdist > snakemake.params.morph_openclose_dist)[0]

# morphological close
maxdist = compute_geodesic_distances(points, faces, bad_v)
bad_v = np.where(maxdist < snakemake.params.morph_openclose_dist)[0]

# toss bad vertices
points[bad_v, :] = np.nan
points, faces = remove_nan_vertices(points, faces)


# write to gifti
write_surface_to_gifti(points, faces, snakemake.output.surf_gii)
