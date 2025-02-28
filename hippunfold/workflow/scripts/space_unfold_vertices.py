import numpy as np
import pyvista as pv
import nibabel as nib
from lib.utils import setup_logger
from scipy.spatial import KDTree

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def read_surface_from_gifti(surf_gii):
    """Load a surface mesh from a GIFTI file."""
    surf = nib.load(surf_gii)
    vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
    faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")
    faces_pv = np.hstack([np.full((faces.shape[0], 1), 3), faces])  # PyVista format
    return pv.PolyData(vertices, faces_pv)


def write_surface_to_gifti(mesh, out_surf_gii):
    """Writes a PyVista mesh to a GIFTI surface file."""
    points = mesh.points
    faces = mesh.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices

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


# Run the processing pipeline
mesh = read_surface_from_gifti(snakemake.input.surf_gii)
mesh_native = read_surface_from_gifti(snakemake.input.native_gii)

# check for bad vertices
# orig_points = mesh.points
# mesh = mesh.clean()
# if not mesh.points.shape == mesh_native.points.shape:
#     logger.info("unfolded surface contained degenerate vertices. removing them for now, but the native surf will also need to be modified accordingly!")
#     tree = KDTree(orig_points)
#     _, retained_indices = tree.query(mesh.points, k=1)  # Find closest match indices
#     mesh_native = mesh_native.extract_points(retained_indices)

"""
Perform spring-based smoothing on a 2D projected mesh, aiming for uniform edge lengths.
Parameters:
    mesh (pv.PolyData): Input 2D mesh.
    iterations (int): Number of smoothing iterations.
    step_size (float): Controls how much vertices move per step.
Returns:
    pv.PolyData: Smoothed mesh.
"""
points = mesh.points
# points += (np.random.rand(points.shape[0],points.shape[1])-.5) *1e-6 # add some noise to avoid perfectly overlapping points
faces = mesh.faces.reshape(-1, 4)[:, 1:4]  # Extract edges as pairs of vertex indices
edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
edges = np.sort(edges, axis=1)  # Ensure ordering for uniqueness
edges = np.unique(edges, axis=0)
edge_vectors = points[edges[:, 1]] - points[edges[:, 0]]
logger.info(edges.shape)

edge_vectors_native = mesh_native.points[edges[:, 1]] - mesh_native.points[edges[:, 0]]
target_length = np.linalg.norm(edge_vectors_native, axis=1)
logger.info(f"mean target edge length: {np.mean(target_length)}")

for i in range(snakemake.params.max_iterations):
    # Compute current edge lengths and stretching forces
    edge_vectors = points[edges[:, 1]] - points[edges[:, 0]]
    current_lengths = np.linalg.norm(edge_vectors, axis=1)
    mask = current_lengths > 1e-6  # Avoid division by zero

    forces = np.zeros_like(edge_vectors)
    forces[mask] = (
        (target_length[mask] - current_lengths[mask])[:, None]
        * edge_vectors[mask]
        / current_lengths[mask][:, None]
    )

    # Accumulate forces at each vertex
    displacement = np.zeros_like(points)
    np.add.at(displacement, edges[:, 0], -snakemake.params.step_size * forces)
    np.add.at(displacement, edges[:, 1], snakemake.params.step_size * forces)

    # Update vertex positions
    points += displacement

    if i % 100 == 0:
        mae = np.mean(np.abs(displacement))
        logger.info(f"iteration {i} mean absolute error: {mae}")
        if mae < 1e-6:
            logger.info("stopping early")
            break

write_surface_to_gifti(mesh, snakemake.output.surf_gii)
