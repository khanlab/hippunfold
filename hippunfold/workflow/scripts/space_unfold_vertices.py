import numpy as np
import pyvista as pv
import nibabel as nib
from lib.utils import setup_logger

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

def spring_based_smoothing(mesh, iterations=10, step_size=0.1):
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
    #points += (np.random.rand(points.shape[0],points.shape[1])-.5) *1e-6 # add some noise to avoid perfectly overlapping points
    faces = mesh.faces.reshape(-1, 4)[:, 1:4]  # Extract edges as pairs of vertex indices
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
    edges = np.sort(edges, axis=1)  # Ensure ordering for uniqueness
    edges = np.unique(edges, axis=0)
    edge_vectors = points[edges[:, 1]] - points[edges[:, 0]]
    logger.info(edges.shape)

    current_lengths = np.linalg.norm(edge_vectors, axis=1)
    target_length = np.zeros_like(current_lengths)
    target_length[:] = np.mean(current_lengths)
    logger.info(f'target edge length: {np.mean(current_lengths)}')

    for _ in range(iterations):
        # Compute current edge lengths and stretching forces
        edge_vectors = points[edges[:, 1]] - points[edges[:, 0]]
        current_lengths = np.linalg.norm(edge_vectors, axis=1)
        mask = current_lengths > 1e-6  # Avoid division by zero

        forces = np.zeros_like(edge_vectors)
        forces[mask] = (
            (target_length[mask] - current_lengths[mask])[:,None]
            * edge_vectors[mask] / current_lengths[mask][:,None]
        )
        logger.info(f"max force size {np.max(forces,axis=0)}")
        logger.info(f"min force size {np.min(forces,axis=0)}")

        # Accumulate forces at each vertex
        displacement = np.zeros_like(points)
        np.add.at(displacement, edges[:, 0], -step_size * forces)
        np.add.at(displacement, edges[:, 1], step_size * forces)

        # Update vertex positions
        points += displacement

        if np.sum(np.abs(displacement)) < 1e-6:
            logger.info('stopping early')
            break

    return mesh


# Run the processing pipeline
surface = read_surface_from_gifti(snakemake.input.surf_gii)
surface_smoothed = spring_based_smoothing(
    surface, iterations=snakemake.params.max_iterations, step_size=snakemake.params.step_size
)
write_surface_to_gifti(surface_smoothed, snakemake.output.surf_gii)
