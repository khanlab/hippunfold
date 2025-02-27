import numpy as np
import nibabel as nib
from scipy.sparse import diags, linalg, lil_matrix, csr_matrix
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

def cotangent_laplacian(vertices, faces):
    """Computes the cotangent-weighted Laplacian matrix using a fully vectorized approach."""
    n = len(vertices)

    # Extract vertex indices for each triangle
    i1, i2, i3 = faces[:, 0], faces[:, 1], faces[:, 2]

    # Get corresponding vertex positions
    v1, v2, v3 = vertices[i1], vertices[i2], vertices[i3]

    # Compute edge vectors
    e1 = v2 - v1
    e2 = v3 - v1
    e3 = v3 - v2

    # Compute cotangents of angles using cross and dot product
    cross_e1e2 = np.cross(e1, e2)  # Normal of triangle
    norm_cross = np.linalg.norm(cross_e1e2, axis=1)  # Triangle area factor
    norm_cross[norm_cross<1e-6] = 1e-6

    # Cotangent weights
    cot_alpha = np.einsum('ij,ij->i', e1, e2) / norm_cross
    cot_beta = np.einsum('ij,ij->i', -e2, e3) / norm_cross
    cot_gamma = np.einsum('ij,ij->i', -e3, e1) / norm_cross

    # Construct sparse Laplacian matrix
    I = np.hstack([i1, i2, i2, i3, i3, i1])
    J = np.hstack([i2, i1, i3, i2, i1, i3])
    W = np.hstack([cot_alpha, cot_alpha, cot_beta, cot_beta, cot_gamma, cot_gamma])

    L = csr_matrix((W, (I, J)), shape=(n, n))
    D = L.sum(axis=1).A1
    L = L - diags(D)  # Ensure row sum is zero

    return L


def solve_laplace_beltrami_open_mesh(vertices, faces, boundary_conditions=None):
    """
    Solve the Laplace-Beltrami equation on a 3D open-faced surface mesh. No islands please!

    Parameters:
        vertices (np.ndarray): Array of shape (n_vertices, 3) containing vertex coordinates.
        faces (np.ndarray): Array of shape (n_faces, 3) containing indices of vertices forming each triangular face.
        boundary_conditions (dict, optional): Dictionary where keys are vertex indices with fixed values.

    Returns:
        solution (np.ndarray): Array of shape (n_vertices,) with the solution values.
    """
    logger.info("solve_laplace_beltrami_open_mesh")
    n_vertices = vertices.shape[0]
    logger.info(f"n_vertices: {n_vertices}")
    # Step 1: Compute cotangent weights
    logger.info("Computing cotangent weights")
    laplacian = cotangent_laplacian(vertices, faces)
    # Step 2: Handle boundaries for open meshes
    logger.info("Handle boundaries for open meshes")
    if boundary_conditions is None:
        boundary_conditions = {}
    boundary_indices = np.array(list(boundary_conditions.keys()))
    boundary_values = np.array(list(boundary_conditions.values()))
    free_indices = np.setdiff1d(np.arange(n_vertices), boundary_indices)
    logger.info("Setting boundary conditions")
    laplacian[boundary_indices, :] = 0  # Zero out entire rows
    laplacian[boundary_indices, boundary_indices] = 1  # Set diagonal entries to 1
    b = np.zeros(n_vertices)
    b[boundary_indices] = boundary_values
    # Step 3: Solve the Laplace-Beltrami equation
    logger.info("Solve the Laplace-Beltrami equation")
    solution = np.zeros(n_vertices)
    if len(free_indices) > 0:
        free_laplacian = laplacian[free_indices][:, free_indices]
        free_b = (
            b[free_indices]
            - laplacian[free_indices][:, boundary_indices] @ boundary_values
        )
        solution[boundary_indices] = boundary_values
        try:
            logger.info("about to solve")
            solution[free_indices] = linalg.spsolve(free_laplacian, free_b)
            logger.info("done solve")
        except linalg.MatrixRankWarning:
            logger.info("Warning: Laplacian matrix is singular or ill-conditioned.")
            solution[free_indices] = np.zeros(len(free_indices))
    else:
        solution[boundary_indices] = boundary_values
    return solution


logger.info(
    "Loading in surface, boundary mask, and src/sink signed distance transforms"
)

surf = nib.load(snakemake.input.surf_gii)
vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")

src_sink_mask = nib.load(snakemake.input.src_sink_mask).agg_data()
src_indices = np.where(src_sink_mask == 1)[0]
sink_indices = np.where(src_sink_mask == 2)[0]


logger.info(f"# of src boundary vertices: {len(src_indices)}")
logger.info(f"# of sink boundary vertices: {len(sink_indices)}")


src_vals = [0 for i in range(len(src_indices))]
sink_vals = [1 for i in range(len(sink_indices))]

boundary_conditions = dict(
    zip(list(src_indices) + list(sink_indices), src_vals + sink_vals)
)


coords = solve_laplace_beltrami_open_mesh(vertices, faces, boundary_conditions)

data_array = nib.gifti.GiftiDataArray(data=coords.astype(np.float32))
image = nib.gifti.GiftiImage()
image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords)
