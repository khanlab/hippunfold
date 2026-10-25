import nibabel as nib
import numpy as np
import scipy.sparse as sp
from lib.surface import (
    get_boundary_loop,
    monotone_free_arcs,
    read_surface_from_gifti,
    write_metric_gii,
)
from lib.utils import setup_logger
from scipy.sparse import coo_matrix, diags

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def cotangent_laplacian(vertices, faces, clamp_negative=True):
    """
    Cotangent Laplacian. Each edge weight is half the sum of the cotangents of
    the angles opposite it. Obtuse angles make weights negative, which can fold
    the solution over (non-monotonic coordinates); clamp_negative sets negative
    edge weights to zero.
    """
    n_vertices = vertices.shape[0]
    rows, cols, vals = [], [], []
    for k in range(3):
        # edge (i, j) opposite corner o
        i, j, o = faces[:, (k + 1) % 3], faces[:, (k + 2) % 3], faces[:, k]
        u, w = vertices[i] - vertices[o], vertices[j] - vertices[o]
        cross = np.linalg.norm(np.cross(u, w), axis=1)
        safe = cross > 1e-12
        cot = np.zeros(len(faces))
        cot[safe] = np.einsum("ij,ij->i", u[safe], w[safe]) / cross[safe]
        rows += [i, j]
        cols += [j, i]
        vals += [cot / 2, cot / 2]
    weights = coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_vertices, n_vertices),
    ).tocsr()  # duplicate entries (an edge's two faces) are summed
    if clamp_negative:
        n_negative = int(np.sum(weights.data < 0)) // 2
        logger.info(f"Clamping {n_negative} negative edge weights to zero")
        weights.data = np.maximum(weights.data, 0)
    diagonal = weights.sum(axis=1).A1
    # Ensure no zero entries in diagonal to avoid singular matrix issues
    diagonal[diagonal < 1e-12] = 1e-12
    laplacian = diags(diagonal) - weights
    return laplacian


def solve_laplace_beltrami_open_mesh(
    vertices, faces, boundary_conditions=None, clamp_negative_weights=True
):
    """
    Solve the Laplace-Beltrami equation on a 3D open-faced surface mesh. No islands please!

    Parameters:
        vertices (np.ndarray): Array of shape (n_vertices, 3) containing vertex coordinates.
        faces (np.ndarray): Array of shape (n_faces, 3) containing indices of vertices forming each triangular face.
        boundary_conditions (dict, optional): Dictionary where keys are vertex indices with fixed values.

    Returns:
        solution (np.ndarray): Array of shape (n_vertices,) with the solution values.
    """
    n_vertices = vertices.shape[0]
    logger.info("solve_laplace_beltrami_open_mesh")
    laplacian = cotangent_laplacian(vertices, faces, clamp_negative_weights)

    # Step 2: Handle boundaries for open meshes
    logger.info("Handle boundaries for open meshes")
    if boundary_conditions is None:
        boundary_conditions = {}
    boundary_indices = np.array(list(boundary_conditions.keys()))
    boundary_values = np.array(list(boundary_conditions.values()))
    free_indices = np.setdiff1d(np.arange(n_vertices), boundary_indices)
    logger.info("Setting boundary conditions")
    for i in boundary_indices:
        start, end = laplacian.indptr[i], laplacian.indptr[i + 1]
        laplacian.data[start:end] = 0  # Zero out row
        laplacian[i, i] = 1  # Set diagonal to 1
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

            # Apply Tikhonov regularization
            eps = 1e-8 * laplacian.diagonal().mean()
            free_laplacian += eps * sp.eye(free_laplacian.shape[0], format="csr")

            solution[free_indices] = sp.linalg.spsolve(free_laplacian, free_b)

            logger.info("done solve")
        except sp.linalg.MatrixRankWarning:
            logger.info("Warning: Laplacian matrix is singular or ill-conditioned.")
            solution[free_indices] = np.zeros(len(free_indices))
    else:
        solution[boundary_indices] = boundary_values
    return solution


logger.info(
    "Loading in surface, boundary mask, and src/sink signed distance transforms"
)


surf, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
vertices = surf.points
faces = surf.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices

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


coords = solve_laplace_beltrami_open_mesh(
    vertices,
    faces,
    boundary_conditions,
    clamp_negative_weights=snakemake.params.clamp_negative_weights,
)

if snakemake.params.monotone_boundary:
    # make the coordinate monotone along the free arcs of the boundary loop, then
    # re-solve with the whole loop fixed. Where it already is, this gives the same
    # solution; elsewhere it removes folds of the unfolded mapping at the boundary.
    src_mask = np.zeros(len(vertices), dtype=bool)
    src_mask[src_indices] = True
    sink_mask = np.zeros(len(vertices), dtype=bool)
    sink_mask[sink_indices] = True
    loop = get_boundary_loop(faces, np.unique(faces))
    loop_values, n_changed = monotone_free_arcs(loop, coords, src_mask, sink_mask)
    logger.info(
        f"{n_changed} of {len(loop)} boundary loop vertices made monotone; "
        "re-solving with the whole loop fixed"
    )
    coords = solve_laplace_beltrami_open_mesh(
        vertices,
        faces,
        {**dict(zip(loop.tolist(), loop_values)), **boundary_conditions},
        clamp_negative_weights=snakemake.params.clamp_negative_weights,
    )

# make fully uniform distribution
domain = np.setdiff1d(
    np.arange(len(coords)), np.concatenate([src_indices, sink_indices])
)  # ignore src/sink indices
sorted_indices = np.argsort(coords[domain])
uniform_values = np.linspace(0, 1, len(domain) + 2)[
    1:-1
]  # ignore first and last values
coords[domain[sorted_indices]] = uniform_values


write_metric_gii(coords, snakemake.output.coords, metadata=metadata)
