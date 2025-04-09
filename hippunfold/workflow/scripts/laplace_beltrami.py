import numpy as np
import nibabel as nib
import scipy.sparse as sp
from lib.utils import setup_logger
import heapq
import sys
from scipy.sparse import coo_matrix, diags

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)
sys.stdout = open(log_file, "a")
sys.stderr = open(log_file, "a")


def cotangent_laplacian(vertices, faces):
    """
    Computes the cotangent Laplacian for a triangular mesh.

    Parameters:
        vertices (np.ndarray): (n_vertices, 3) array of vertex positions.
        faces (np.ndarray): (n_faces, 3) array of mesh faces.

    Returns:
        scipy.sparse.csr_matrix: The Laplacian matrix (normalized if specified).
    """
    n_vertices = vertices.shape[0]

    # Step 1: Compute cotangent weights
    logger.info("Computing cotangent weights")
    weights = lil_matrix((n_vertices, n_vertices))

    for tri in faces:
        v0, v1, v2 = vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]
        e0 = v1 - v2
        e1 = v2 - v0
        e2 = v0 - v1

        # Compute cross products and norms
        cross0 = np.cross(e1, -e2)
        cross1 = np.cross(e2, -e0)
        cross2 = np.cross(e0, -e1)

        norm0 = np.linalg.norm(cross0)
        norm1 = np.linalg.norm(cross1)
        norm2 = np.linalg.norm(cross2)

        # Compute cotangent weights safely
        cot0 = np.dot(e1, -e2) / norm0 if norm0 > 1e-12 else 0.0
        cot1 = np.dot(e2, -e0) / norm1 if norm1 > 1e-12 else 0.0
        cot2 = np.dot(e0, -e1) / norm2 if norm2 > 1e-12 else 0.0

        weights[tri[0], tri[1]] += cot2 / 2
        weights[tri[1], tri[0]] += cot2 / 2
        weights[tri[1], tri[2]] += cot0 / 2
        weights[tri[2], tri[1]] += cot0 / 2
        weights[tri[2], tri[0]] += cot1 / 2
        weights[tri[0], tri[2]] += cot1 / 2

    logger.info("Converting to CSR format")
    weights = weights.tocsr()
    diagonal = weights.sum(axis=1).A1

    # Ensure no zero entries in diagonal to avoid singular matrix issues
    diagonal[diagonal < 1e-12] = 1e-12
    laplacian = diags(diagonal) - weights

    return laplacian


def fast_marching_mesh(vertices, faces, source_indices, normalize=True):
    """
    Computes geodesic distances using an optimized Fast Marching Method (FMM) on a triangular mesh.

    Parameters:
        vertices (np.ndarray): (n_vertices, 3) array of vertex positions.
        faces (np.ndarray): (n_faces, 3) array of mesh faces.
        source_indices (list or np.ndarray): List of source vertex indices (geodesic distance = 0).

    Returns:
        np.ndarray: Geodesic distances at each vertex.
    """
    n_vertices = len(vertices)

    # Step 1: Initialize distance field
    phi = np.full(n_vertices, np.inf)
    phi[source_indices] = 0  # Source vertices are at 0 distance

    # Step 2: Construct adjacency list for fast neighbor lookup
    adjacency = {i: set() for i in range(n_vertices)}
    for face in faces:
        for i in range(3):
            adjacency[face[i]].add(face[(i + 1) % 3])
            adjacency[face[i]].add(face[(i + 2) % 3])

    # Step 3: Priority queue for fast marching
    heap = [(0.0, idx) for idx in source_indices]
    heapq.heapify(heap)

    # Step 4: Fast Marching Loop
    visited = set(source_indices)

    while heap:
        dist, v_idx = heapq.heappop(heap)

        for neighbor in adjacency[v_idx]:
            if neighbor in visited:
                continue  # Skip if already processed

            # Compute geodesic update
            new_dist = np.linalg.norm(vertices[v_idx] - vertices[neighbor]) + dist

            if new_dist < phi[neighbor]:  # Only update if we found a shorter path
                phi[neighbor] = new_dist
                heapq.heappush(heap, (new_dist, neighbor))
                visited.add(neighbor)

    # Step 5: Normalize distances to [0,1]
    if normalize:
        phi_min, phi_max = np.min(phi), np.max(phi)
        phi = (
            (phi - phi_min) / (phi_max - phi_min)
            if phi_max - phi_min > 1e-12
            else np.zeros_like(phi)
        )

    return phi


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
    n_vertices = vertices.shape[0]
    logger.info("solve_laplace_beltrami_open_mesh")
    laplacian = cotangent_laplacian(vertices, faces)

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

surf = nib.load(snakemake.input.surf_gii)
vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")

src_sink_mask = nib.load(snakemake.input.src_sink_mask).agg_data()
src_indices = np.where(src_sink_mask == 1)[0]
sink_indices = np.where(src_sink_mask == 2)[0]

# get structure metadata from src/sink mask
structure_metadata = nib.load(snakemake.input.src_sink_mask).meta[
    "AnatomicalStructurePrimary"
]


logger.info(f"# of src boundary vertices: {len(src_indices)}")
logger.info(f"# of sink boundary vertices: {len(sink_indices)}")

src_vals = [0 for i in range(len(src_indices))]
sink_vals = [1 for i in range(len(sink_indices))]

boundary_conditions = dict(
    zip(list(src_indices) + list(sink_indices), src_vals + sink_vals)
)


if snakemake.params.method == "fastmarching":

    logger.info("calculating forward phi")
    forward_phi = fast_marching_mesh(vertices, faces, src_indices)
    logger.info(f"min: {forward_phi.min()}, max: {forward_phi.max()}")
    logger.info("calculating backward phi")
    backward_phi = fast_marching_mesh(vertices, faces, sink_indices)
    logger.info(f"min: {backward_phi.min()}, max: {backward_phi.max()}")

    # Average both distance maps
    logger.info("averaging forward and backward marching")
    coords = 0.5 * (
        forward_phi + (1 - backward_phi)
    )  # 1-backward_phi to invert its direction
    logger.info(f"min: {coords.min()}, max: {coords.max()}")

    coords = forward_phi

elif snakemake.params.method == "laplace":

    coords = solve_laplace_beltrami_open_mesh(vertices, faces, boundary_conditions)


data_array = nib.gifti.GiftiDataArray(data=coords.astype(np.float32))
image = nib.gifti.GiftiImage()

# set structure metadata
image.meta["AnatomicalStructurePrimary"] = structure_metadata


image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords)
