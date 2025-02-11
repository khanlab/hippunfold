import numpy as np
import nibabel as nib
from scipy.sparse import diags, linalg, lil_matrix
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

import numpy as np

import numpy as np
from scipy.signal import argrelextrema


def get_terminal_indices_firstminima(
    sdt, min_vertices, boundary_mask, bins=100, smoothing_window=5
):
    """
    Gets the terminal (src/sink) vertex indices by determining an adaptive threshold
    using the first local minimum of the histogram of `sdt` values.

    Parameters:
    - sdt: Signed distance transform array.
    - min_vertices: The minimum number of vertices required.
    - boundary_mask: Boolean or binary mask indicating boundary regions.
    - bins: Number of bins to use in the histogram (default: 100).
    - smoothing_window: Window size for moving average smoothing (default: 5).

    Returns:
    - indices: List of terminal vertex indices.

    Raises:
    - ValueError: If the minimum number of vertices cannot be found.
    """

    # Extract SDT values within the boundary mask
    sdt_values = sdt[boundary_mask == 1]

    # Compute histogram
    hist, bin_edges = np.histogram(sdt_values, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Smooth the histogram using a simple moving average
    smoothed_hist = np.convolve(
        hist, np.ones(smoothing_window) / smoothing_window, mode="same"
    )

    # Find local minima
    minima_indices = argrelextrema(smoothed_hist, np.less)[0]

    if len(minima_indices) == 0:
        raise ValueError("No local minima found in the histogram.")

    # Select the first local minimum after the first peak
    first_minimum_bin = bin_centers[minima_indices[0]]

    # Select indices where SDT is below this threshold
    indices = np.where((sdt < first_minimum_bin) & (boundary_mask == 1))[0].tolist()

    if len(indices) >= min_vertices:
        return indices

    raise ValueError(
        f"Unable to find minimum of {min_vertices} vertices on boundary within the first local minimum threshold."
    )


def get_terminal_indices_percentile(
    sdt, min_percentile, max_percentile, min_vertices, boundary_mask
):
    """
    Gets the terminal (src/sink) vertex indices by sweeping a percentile-based threshold
    of the signed distance transform (sdt), ensuring at least `min_vertices` are selected.

    Instead of a fixed distance range, this function dynamically determines the threshold
    by scanning from `min_percentile` to `max_percentile`.

    Parameters:
    - sdt: Signed distance transform array.
    - min_percentile: Starting percentile for thresholding (0-100).
    - max_percentile: Maximum percentile for thresholding (0-100).
    - min_vertices: The minimum number of vertices required.
    - boundary_mask: Boolean or binary mask indicating boundary regions.

    Returns:
    - indices: List of terminal vertex indices.

    Raises:
    - ValueError: If the minimum number of vertices cannot be found.
    """

    for percentile in np.arange(min_percentile, max_percentile, 0.5):
        dist_threshold = np.percentile(sdt[boundary_mask == 1], percentile)
        indices = np.where((sdt < dist_threshold) & (boundary_mask == 1))[0].tolist()

        if len(indices) >= min_vertices:
            logger.info(
                f"Using {percentile}-th percentile to obtain sdt threshold of {dist_threshold}, with {len(indices)} vertices"
            )
            return indices

    raise ValueError(
        f"Unable to find minimum of {min_vertices} vertices on boundary within the {max_percentile}th percentile of distances"
    )


def get_terminal_indices_threshold(
    sdt, min_dist, max_dist, min_vertices, boundary_mask
):
    """
    Gets the terminal (src/sink) vertex indices based on distance to the src/sink mask,
    a boundary mask, and a minumum number of vertices. The distance from the mask is
    swept from min_dist to max_dist, until the min_vertices is achieved, else an
    exception is thrown."""

    for dist in np.linspace(min_dist, max_dist, 20):
        indices = np.where((sdt < dist) & (boundary_mask == 1))[0].tolist()
        if len(indices) >= min_vertices:
            return indices
    raise ValueError(
        f"Unable to find minimum of {min_vertices} vertices on boundary, within {max_dist}mm of the terminal mask"
    )


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
        # Avoid division by zero
        cot0 = np.dot(e1, -e2) / norm0 if norm0 > 1e-12 else 0.0
        cot1 = np.dot(e2, -e0) / norm1 if norm1 > 1e-12 else 0.0
        cot2 = np.dot(e0, -e1) / norm2 if norm2 > 1e-12 else 0.0
        weights[tri[0], tri[1]] += cot2 / 2
        weights[tri[1], tri[0]] += cot2 / 2
        weights[tri[1], tri[2]] += cot0 / 2
        weights[tri[2], tri[1]] += cot0 / 2
        weights[tri[2], tri[0]] += cot1 / 2
        weights[tri[0], tri[2]] += cot1 / 2
    logger.info("weights.tocsr()")
    weights = weights.tocsr()
    # Step 2: Handle boundaries for open meshes
    logger.info("Handle boundaries for open meshes")
    diagonal = weights.sum(axis=1).A1
    # Ensure no zero entries in diagonal to avoid singular matrix issues
    diagonal[diagonal < 1e-12] = 1e-12
    laplacian = diags(diagonal) - weights
    if boundary_conditions is None:
        boundary_conditions = {}
    boundary_indices = list(boundary_conditions.keys())
    boundary_values = np.array(list(boundary_conditions.values()))
    free_indices = np.setdiff1d(np.arange(n_vertices), boundary_indices)
    logger.info("Setting boundary conditions")
    b = np.zeros(n_vertices)
    for idx, value in boundary_conditions.items():
        laplacian[idx, :] = 0
        laplacian[idx, idx] = 1
        b[idx] = value
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

boundary_mask = nib.load(snakemake.input.boundary).agg_data()
src_sdt = nib.load(snakemake.input.src_sdt).agg_data()
sink_sdt = nib.load(snakemake.input.sink_sdt).agg_data()

if snakemake.params.threshold_method == "percentile":
    src_indices = get_terminal_indices_percentile(
        src_sdt,
        snakemake.params.min_dist_percentile,
        snakemake.params.max_dist_percentile,
        snakemake.params.min_terminal_vertices,
        boundary_mask,
    )
    sink_indices = get_terminal_indices_percentile(
        sink_sdt,
        snakemake.params.min_dist_percentile,
        snakemake.params.max_dist_percentile,
        snakemake.params.min_terminal_vertices,
        boundary_mask,
    )


elif snakemake.params.threshold_method == "firstminima":
    src_indices = get_terminal_indices_firstminima(
        src_sdt,
        snakemake.params.min_terminal_vertices,
        boundary_mask,
    )
    sink_indices = get_terminal_indices_firstminima(
        sink_sdt,
        snakemake.params.min_terminal_vertices,
        boundary_mask,
    )


logger.info(f"# of src boundary vertices: {len(src_indices)}")
logger.info(f"# of sink boundary vertices: {len(sink_indices)}")

src_vals = [0 for i in range(len(src_indices))]
sink_vals = [1 for i in range(len(sink_indices))]

boundary_conditions = dict(zip(src_indices + sink_indices, src_vals + sink_vals))


coords = solve_laplace_beltrami_open_mesh(vertices, faces, boundary_conditions)

data_array = nib.gifti.GiftiDataArray(data=coords.astype(np.float32))
image = nib.gifti.GiftiImage()
image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords)
