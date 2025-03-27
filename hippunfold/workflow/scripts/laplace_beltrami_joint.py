import numpy as np
import nibabel as nib
from scipy.sparse import diags, linalg, lil_matrix
from scipy.spatial import cKDTree
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def cotangent_laplacian(vertices, faces):
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
    diagonal = weights.sum(axis=1).A1
    # Ensure no zero entries in diagonal to avoid singular matrix issues
    diagonal[diagonal < 1e-12] = 1e-12
    laplacian = diags(diagonal) - weights
    return laplacian


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

surf1 = nib.load(snakemake.input.surf_gii_hipp)
vertices1 = surf.agg_data("NIFTI_INTENT_POINTSET")
faces1 = surf.agg_data("NIFTI_INTENT_TRIANGLE")

surf2 = nib.load(snakemake.input.surf_gii_dentate)
vertices2 = surf.agg_data("NIFTI_INTENT_POINTSET")
faces2 = surf.agg_data("NIFTI_INTENT_TRIANGLE")

src_sink_mask1 = nib.load(snakemake.input.src_sink_mask_hipp).agg_data()
src_indices1 = np.where(src_sink_mask == 1)[0]
sink_indices1 = np.where(src_sink_mask == 2)[0]

src_sink_mask1 = nib.load(snakemake.input.src_sink_mask_dentate).agg_data()
src_indices1 = np.where(src_sink_mask == 1)[0]
sink_indices1 = np.where(src_sink_mask == 2)[0]

## stack
vertices = np.vstack(vertices1, vertices2)
faces = np.vstack(faces1, faces2+len(vertices1))

src_indices = np.vstack(src_indices1, src_indices2+len(vertices1))
sink_indices = np.vstack(sink_indices1, sink_indices2+len(vertices1))

## bridgeheads
bridge_hipp = nib.load(snakemake.input.bridgehead_hipp).agg_data()
bridge_hipp_indices = np.where(bridge_hipp == 2)[0]
bridge_dentate = nib.load(snakemake.input.bridgehead_dentate).agg_data()
bridge_dentate_indices = np.where(bridge_dentate == 1)[0] + len(vertices1)

## stitch
tree_dentate = cKDTree(vertices[bridge_dentate_indices])
for e in bridge_hipp_indices:
    nearestvert = tree_dentate.query(vertices[e,:])
    faces = faces.vstack(faces, [e, nearestvert, e])
tree_hipp = cKDTree(vertices[bridge_dentate_indices])
for e in bridge_dentate_indices:
    nearestvert = tree_hipp.query(vertices[e,:])
    faces = faces.vstack(faces, [e, nearestvert, e])


logger.info(f"# of src boundary vertices: {len(src_indices)}")
logger.info(f"# of sink boundary vertices: {len(sink_indices)}")


src_vals = [0 for i in range(len(src_indices))]
sink_vals = [1 for i in range(len(sink_indices))]

boundary_conditions = dict(
    zip(list(src_indices) + list(sink_indices), src_vals + sink_vals)
)


coords = solve_laplace_beltrami_open_mesh(vertices, faces, boundary_conditions)


## unstitch
coords_hipp = coords[:len(vertices)]
coords_dentate = coords[len(vertices)+1:]

## save

# get structure metadata from src/sink mask
structure_metadata = nib.load(snakemake.input.src_sink_mask_hipp).meta[
    "AnatomicalStructurePrimary"
]

data_array = nib.gifti.GiftiDataArray(data=coords_hipp.astype(np.float32))
image = nib.gifti.GiftiImage()

# set structure metadata
image.meta["AnatomicalStructurePrimary"] = structure_metadata


image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords_hipp)

## save

structure_metadata = nib.load(snakemake.input.src_sink_mask_dentate).meta[
    "AnatomicalStructurePrimary"
]

data_array = nib.gifti.GiftiDataArray(data=coords_dentate.astype(np.float32))
image = nib.gifti.GiftiImage()

# set structure metadata
image.meta["AnatomicalStructurePrimary"] = structure_metadata


image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords_dentate)
