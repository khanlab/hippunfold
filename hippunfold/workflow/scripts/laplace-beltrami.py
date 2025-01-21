import numpy as np
from scipy.sparse import diags, linalg, lil_matrix
import nibabel as nib
from scipy.interpolate import NearestNDInterpolator
from collections import defaultdict
import networkx as nx


def largest_connected_component(vertices: np.ndarray, faces: np.ndarray):
    """
    Returns the vertices, faces, and original indices of the largest connected component of a 3D surface mesh.

    Parameters:
    vertices (np.ndarray): Array of vertex coordinates (N x 3).
    faces (np.ndarray): Array of face indices (M x 3).

    Returns:
    tuple: (new_vertices, new_faces, largest_component), where new_vertices are the vertex coordinates of the largest component,
           new_faces are the face indices adjusted to the new vertex order, and largest_component contains the indices of the original vertices.
    """
    # Build adjacency graph from face connectivity
    G = nx.Graph()
    for face in faces:
        G.add_edges_from(
            [(face[i], face[j]) for i in range(3) for j in range(i + 1, 3)]
        )

    # Find connected components
    components = list(nx.connected_components(G))

    # Select the largest connected component
    largest_component = max(components, key=len)
    largest_component = np.array(list(largest_component))

    # Create a mapping from old vertex indices to new ones
    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(largest_component)}

    # Filter the vertices and faces
    new_vertices = vertices[largest_component]
    new_faces = np.array(
        [
            [index_map[v] for v in face]
            for face in faces
            if all(v in index_map for v in face)
        ]
    )

    return new_vertices, new_faces, largest_component


def find_boundary_vertices(vertices, faces):
    """
    Find boundary vertices of a 3D mesh.

    Args:
        vertices (list of tuples): List of 3D points (x, y, z).
        faces (list of tuples): List of triangular faces, where each face
                                is a tuple of three vertex indices (v1, v2, v3).

    Returns:
        list: List of vertex indices that are boundary vertices, sorted in ascending order.
    """
    edge_count = defaultdict(int)
    # Step 1: Count edge occurrences
    for face in faces:
        # Extract edges from the face, ensure consistent ordering (min, max)
        edges = [
            tuple(sorted((face[0], face[1]))),
            tuple(sorted((face[1], face[2]))),
            tuple(sorted((face[2], face[0]))),
        ]
        for edge in edges:
            edge_count[edge] += 1
    # Step 2: Identify boundary edges
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
    # Step 3: Collect boundary vertices
    boundary_vertices = set()
    for edge in boundary_edges:
        boundary_vertices.update(edge)
    # Convert the set to a sorted list (array)
    return sorted(boundary_vertices)


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
    # Step 1: Compute cotangent weights
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
    weights = weights.tocsr()
    # Step 2: Handle boundaries for open meshes
    diagonal = weights.sum(axis=1).A1
    # Ensure no zero entries in diagonal to avoid singular matrix issues
    diagonal[diagonal < 1e-12] = 1e-12
    laplacian = diags(diagonal) - weights
    if boundary_conditions is None:
        boundary_conditions = {}
    boundary_indices = list(boundary_conditions.keys())
    boundary_values = np.array(list(boundary_conditions.values()))
    free_indices = np.setdiff1d(np.arange(n_vertices), boundary_indices)
    b = np.zeros(n_vertices)
    for idx, value in boundary_conditions.items():
        laplacian[idx, :] = 0
        laplacian[idx, idx] = 1
        b[idx] = value
    # Step 3: Solve the Laplace-Beltrami equation
    solution = np.zeros(n_vertices)
    if len(free_indices) > 0:
        free_laplacian = laplacian[free_indices][:, free_indices]
        free_b = (
            b[free_indices]
            - laplacian[free_indices][:, boundary_indices] @ boundary_values
        )
        solution[boundary_indices] = boundary_values
        try:
            solution[free_indices] = linalg.spsolve(free_laplacian, free_b)
        except linalg.MatrixRankWarning:
            print("Warning: Laplacian matrix is singular or ill-conditioned.")
            solution[free_indices] = np.zeros(len(free_indices))
    else:
        solution[boundary_indices] = boundary_values
    return solution


surf = nib.load(snakemake.input.surf_gii)
vertices_orig = surf.agg_data("NIFTI_INTENT_POINTSET")
faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")

vertices, faces, i_concomp = largest_connected_component(vertices_orig, faces)


# get source/sink vertices by nearest neighbour (only of edge vertices)
boundary_vertices = np.array(find_boundary_vertices(vertices, faces))
seg = nib.load(snakemake.input.seg)

src_AP = np.array(np.where(seg.get_fdata() == snakemake.params.src_labels["AP"]["src"]))
sink_AP = np.array(
    np.where(seg.get_fdata() == snakemake.params.src_labels["AP"]["sink"])
)
src_PD = np.array(np.where(seg.get_fdata() == snakemake.params.src_labels["PD"]["src"]))
sink_PD = np.array(
    np.where(seg.get_fdata() == snakemake.params.src_labels["PD"]["sink"])
)

# apply affine
src_AP = (seg.affine @ np.vstack([src_AP, np.ones([1, src_AP.shape[1]])]))[:3, :]
sink_AP = (seg.affine @ np.vstack([sink_AP, np.ones([1, sink_AP.shape[1]])]))[:3, :]
src_PD = (seg.affine @ np.vstack([src_PD, np.ones([1, src_PD.shape[1]])]))[:3, :]
sink_PD = (seg.affine @ np.vstack([sink_PD, np.ones([1, sink_PD.shape[1]])]))[:3, :]

vals = np.hstack(
    [
        np.ones([src_AP.shape[1]]) * 10,
        np.ones([sink_AP.shape[1]]) * 11,
        np.ones([src_PD.shape[1]]) * 20,
        np.ones([sink_PD.shape[1]]) * 21,
    ]
)
interpol = NearestNDInterpolator(np.hstack([src_AP, sink_AP, src_PD, sink_PD]).T, vals)
boundary_values = np.array(interpol(vertices[boundary_vertices, :])).astype(int)

APinds = np.array(np.where(boundary_values < 12)[0]).astype(int)
boundary_conditions = dict(zip(boundary_vertices[APinds], boundary_values[APinds] - 10))
APcoords = np.nans((len(vertices_orig)))
APcoords[i_concomp] = solve_laplace_beltrami_open_mesh(
    vertices, faces, boundary_conditions
)
PDinds = np.array(np.where(boundary_values > 12)[0]).astype(int)
boundary_conditions = dict(zip(boundary_vertices[PDinds], boundary_values[PDinds] - 20))
PDcoords = np.nans((len(vertices_orig)))
PDcoords[i_concomp] = solve_laplace_beltrami_open_mesh(
    vertices, faces, boundary_conditions
)

data_array = nib.gifti.GiftiDataArray(data=APcoords.astype(np.float32))
image = nib.gifti.GiftiImage()
image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords_AP)

data_array = nib.gifti.GiftiDataArray(data=PDcoords.astype(np.float32))
image = nib.gifti.GiftiImage()
image.add_gifti_data_array(data_array)
nib.save(image, snakemake.output.coords_PD)
