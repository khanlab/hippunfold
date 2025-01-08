import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree


def replace_outlier_vertices(vertices, faces, max_edge_length=0.3):
    """
    Regularize a mesh by replacing outlier vertices with the mean of neighbors.

    Parameters:
    - vertices (numpy.ndarray): Vertex coordinates of shape (N, 3).
    - faces (numpy.ndarray): Face indices of shape (M, 3).
    - max_edge_length (float): Threshold for the maximum allowable edge length.

    Returns:
    - numpy.ndarray: Regularized vertex coordinates.
    """

    def compute_edges(vertices, faces):
        """Compute all edges and their lengths in the mesh."""
        edges = set()
        for face in faces:
            for i, j in zip(face, np.roll(face, 1)):
                edges.add(tuple(sorted((i, j))))
        edges = np.array(list(edges))
        lengths = np.linalg.norm(vertices[edges[:, 0]] - vertices[edges[:, 1]], axis=1)
        return edges, lengths

    def move_vertex(vertices, edge, midpoint):
        """Move the vertices of an edge to their midpoint."""
        vertices[edge[0]] = midpoint
        vertices[edge[1]] = midpoint

    vertices = vertices.copy()

    while True:
        edges, lengths = compute_edges(vertices, faces)

        # Find edges longer than the threshold
        long_edges = edges[lengths > max_edge_length]
        if len(long_edges) == 0:
            break

        # Move vertices of the longest edge to its midpoint
        for edge in long_edges:
            midpoint = vertices[edge].mean(axis=0)
            move_vertex(vertices, edge, midpoint)

    return vertices


# Snakemake input/output paths
input_gii_path = snakemake.input[0]
output_gii_path = snakemake.output[0]

# Load the GIFTI file
gii = nib.load(input_gii_path)

# Extract vertices and faces from the GIFTI file
vertices = gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data
faces = gii.get_arrays_from_intent("NIFTI_INTENT_TRIANGLE")[0].data

# Regularize the mesh
regularized_vertices = replace_outlier_vertices(vertices, faces)

# Update the vertices in the GIFTI object
gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data = regularized_vertices

# Save the modified GIFTI file
nib.save(gii, output_gii_path)
