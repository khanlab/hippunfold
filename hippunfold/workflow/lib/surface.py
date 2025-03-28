import pyvista as pv
import nibabel as nib
import numpy as np


def read_surface_from_gifti(surf_gii):
    """Load a surface mesh from a GIFTI file."""
    surf = nib.load(surf_gii)
    vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
    faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")
    faces_pv = np.hstack([np.full((faces.shape[0], 1), 3), faces])  # PyVista format

    # Find the first darray that represents vertices (NIFTI_INTENT_POINTSET)
    vertices_darray = next(
        (
            darray
            for darray in surf.darrays
            if darray.intent == intent_codes["NIFTI_INTENT_POINTSET"]
        ),
        None,
    )

    # Extract metadata as a dictionary (return empty dict if no metadata)
    metadata = dict(vertices_darray.meta) if vertices_darray else {}

    return pv.PolyData(vertices, faces_pv), metadata

def compute_edge_lengths(surface):

    # Extract edges
    edges = surface.extract_all_edges()

    # Extract individual edge segments
    edge_lengths = []
    lines = edges.lines.reshape(-1, 3)  # Each row: [2, point1, point2]

    for line in lines:
        _, p1, p2 = line  # First entry is always "2" (pairs of points)
        length = np.linalg.norm(edges.points[p1] - edges.points[p2])
        edge_lengths.append(length)

    edge_lengths = np.array(edge_lengths)

    return edge_lengths


def write_surface_to_gifti(mesh, out_surf_gii, metadata=None):
    """
    Writes a PyVista mesh to a GIFTI surface file.
    """
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


def read_metric_from_gii(metric_gifti):
    return nib.load(metric_gifti).darrays[0].data


def write_label_gii(label_scalars, out_label_gii, label_dict={}, metadata=None):

    # Create a GIFTI label data array
    gii_data = gifti.GiftiDataArray(label_scalars, intent="NIFTI_INTENT_LABEL")

    # Create a Label Table (LUT)
    label_table = gifti.GiftiLabelTable()

    for label_name, label_dict in label_table.items():

        lbl = gifti.GiftiLabel(**label_dict)
        lbl.label = label_name
        label_table.labels.append(lbl)

    # Assign label table to GIFTI image
    gii_img = gifti.GiftiImage(darrays=[gii_data], labeltable=label_table)

    # set structure metadata
    if metadata is not None:
        gii_img.meta["AnatomicalStructurePrimary"] = metadata[
            "AnatomicalStructurePrimary"
        ]

    # Save the label file
    gii_img.to_filename(out_label_gii)


def apply_affine_transform(mesh, affine_matrix, inverse=False):
    """
    Applies an affine transformation to a PyVista mesh.
    """
    if inverse:
        affine_matrix = np.linalg.inv(affine_matrix)

    transformed_points = (
        np.hstack([mesh.points, np.ones((mesh.n_points, 1))]) @ affine_matrix.T
    )[:, :3]
    transformed_mesh = pv.PolyData(transformed_points, mesh.faces)
    return transformed_mesh


def remove_nan_vertices(mesh):
    """
    Removes NaN vertices from a PyVista mesh and updates faces accordingly.
    """
    valid_mask = ~np.isnan(mesh.points).any(axis=1)
    new_indices = np.full(mesh.n_points, -1, dtype=int)
    new_indices[valid_mask] = np.arange(valid_mask.sum())

    faces = mesh.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices
    valid_faces_mask = np.all(valid_mask[faces], axis=1)
    new_faces = new_indices[faces[valid_faces_mask]]
    new_faces_pv = np.hstack([np.full((new_faces.shape[0], 1), 3), new_faces])

    cleaned_mesh = pv.PolyData(mesh.points[valid_mask], new_faces_pv)
    return cleaned_mesh


def find_boundary_vertices(mesh):
    """
    Find boundary vertices of a 3D mesh.

    Args:
        mesh
    Returns:
        list: List of vertex indices that are boundary vertices, sorted in ascending order.
    """
    vertices = mesh.points
    faces = mesh.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices

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
    return np.array(sorted(boundary_vertices), dtype=np.int32)
