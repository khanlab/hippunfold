from collections import defaultdict

import nibabel as nib
import numpy as np
import pyvista as pv
from nibabel.gifti.gifti import intent_codes


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


def write_points_faces_to_gifti(points, faces, out_surf_gii, metadata=None):

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


def write_surface_to_gifti(mesh, out_surf_gii, metadata=None):
    """
    Writes a PyVista mesh to a GIFTI surface file.
    """
    points = mesh.points
    faces = mesh.faces.reshape((-1, 4))[:, 1:4]  # Extract triangle indices

    write_points_faces_to_gifti(points, faces, out_surf_gii, metadata)


def read_metric_from_gii(metric_gifti):
    return nib.load(metric_gifti).darrays[0].data


def write_metric_gii(scalars, out_metric_gii, metadata=None):

    # save the coordinates to a gifti file
    data_array = nib.gifti.GiftiDataArray(data=scalars.astype(np.float32))
    image = nib.gifti.GiftiImage()

    # set structure metadata
    image.meta["AnatomicalStructurePrimary"] = metadata["AnatomicalStructurePrimary"]

    image.add_gifti_data_array(data_array)
    nib.save(image, out_metric_gii)


def write_label_gii(label_scalars, out_label_gii, label_dict=None, metadata=None):

    # Create a GIFTI label data array
    if label_dict is None:
        label_dict = {}
    gii_data = nib.gifti.GiftiDataArray(label_scalars, intent="NIFTI_INTENT_LABEL")

    # Create a Label Table (LUT)
    label_table = nib.gifti.GiftiLabelTable()

    for label_name, label_kwargs in label_dict.items():
        lbl = nib.gifti.GiftiLabel(**label_kwargs)
        lbl.label = label_name
        label_table.labels.append(lbl)

    # Assign label table to GIFTI image
    gii_img = nib.gifti.GiftiImage(darrays=[gii_data], labeltable=label_table)

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


def vertex_normals(points, faces):
    """area-weighted vertex normals from the face winding (not normalised)"""
    face_normals = np.cross(
        points[faces[:, 1]] - points[faces[:, 0]],
        points[faces[:, 2]] - points[faces[:, 0]],
    )
    normals = np.zeros(points.shape)
    for k in range(3):
        np.add.at(normals, faces[:, k], face_normals)
    return normals


def orient_faces_up_gradient(points, faces, volume, affine, step=None):
    """
    Reverse the winding of all faces if the vertex normals mostly point toward
    lower values of volume, so that normals point toward increasing values.

    points are in world coordinates, and volume is sampled through affine, so
    the result does not depend on the handedness of affine (a negative
    determinant mirrors the mesh, which reverses the normals of a mesh
    contoured in voxel space). step is the distance (mm) sampled either side
    of each vertex; by default the mean voxel size. Returns (faces, flipped).
    """
    if step is None:
        step = np.mean(np.linalg.norm(affine[:3, :3], axis=0))
    normals = vertex_normals(points, faces)
    normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1e-12
    to_voxel = np.linalg.inv(affine)

    def sample(p):
        ijk = np.floor(p @ to_voxel[:3, :3].T + to_voxel[:3, 3]).astype(int)
        inside = np.all((ijk >= 0) & (ijk < volume.shape), axis=1)
        values = np.full(len(p), np.nan)
        values[inside] = volume[tuple(ijk[inside].T)]
        return values

    diff = sample(points + step * normals) - sample(points - step * normals)
    flipped = bool(np.sum(diff < 0) > np.sum(diff > 0))
    if flipped:
        faces = faces[:, ::-1]
    return faces, flipped


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


# terminal labels, in the order of the distance columns
AP_SRC, AP_SINK, PD_SRC, PD_SINK = 0, 1, 2, 3
TERMINAL_NAMES = ["AP src", "AP sink", "PD src", "PD sink"]

# Valid cyclic orders of terminals around the boundary loop, starting at AP src.
# AP and PD terminals alternate, so each AP terminal is opposite the other, and
# likewise for PD. The two orders are mirror images (the loop can be walked
# either way, and hemispheres are mirrored), so both are allowed.
TERMINAL_ORDERS = (
    (AP_SRC, PD_SRC, AP_SINK, PD_SINK),
    (AP_SRC, PD_SINK, AP_SINK, PD_SRC),
)

# Expected order along the loop when walked in face-winding order, per
# hemisphere: the two hemispheres are mirror images. Used as a QC check, not a
# constraint, so a flipped surface or swapped terminals is reported.
EXPECTED_TERMINAL_ORDER = {"L": TERMINAL_ORDERS[0], "R": TERMINAL_ORDERS[1]}


def _mesh_edges(faces):
    """Directed edges of each triangle, in face winding order."""
    return np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])


def get_boundary_loop(faces, boundary_vertices):
    """
    Get the outer boundary loop of a mesh, as an ordered cycle of vertices.

    Only boundary edges (edges in exactly one face) with both vertices in
    boundary_vertices are used, and the largest loop is returned. The cycle
    follows the face winding, so its direction is consistent with the mesh
    orientation.
    """
    import networkx as nx

    directed = _mesh_edges(faces)
    undirected, counts = np.unique(
        np.sort(directed, axis=1), axis=0, return_counts=True
    )
    boundary_edges = undirected[counts == 1]
    in_set = np.isin(boundary_edges, boundary_vertices).all(axis=1)

    graph = nx.Graph()
    graph.add_edges_from(map(tuple, boundary_edges[in_set]))
    if graph.number_of_nodes() == 0:
        raise ValueError("No boundary edges found among boundary vertices")
    largest = graph.subgraph(max(nx.connected_components(graph), key=len))
    # a manifold boundary is a single cycle; otherwise take the longest one
    loop = max(nx.cycle_basis(largest), key=len)

    # orient the loop to follow the face winding
    directed_set = set(map(tuple, directed))
    if (loop[0], loop[1]) not in directed_set:
        loop = loop[::-1]
    return np.array(loop)


def _fit_cyclic_arcs(costs, order, nmin):
    """
    Split a cyclic sequence into contiguous arcs, one per terminal in order.

    costs is (m, 4), the cost of assigning each loop vertex to each terminal.
    Returns (total cost, labels) for the minimum cost split with every arc at
    least nmin long. Each split is found by dynamic programming over the
    sequence unrolled at each possible start of the first arc.
    """
    m = costs.shape[0]
    n_arcs = len(order)
    best_cost, best = np.inf, None

    for start in range(m):
        rolled = np.roll(costs, -start, axis=0)
        # prefix[k][i] is the cost of assigning the first i vertices to order[k]
        prefix = [np.concatenate([[0.0], np.cumsum(rolled[:, k])]) for k in order]

        # dp[j][b] is the min cost of arcs 0..j covering the first b vertices
        dp = [np.where(np.arange(m + 1) >= nmin, prefix[0], np.inf)]
        for j in range(1, n_arcs):
            prev_minus = dp[-1] - prefix[j]
            running_min = np.minimum.accumulate(prev_minus)
            cur = np.full(m + 1, np.inf)
            cur[nmin:] = prefix[j][nmin:] + running_min[: m + 1 - nmin]
            dp.append(cur)

        if dp[-1][m] < best_cost:
            best_cost, best = dp[-1][m], (start, dp, prefix)

    if best is None:
        return np.inf, None

    # backtrack the arc boundaries for the best start
    start, dp, prefix = best
    cuts = [m]
    for j in range(n_arcs - 1, 0, -1):
        end = cuts[-1]
        candidates = dp[j - 1][: end - nmin + 1] - prefix[j][: end - nmin + 1]
        cuts.append(int(np.argmin(candidates)))
    cuts = [0] + cuts[::-1]

    rolled_labels = np.empty(m, dtype=int)
    for j, label in enumerate(order):
        rolled_labels[cuts[j] : cuts[j + 1]] = label
    return best_cost, np.roll(rolled_labels, start)


def assign_boundary_terminals(
    distances, boundary_vertices, surface, nmin, ap_length_penalty=0.0, logger=None
):
    """
    Assign AP/PD src/sink terminal labels to boundary vertices.

    The terminals are fit to the outer boundary loop of the (topological disk)
    surface as four contiguous arcs, in alternating AP/PD order, so that the
    AP src and sink are on opposite sides, as are the PD src and sink. The arcs
    minimise the summed distance of each loop vertex to its assigned terminal,
    plus ap_length_penalty per vertex in the AP arcs, which favours short,
    point-like AP terminals. Every arc has at least nmin vertices.

    The remaining boundary vertices (the band of interior vertices adjacent to
    the loop) take the label of their nearest loop vertex, so each label is
    contiguous.

    Parameters:
        distances (np.ndarray): (n, 4) signed distance from each boundary vertex
            to the AP src, AP sink, PD src and PD sink.
        boundary_vertices (np.ndarray): (n,) mesh indices of the boundary vertices.
        surface: mesh with pyvista-style faces.
        nmin (int): minimum number of loop vertices per terminal.
        ap_length_penalty (float): cost per vertex of AP src/sink arcs, in the
            units of distances.

    Returns:
        labels (np.ndarray): (n,) terminal labels, 0-3, for boundary_vertices.
        info (dict): diagnostics, with the boundary loop (mesh indices, in
            face-winding order), the fitted and nearest-terminal labels along
            it, the fitted terminal order, and the number of band vertices not
            connected to the loop.
    """
    from collections import deque

    faces = surface.faces.reshape((-1, 4))[:, 1:4]
    vertex_to_idx = {v: i for i, v in enumerate(boundary_vertices)}

    loop = get_boundary_loop(faces, boundary_vertices)
    loop_idx = np.array([vertex_to_idx[v] for v in loop])
    m = len(loop)
    if m < 4 * nmin:
        raise ValueError(
            f"Boundary loop has {m} vertices, too few for 4 terminals of at least "
            f"{nmin} vertices"
        )

    costs = np.asarray(distances, dtype=float)[loop_idx].copy()
    costs[:, [AP_SRC, AP_SINK]] += ap_length_penalty

    fits = [_fit_cyclic_arcs(costs, order, nmin) for order in TERMINAL_ORDERS]
    best = int(np.argmin([cost for cost, _ in fits]))
    loop_labels = fits[best][1]
    unconstrained = np.argmin(distances[loop_idx], axis=1)

    if logger:
        names = " -> ".join(TERMINAL_NAMES[k] for k in TERMINAL_ORDERS[best])
        logger.info(
            f"Boundary loop has {m} of {len(boundary_vertices)} boundary vertices"
        )
        logger.info(f"Terminal order along loop (face winding): {names}")
        logger.info(
            "Loop vertices per terminal: "
            + ", ".join(
                f"{TERMINAL_NAMES[k]}={np.sum(loop_labels == k)}" for k in range(4)
            )
        )
        logger.info(
            f"{np.sum(loop_labels != unconstrained)} loop vertices differ from the "
            "nearest-terminal labels"
        )

    # label the rest of the band by nearest loop vertex (multi-source BFS)
    neighbors = [[] for _ in boundary_vertices]
    for u, v in _mesh_edges(faces):
        if u in vertex_to_idx and v in vertex_to_idx:
            neighbors[vertex_to_idx[u]].append(vertex_to_idx[v])
            neighbors[vertex_to_idx[v]].append(vertex_to_idx[u])

    labels = np.full(len(boundary_vertices), -1, dtype=int)
    labels[loop_idx] = loop_labels
    queue = deque(loop_idx)
    while queue:
        i = queue.popleft()
        for j in neighbors[i]:
            if labels[j] < 0:
                labels[j] = labels[i]
                queue.append(j)

    unreached = labels < 0
    if np.any(unreached):
        # not connected to the loop, so fall back to the nearest terminal
        if logger:
            logger.warning(
                f"{np.sum(unreached)} boundary vertices are not connected to the "
                "boundary loop, assigning nearest terminal"
            )
        labels[unreached] = np.argmin(distances[unreached], axis=1)

    info = {
        "loop": loop,
        "loop_labels": loop_labels,
        "loop_nearest_labels": unconstrained,
        "order": TERMINAL_ORDERS[best],
        "n_unreached": int(np.sum(unreached)),
    }
    return labels, info


def isotonic_increasing(y):
    """Least-squares non-decreasing fit to y (pool adjacent violators)."""
    values, weights = [], []
    for yi in np.asarray(y, dtype=float):
        values.append(yi)
        weights.append(1)
        while len(values) > 1 and values[-2] > values[-1]:
            w = weights[-2] + weights[-1]
            v = (values[-2] * weights[-2] + values[-1] * weights[-1]) / w
            values[-2:] = [v]
            weights[-2:] = [w]
    return np.repeat(values, weights)


def monotone_free_arcs(loop, values, src, sink, eps=1e-9):
    """
    Values along the boundary loop, made strictly monotone on each free arc
    (the vertices between the src and sink arcs, which are not fixed in the
    Laplace-Beltrami solve), increasing from the src end to the sink end.

    With non-negative weights, the unfolded (AP, PD) map cannot fold if the
    whole boundary maps monotonically onto the boundary of the unit square;
    the free arcs of each coordinate are where it can fail to.

    Parameters:
        loop (np.ndarray): boundary loop vertices, in order.
        values (np.ndarray): coordinate at each mesh vertex.
        src, sink (np.ndarray): boolean masks of the src and sink vertices.
        eps (float): step between values pooled by the isotonic fit, to make
            them strictly monotone.

    Returns:
        loop_values (np.ndarray): corrected values at the loop vertices.
        n_changed (int): number of loop vertices changed (pooled) by the fit.
    """
    loop = np.asarray(loop)
    loop_values = np.asarray(values, dtype=float)[loop].copy()
    fixed = src[loop] | sink[loop]
    if fixed.all() or not fixed.any():
        return loop_values, 0
    m = len(loop)
    # walk the loop from just after a fixed vertex, collecting free arcs
    start = int(np.flatnonzero(fixed & ~np.roll(fixed, -1))[0]) + 1
    arcs, arc = [], []
    for i in np.roll(np.arange(m), -start):
        if fixed[i]:
            if arc:
                arcs.append(arc)
            arc = []
        else:
            arc.append(i)
    if arc:
        arcs.append(arc)

    n_changed = 0
    for arc in arcs:
        before, after = loop[(arc[0] - 1) % m], loop[(arc[-1] + 1) % m]
        if src[before] and sink[after]:
            sign = 1.0
        elif sink[before] and src[after]:
            sign = -1.0
        else:
            continue  # both ends on the same terminal: no expected direction
        y = loop_values[arc]
        fit = sign * isotonic_increasing(sign * y)
        # break the ties of pooled blocks, so the arc is strictly monotone
        starts = np.flatnonzero(np.r_[True, fit[1:] != fit[:-1]])
        ends = np.r_[starts[1:], len(fit)]
        for a, b in zip(starts, ends):
            if b - a > 1:
                fit[a:b] += sign * eps * (np.arange(b - a) - (b - a - 1) / 2)
                n_changed += b - a
        loop_values[arc] = np.clip(fit, eps, 1 - eps)
    return loop_values, n_changed
