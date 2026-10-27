"""
Repair of the native isosurface mesh, so that it is a clean disk for the
Laplace-Beltrami solve.

- fill_holes closes every boundary loop except the outer boundary (up to a
  size limit), with a fan of triangles to the loop centroid. Each new face
  reverses the boundary edge of the face it adjoins, so the winding is
  consistent by construction (vtkFillHolesFilter can produce misoriented and
  non-manifold faces).
- repair_mesh removes defects while leaving the remaining geometry unchanged:
  faces with repeated vertices, duplicate faces, extra faces on non-manifold
  edges, small components; it splits non-manifold (bowtie) vertices, flips
  misoriented faces to the majority winding and removes unreferenced vertices.

Topology helpers are shared with the surface QC (lib/qc_surface.py).
"""

import numpy as np
from lib.qc_surface import (
    _components,
    corner_groups,
    face_geometry,
    orientation,
    unique_edges,
)


def _edge_data(faces):
    edges, inverse, directed = unique_edges(faces)
    edge_counts = np.bincount(inverse, minlength=len(edges))
    return edges, inverse, directed, edge_counts


def fill_holes(points, faces, max_radius):
    """
    Fill holes: boundary loops other than the outer one (longest perimeter)
    whose radius (max distance of a loop vertex to the loop centroid) is at
    most max_radius. Returns points, faces and the radius of each filled hole.
    """
    faces = np.asarray(faces, dtype=np.int64)
    points = np.asarray(points, dtype=float)
    _, inverse, directed, edge_counts = _edge_data(faces)
    bdir = directed[edge_counts[inverse] == 1]  # boundary edges, face winding
    if len(bdir) == 0:
        return points, faces, []

    verts, local = np.unique(bdir, return_inverse=True)
    local = local.reshape(-1, 2)
    n_loops, comp = _components(len(verts), local[:, 0], local[:, 1])
    edge_loop = comp[local[:, 0]]
    lengths = np.linalg.norm(points[bdir[:, 0]] - points[bdir[:, 1]], axis=1)
    perimeter = np.bincount(edge_loop, weights=lengths, minlength=n_loops)
    outer = int(np.argmax(perimeter))

    new_points, new_faces, radii = [], [], []
    for k in range(n_loops):
        if k == outer:
            continue
        loop_edges = bdir[edge_loop == k]
        loop_verts = np.unique(loop_edges)
        centroid = points[loop_verts].mean(axis=0)
        radius = float(np.linalg.norm(points[loop_verts] - centroid, axis=1).max())
        if radius > max_radius:
            continue
        c = len(points) + len(new_points)
        new_points.append(centroid)
        # the face on edge (a, b) winds a -> b, so the fill face winds b -> a
        new_faces.append(
            np.c_[loop_edges[:, 1], loop_edges[:, 0], np.full(len(loop_edges), c)]
        )
        radii.append(radius)
    if not new_faces:
        return points, faces, []
    return (
        np.vstack([points, new_points]),
        np.vstack([faces, *new_faces]),
        radii,
    )


def remove_bad_faces(points, faces):
    """Drop faces with repeated vertices, and duplicates of a face."""
    repeated = (faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2])
    repeated |= faces[:, 0] == faces[:, 2]
    faces = faces[~repeated]
    _, first = np.unique(np.sort(faces, axis=1), axis=0, return_index=True)
    n_dup = len(faces) - len(first)
    return faces[np.sort(first)], {
        "repeated_vertex_faces": int(repeated.sum()),
        "duplicate_faces": n_dup,
    }


def remove_nonmanifold_edge_faces(points, faces):
    """
    On edges with more than two faces, keep two faces and drop the others.
    Dangling faces (fins, flaps), whose other edges are on the boundary, are
    dropped first, then the smallest.
    """
    _, inverse, _, edge_counts = _edge_data(faces)
    bad_edges = np.flatnonzero(edge_counts > 2)
    if len(bad_edges) == 0:
        return faces, 0
    area = face_geometry(points, faces)[0]
    n_boundary = (edge_counts[inverse] == 1).reshape(-1, 3).sum(axis=1)
    face_of = np.repeat(np.arange(len(faces)), 3)
    drop = set()
    for e in bad_edges:
        incident = face_of[inverse == e]
        incident = incident[np.lexsort((-area[incident], n_boundary[incident]))]
        drop.update(incident[2:].tolist())
    keep = np.setdiff1d(np.arange(len(faces)), list(drop))
    return faces[keep], len(drop)


def split_nonmanifold_vertices(points, faces):
    """
    Give each extra fan of faces around a non-manifold (bowtie) vertex its own
    copy of the vertex.
    """
    _, inverse, _, edge_counts = _edge_data(faces)
    comp = corner_groups(faces, inverse, edge_counts)
    corners = faces.ravel()
    groups, group_idx = np.unique(np.c_[corners, comp], axis=0, return_inverse=True)
    group_idx = group_idx.ravel()
    # the first group of each vertex keeps the vertex, the others get copies
    first = np.r_[True, groups[1:, 0] != groups[:-1, 0]]
    if first.all():
        return points, faces, 0
    new_vertex = groups[:, 0].copy()
    extra = np.flatnonzero(~first)
    new_vertex[extra] = len(points) + np.arange(len(extra))
    points = np.vstack([points, points[groups[extra, 0]]])
    faces = new_vertex[group_idx].reshape(-1, 3)
    return points, faces, len(np.unique(groups[extra, 0]))


def keep_largest_component(points, faces):
    """Keep the faces of the largest connected component (by faces)."""
    # link each corner to the next one in its face
    _, comp = _components(len(points), faces.ravel(), np.roll(faces, 1, axis=1).ravel())
    face_comp = comp[faces[:, 0]]
    sizes = np.bincount(face_comp)
    keep = face_comp == np.argmax(sizes)
    return faces[keep], int((~keep).sum())


def orient_faces(points, faces):
    """Flip faces wound opposite to the majority of their component."""
    edges, inverse, directed, edge_counts = _edge_data(faces)
    orientable, _, misoriented = orientation(
        faces, edges, inverse, directed, edge_counts
    )
    if not orientable or len(misoriented) == 0:
        return faces, 0, orientable
    faces = faces.copy()
    faces[misoriented] = faces[misoriented, ::-1]
    return faces, len(misoriented), orientable


def remove_unreferenced_vertices(points, faces):
    used = np.unique(faces)
    if len(used) == len(points):
        return points, faces, 0
    remap = np.full(len(points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return points[used], remap[faces], len(points) - len(used)


def repair_mesh(points, faces, max_iterations=5):
    """
    Remove topological defects, repeating until the mesh no longer changes.
    Returns points, faces and counts of what was changed.
    """
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    counts = {
        "repeated_vertex_faces": 0,
        "duplicate_faces": 0,
        "nonmanifold_edge_faces": 0,
        "nonmanifold_vertices_split": 0,
        "small_component_faces": 0,
        "misoriented_faces_flipped": 0,
        "unreferenced_vertices": 0,
    }
    orientable = True
    for _ in range(max_iterations):
        changed = 0
        faces, c = remove_bad_faces(points, faces)
        for k, v in c.items():
            counts[k] += v
            changed += v
        faces, n = remove_nonmanifold_edge_faces(points, faces)
        counts["nonmanifold_edge_faces"] += n
        changed += n
        points, faces, n = split_nonmanifold_vertices(points, faces)
        counts["nonmanifold_vertices_split"] += n
        changed += n
        faces, n = keep_largest_component(points, faces)
        counts["small_component_faces"] += n
        changed += n
        faces, n, orientable = orient_faces(points, faces)
        counts["misoriented_faces_flipped"] += n
        changed += n
        points, faces, n = remove_unreferenced_vertices(points, faces)
        counts["unreferenced_vertices"] += n
        if not changed:
            break
    counts["orientable"] = orientable
    return points, faces, counts
