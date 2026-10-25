"""
Triangle-quality improvement of a surface mesh, keeping its topology and
shape.

- flip_edges: flips interior edges where that raises the smallest angle of
  the two adjoining triangles, if the quad they form is nearly flat (so the
  surface hardly changes). This removes caps (triangles with an angle near
  180 degrees), which make up most slivers of the decimated isosurface.
- relax_vertices: moves interior vertices towards the area-weighted centroid
  of their faces, within the tangent plane, then back onto the original
  surface. Boundary vertices are fixed.

- boundary chords, interior edges that join two boundary vertices (e.g.
  the inner edge of an ear, a face with two boundary edges), are flipped
  where possible, and otherwise split at their midpoint by
  split_boundary_chords. Otherwise a face with all three vertices on one
  terminal arc of the boundary has a fixed coordinate at every vertex, and
  collapses to zero area in the unfolded space.

improve_quality removes boundary chords (by flips, or else splits), then
alternates flips and relaxation, and reports how far vertices moved.
"""

import numpy as np
from lib.qc_surface import boundary_chords, face_geometry, unique_edges
from scipy.spatial import cKDTree


def _face_normals(points, faces):
    tri = points[faces]
    return np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])


def _unit(v):
    return v / (np.linalg.norm(v, axis=-1, keepdims=True) + 1e-30)


def _min_angle(a, b, c):
    """Smallest angle (degrees) of triangles (a, b, c), vectorised."""
    angles = []
    for p, q, r in ((a, b, c), (b, c, a), (c, a, b)):
        u, v = _unit(q - p), _unit(r - p)
        angles.append(np.degrees(np.arccos(np.clip(np.sum(u * v, axis=-1), -1, 1))))
    return np.min(angles, axis=0)


def flip_edges(
    points, faces, max_flip_dist, min_gain=1.0, min_chord_angle=10.0, max_passes=20
):
    """
    Flip interior edges that raise the smallest angle of their two triangles
    by at least min_gain degrees, and boundary chords (interior edges joining
    two boundary vertices), first and whenever the new triangles' smallest
    angle is at least min_chord_angle. A flip is allowed only if the quad is
    nearly flat (the two diagonals are at most max_flip_dist apart), both new
    faces keep the orientation of the old ones, the new edge does not exist
    yet and is not a chord, and no vertex is left with fewer than 3 edges (2
    on the boundary). Returns the new faces and the number of flips.
    """
    faces = np.asarray(faces, dtype=np.int64).copy()
    n_points = len(points)
    total = 0
    for _ in range(max_passes):
        edges, inverse, directed = unique_edges(faces)
        counts = np.bincount(inverse, minlength=len(edges))
        interior = counts[inverse] == 2
        order = np.argsort(inverse[interior], kind="stable")
        pairs = np.flatnonzero(interior)[order].reshape(-1, 2)  # face-edge ids
        f1, f2 = pairs[:, 0] // 3, pairs[:, 1] // 3
        a, b = directed[pairs[:, 0], 0], directed[pairs[:, 0], 1]
        c = faces[f1, (pairs[:, 0] % 3 + 2) % 3]
        d = faces[f2, (pairs[:, 1] % 3 + 2) % 3]

        # the two faces must be wound consistently (edge in opposite directions)
        consistent = directed[pairs[:, 1], 0] == b

        pa, pb, pc, pd = points[a], points[b], points[c], points[d]
        before = np.minimum(_min_angle(pa, pb, pc), _min_angle(pb, pa, pd))
        # new faces (a, d, c) and (b, c, d) keep the winding of (a, b, c), (b, a, d)
        after = np.minimum(_min_angle(pa, pd, pc), _min_angle(pb, pc, pd))
        gain = after - before

        normal = _face_normals(points, faces[f1]) + _face_normals(points, faces[f2])
        n1 = _unit(np.cross(pd - pa, pc - pa))
        n2 = _unit(np.cross(pc - pb, pd - pb))
        keeps_orientation = (np.sum(n1 * normal, axis=1) > 0) & (
            np.sum(n2 * normal, axis=1) > 0
        )
        cross = np.cross(pb - pa, pd - pc)
        dist = np.abs(np.sum(cross * (pc - pa), axis=1)) / (
            np.linalg.norm(cross, axis=1) + 1e-30
        )

        boundary = np.zeros(n_points, dtype=bool)
        boundary[edges[counts == 1].ravel()] = True
        degree = np.bincount(edges.ravel(), minlength=n_points)
        min_degree = np.where(boundary, 2, 3)
        ok_degree = (degree[a] - 1 >= min_degree[a]) & (degree[b] - 1 >= min_degree[b])

        # never join two boundary vertices by an interior edge (a chord), but
        # remove chords where possible
        no_chord = ~(boundary[c] & boundary[d])
        is_chord = boundary[a] & boundary[b]
        wanted = np.where(is_chord, after >= min_chord_angle, gain >= min_gain)
        priority = np.where(is_chord, np.inf, gain)

        candidates = np.flatnonzero(
            wanted
            & consistent
            & no_chord
            & keeps_orientation
            & (dist <= max_flip_dist)
            & ok_degree
        )
        if len(candidates) == 0:
            break
        existing = set((edges[:, 0] * n_points + edges[:, 1]).tolist())
        touched = np.zeros(len(faces), dtype=bool)
        n_flips = 0
        for i in candidates[np.argsort(-priority[candidates], kind="stable")]:
            if touched[f1[i]] or touched[f2[i]]:
                continue
            key = min(c[i], d[i]) * n_points + max(c[i], d[i])
            if key in existing:
                continue
            faces[f1[i]] = (a[i], d[i], c[i])
            faces[f2[i]] = (b[i], c[i], d[i])
            touched[[f1[i], f2[i]]] = True
            existing.add(key)
            n_flips += 1
        total += n_flips
        if n_flips == 0:
            break
    return faces, total


def split_boundary_chords(points, faces, max_passes=10):
    """
    Split interior edges joining two boundary vertices at their midpoint, so
    every face has a vertex off the boundary. Each split replaces the two
    faces of the edge by four, keeping their winding. Returns points, faces
    and the number of splits.
    """
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=np.int64).copy()
    total = 0
    for _ in range(max_passes):
        edges, inverse, directed = unique_edges(faces)
        counts = np.bincount(inverse, minlength=len(edges))
        boundary = np.zeros(len(points), dtype=bool)
        boundary[edges[counts == 1].ravel()] = True
        chord = (counts == 2) & boundary[edges[:, 0]] & boundary[edges[:, 1]]
        if not chord.any():
            break
        face_edges = np.flatnonzero(chord[inverse])
        touched = np.zeros(len(faces), dtype=bool)
        new_points, replace, extra = [], {}, []
        for e in np.unique(inverse[face_edges]):
            fe1, fe2 = face_edges[inverse[face_edges] == e]
            f1, f2 = fe1 // 3, fe2 // 3
            if touched[f1] or touched[f2]:
                continue  # a face with two chords: split the other next pass
            a, b = directed[fe1]
            c = faces[f1, (fe1 % 3 + 2) % 3]
            d = faces[f2, (fe2 % 3 + 2) % 3]
            m = len(points) + len(new_points)
            new_points.append(0.5 * (points[a] + points[b]))
            # (a, b, c) -> (a, m, c), (m, b, c); (b, a, d) -> (b, m, d), (m, a, d)
            replace[f1] = (a, m, c)
            replace[f2] = (b, m, d)
            extra += [(m, b, c), (m, a, d)]
            touched[[f1, f2]] = True
        for f, tri in replace.items():
            faces[f] = tri
        points = np.vstack([points, new_points])
        faces = np.vstack([faces, extra])
        total += len(new_points)
    return points, faces, total


def closest_points_on_triangles(p, a, b, c):
    """Closest point to each p on triangle (a, b, c), vectorised (Ericson)."""
    p, a, b, c = np.broadcast_arrays(p, a, b, c)
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = np.sum(ab * ap, -1), np.sum(ac * ap, -1)
    bp = p - b
    d3, d4 = np.sum(ab * bp, -1), np.sum(ac * bp, -1)
    cp = p - c
    d5, d6 = np.sum(ab * cp, -1), np.sum(ac * cp, -1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    denom = va + vb + vc
    with np.errstate(divide="ignore", invalid="ignore"):
        v = vb / denom
        w = vc / denom
        result = a + ab * v[..., None] + ac * w[..., None]  # interior

    def where(mask, value):
        result[mask] = value[mask]

    # regions in reverse order of Ericson's tests, so the first test that
    # applies is written last
    with np.errstate(divide="ignore", invalid="ignore"):
        t_bc = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        where(
            (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0), b + (c - b) * t_bc[..., None]
        )
        t_ac = d2 / (d2 - d6)
        where((vb <= 0) & (d2 >= 0) & (d6 <= 0), a + ac * t_ac[..., None])
        where((d6 >= 0) & (d5 <= d6), c)
        t_ab = d1 / (d1 - d3)
        where((vc <= 0) & (d1 >= 0) & (d3 <= 0), a + ab * t_ab[..., None])
    where((d3 >= 0) & (d4 <= d3), b)
    where((d1 <= 0) & (d2 <= 0), a)
    # degenerate (zero-area) triangles: the nearest vertex
    bad = ~np.all(np.isfinite(result), axis=-1)
    if bad.any():
        corners = np.stack(np.broadcast_arrays(a, b, c), axis=-2)
        nearest = np.argmin(
            np.linalg.norm(
                corners - np.broadcast_to(p, a.shape)[..., None, :], axis=-1
            ),
            axis=-1,
        )
        vertex = np.take_along_axis(corners, nearest[..., None, None], axis=-2)[
            ..., 0, :
        ]
        result[bad] = vertex[bad]
    return result


class SurfaceProjector:
    """Projects points onto the nearest point of a fixed reference surface."""

    def __init__(self, points, faces, k=12):
        self.tri = np.asarray(points, dtype=float)[faces]
        self.tree = cKDTree(self.tri.mean(axis=1))
        self.k = min(k, len(faces))

    def __call__(self, p):
        _, idx = self.tree.query(p, k=self.k)
        idx = idx.reshape(len(p), -1)
        t = self.tri[idx]  # (n, k, 3, 3)
        q = closest_points_on_triangles(
            p[:, None, :], t[:, :, 0], t[:, :, 1], t[:, :, 2]
        )
        best = np.argmin(np.linalg.norm(q - p[:, None, :], axis=2), axis=1)
        return q[np.arange(len(p)), best]


def _worse_faces(points, faces, normals, min_angle):
    """Faces flipped relative to normals, or with an angle below min_angle."""
    flipped = np.sum(_face_normals(points, faces) * normals, axis=1) <= 0
    return flipped | (face_geometry(points, faces)[1] < min_angle)


def relax_vertices(points, faces, project, iterations=3, step=0.5, min_angle=5.0):
    """
    Tangential relaxation of interior vertices towards the area-weighted
    centroid of their faces, projected back onto the reference surface.
    Moves that would flip a face, or make its smallest angle drop below both
    its current value and min_angle (degrees, the sliver threshold of the QC),
    are undone.
    """
    points = np.asarray(points, dtype=float).copy()
    n = len(points)
    edges, inverse, _ = unique_edges(faces)
    counts = np.bincount(inverse, minlength=len(edges))
    interior = np.ones(n, dtype=bool)
    interior[edges[counts == 1].ravel()] = False
    corner_vertex = faces.ravel()

    for _ in range(iterations):
        fn = _face_normals(points, faces)
        area = 0.5 * np.linalg.norm(fn, axis=1)
        centroid = points[faces].mean(axis=1)
        weight = np.repeat(area, 3)
        total = np.bincount(corner_vertex, weights=weight, minlength=n)
        target = np.stack(
            [
                np.bincount(
                    corner_vertex,
                    weights=weight * np.repeat(centroid[:, k], 3),
                    minlength=n,
                )
                for k in range(3)
            ],
            axis=1,
        ) / (total[:, None] + 1e-30)
        vn = _unit(
            np.stack(
                [
                    np.bincount(
                        corner_vertex, weights=np.repeat(fn[:, k], 3), minlength=n
                    )
                    for k in range(3)
                ],
                axis=1,
            )
        )
        move = target - points
        move -= np.sum(move * vn, axis=1, keepdims=True) * vn
        new = points.copy()
        new[interior] = project(points[interior] + step * move[interior])

        # undo moves around faces that would flip, or whose smallest angle
        # would drop below both its current value and min_angle (e.g. a vertex
        # next to the boundary projected onto a boundary edge)
        floor = np.minimum(face_geometry(points, faces)[1], min_angle)
        bad = _worse_faces(new, faces, fn, floor)
        while bad.any():
            revert = np.unique(faces[bad])
            new[revert] = points[revert]
            bad = _worse_faces(new, faces, fn, floor)
        points = new
    return points


def improve_quality(
    points, faces, rounds=3, relax_iterations=3, step=0.5, flip_tolerance=0.1
):
    """
    Alternate edge flips and tangential relaxation. flip_tolerance is the
    maximum distance between the diagonals of a flipped quad, as a fraction of
    the median edge length. Returns points, faces and a summary.
    """
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    original_points = points.copy()
    project = SurfaceProjector(points, faces)
    median_edge = float(np.median(face_geometry(points, faces)[2]))
    max_flip_dist = flip_tolerance * median_edge

    # remove boundary chords by flipping where possible, else by splitting
    # (midpoints of split edges lie on the surface, so need no projection)
    n_chords = len(boundary_chords(faces))
    faces, n_flips = flip_edges(points, faces, max_flip_dist)
    n_chords_flipped = n_chords - len(boundary_chords(faces))
    points, faces, n_splits = split_boundary_chords(points, faces)
    original = points.copy()

    for _ in range(rounds):
        faces, n = flip_edges(points, faces, max_flip_dist)
        n_flips += n
        points = relax_vertices(points, faces, project, relax_iterations, step)
    faces, n = flip_edges(points, faces, max_flip_dist)
    n_flips += n

    moved = np.linalg.norm(points[: len(original)] - original, axis=1)
    # vertices stay on the original surface; faces are chords of it, so measure
    # how far face centroids, and original vertices, are from the other surface
    tri = points[faces]
    centroids = tri.mean(axis=1)
    to_original = np.linalg.norm(project(centroids) - centroids, axis=1)
    to_new = np.linalg.norm(
        SurfaceProjector(points, faces)(original_points) - original_points, axis=1
    )
    return (
        points,
        faces,
        {
            "boundary_chords_flipped": n_chords_flipped,
            "boundary_chords_split": n_splits,
            "edge_flips": n_flips,
            "max_vertex_move": float(moved.max()),
            "mean_vertex_move": float(moved.mean()),
            "max_surface_distance": float(max(to_original.max(), to_new.max())),
            "mean_surface_distance": float(np.r_[to_original, to_new].mean()),
            "median_edge": median_edge,
        },
    )
