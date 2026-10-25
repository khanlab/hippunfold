import numpy as np
from lib.qc_surface import boundary_chords, compute_topology, face_geometry
from lib.surface_quality import (
    SurfaceProjector,
    closest_points_on_triangles,
    flip_edges,
    improve_quality,
    split_boundary_chords,
)
from test_qc_surface import disk, grid


def capped_disk(n=8, seed=0):
    """Disk with jittered vertices and every other diagonal skewed into caps."""
    u, v, faces = grid(n, n)
    rng = np.random.default_rng(seed)
    interior = (u > 0) & (u < 1) & (v > 0) & (v < 1)
    u = u + interior * rng.uniform(-0.45, 0.45, len(u)) / n
    v = v + interior * rng.uniform(-0.45, 0.45, len(v)) / n
    points = np.c_[u * 10, v * 10, 0.5 * np.sin(2 * u)]
    return points, faces


def test_closest_points_brute_force():
    rng = np.random.default_rng(1)
    a, b, c = rng.normal(size=(3, 500, 3))
    p = rng.normal(size=(500, 3)) * 2
    q = closest_points_on_triangles(p, a, b, c)
    s, t = np.meshgrid(np.linspace(0, 1, 101), np.linspace(0, 1, 101))
    keep = s + t <= 1
    s, t = s[keep], t[keep]
    samples = (
        a[:, None]
        + (b - a)[:, None] * s[None, :, None]
        + (c - a)[:, None] * t[None, :, None]
    )
    brute = np.linalg.norm(samples - p[:, None], axis=2).min(axis=1)
    exact = np.linalg.norm(q - p, axis=1)
    assert np.all(exact <= brute + 1e-9)
    assert np.all(brute - exact < 0.05)


def test_flip_removes_cap():
    # move an interior vertex next to the diagonal of its cell, making a cap
    points, faces, u, v = disk()
    points = points.copy()
    a = np.flatnonzero(np.isclose(u, 3 / 8) & np.isclose(v, 3 / 8))[0]
    c = np.flatnonzero(np.isclose(u, 4 / 8) & np.isclose(v, 4 / 8))[0]
    d = np.flatnonzero(np.isclose(u, 3 / 8) & np.isclose(v, 4 / 8))[0]
    points[d] = 0.55 * points[a] + 0.45 * points[c] + [-0.05, 0.05, 0]
    before = face_geometry(points, faces)[1].min()
    new, n = flip_edges(points, faces, max_flip_dist=0.1)
    assert n >= 1
    assert face_geometry(points, new)[1].min() > before
    assert compute_topology(points, new)["misoriented_faces"].size == 0


def test_flip_never_makes_chords():
    # flipping this quad would join its two boundary vertices
    points = np.array([[0, 0, 0], [2, 0, 0], [1, 0.1, 0], [1, -0.1, 0]], float)
    faces = np.array([[0, 1, 2], [1, 0, 3]])
    _, n = flip_edges(points, faces, max_flip_dist=0.01)
    assert n == 0


def test_flip_respects_planarity():
    points, faces, u, v = disk()
    points = points.copy()
    a = np.flatnonzero(np.isclose(u, 3 / 8) & np.isclose(v, 3 / 8))[0]
    c = np.flatnonzero(np.isclose(u, 4 / 8) & np.isclose(v, 4 / 8))[0]
    d = np.flatnonzero(np.isclose(u, 3 / 8) & np.isclose(v, 4 / 8))[0]
    # same cap, but lifted out of the surface: flipping would change the shape
    points[d] = 0.55 * points[a] + 0.45 * points[c] + [-0.05, 0.05, 2]
    new, _ = flip_edges(points, faces, max_flip_dist=0.1)
    # the cap's diagonal (a, c) is kept (only the two corner chords flip)
    edges = {
        tuple(sorted(e))
        for f in new
        for e in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0]))
    }
    assert tuple(sorted((a, c))) in edges


def test_split_boundary_chords():
    points, faces, _, _ = disk()
    chords = boundary_chords(faces)
    assert len(chords) == 2  # the grid diagonals at two corners make ears
    new_points, new_faces, n = split_boundary_chords(points, faces)
    assert n == 2
    assert len(new_points) == len(points) + 2 and len(new_faces) == len(faces) + 4
    assert len(boundary_chords(new_faces)) == 0
    t0, t1 = compute_topology(points, faces), compute_topology(new_points, new_faces)
    for k in ("chi", "components", "boundary_loops", "genus"):
        assert t0[k] == t1[k]
    assert len(t1["misoriented_faces"]) == 0
    # midpoints lie on the chords, so the surface is unchanged
    assert np.allclose(
        new_points[-2:],
        points[chords].mean(axis=1)[
            np.argsort(
                [
                    np.argmin(np.linalg.norm(points[chords].mean(axis=1) - q, axis=1))
                    for q in new_points[-2:]
                ]
            )
        ],
    )
    assert np.isclose(
        face_geometry(new_points, new_faces)[0].sum(),
        face_geometry(points, faces)[0].sum(),
    )


def test_split_face_with_two_chords():
    # pentagon fanned into 3 faces: the middle face (a, b, c) has all its
    # vertices on the boundary and two interior edges, both chords
    points = np.array(
        [[0, 0, 0], [2, 0, 0], [2, 1.5, 0], [1, 2, 0], [0, 1.5, 0]], float
    )
    a, b, y, c, x = range(5)
    faces = np.array([[a, b, c], [b, y, c], [c, x, a]])
    assert len(boundary_chords(faces)) == 2
    new_points, new_faces, n = split_boundary_chords(points, faces)
    assert n == 2  # one per pass, as both chords share the middle face
    assert len(boundary_chords(new_faces)) == 0
    t = compute_topology(new_points, new_faces)
    assert t["chi"] == 1 and t["boundary_loops"] == 1
    assert len(t["misoriented_faces"]) == 0


def test_improve_quality_keeps_topology_boundary_and_shape():
    points, faces = capped_disk()
    t0 = compute_topology(points, faces)
    new_points, new_faces, info = improve_quality(points, faces)
    t1 = compute_topology(new_points, new_faces)
    for k in ("chi", "components", "boundary_loops", "genus"):
        assert t0[k] == t1[k]
    assert t1["V"] == t0["V"] + info["boundary_chords_split"]
    assert len(t1["misoriented_faces"]) == 0
    assert len(boundary_chords(new_faces)) == 0
    # boundary vertices fixed, all vertices on the original surface
    boundary = np.unique(t0["loops"][0])
    assert np.allclose(new_points[boundary], points[boundary])
    on_surface = SurfaceProjector(points, faces)(new_points)
    assert np.allclose(on_surface, new_points, atol=1e-9)
    assert info["max_surface_distance"] < 0.1
    # quality improves
    before = face_geometry(points, faces)[1]
    after = face_geometry(new_points, new_faces)[1]
    assert np.percentile(after, 1) > np.percentile(before, 1)
    assert info["edge_flips"] > 0


def test_improve_quality_on_good_mesh_is_gentle():
    points, faces, _, _ = disk()
    new_points, new_faces, info = improve_quality(points, faces)
    assert info["max_surface_distance"] < 0.05
    # the two corner ears are fixed by flips; angles stay close to the grid's
    assert info["boundary_chords_flipped"] == 2
    assert info["boundary_chords_split"] == 0
    assert face_geometry(new_points, new_faces)[1].min() > 30


def test_closest_points_degenerate_triangle():
    a = np.array([[0.0, 0, 0], [0, 0, 0]])
    b = np.array([[1.0, 0, 0], [0, 0, 0]])
    c = np.array([[2.0, 0, 0], [0, 0, 0]])  # collinear, and a single point
    p = np.array([[1.9, 1, 0], [1, 1, 1]])
    q = closest_points_on_triangles(p, a, b, c)
    assert np.all(np.isfinite(q))
    assert np.allclose(q[1], 0)
    assert np.linalg.norm(q[0] - p[0]) <= np.linalg.norm(c[0] - p[0]) + 1e-12
