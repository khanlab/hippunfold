import numpy as np
from lib.qc_surface import compute_topology, face_geometry
from lib.surface_quality import (
    SurfaceProjector,
    closest_points_on_triangles,
    flip_edges,
    improve_quality,
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
    # square split along the long diagonal of a flattened quad gives a cap
    points = np.array([[0, 0, 0], [2, 0, 0], [1, 0.1, 0], [1, -0.1, 0]], float)
    faces = np.array([[0, 1, 2], [1, 0, 3]])
    before = face_geometry(points, faces)[1].min()
    new, n = flip_edges(points, faces, max_flip_dist=0.01)
    assert n == 1
    assert face_geometry(points, new)[1].min() > before
    assert compute_topology(points, new)["misoriented_faces"].size == 0


def test_flip_respects_planarity():
    # same quad, folded out of plane: flipping would change the shape
    points = np.array([[0, 0, 0], [2, 0, 0], [1, 0.1, 0.5], [1, -0.1, 0.5]], float)
    faces = np.array([[0, 1, 2], [1, 0, 3]])
    _, n = flip_edges(points, faces, max_flip_dist=0.01)
    assert n == 0


def test_improve_quality_keeps_topology_boundary_and_shape():
    points, faces = capped_disk()
    t0 = compute_topology(points, faces)
    new_points, new_faces, info = improve_quality(points, faces)
    t1 = compute_topology(new_points, new_faces)
    for k in ("V", "F", "chi", "components", "boundary_loops", "genus"):
        assert t0[k] == t1[k]
    assert len(t1["misoriented_faces"]) == 0
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
    assert face_geometry(new_points, new_faces)[1].min() >= (
        face_geometry(points, faces)[1].min() - 1e-6
    )
