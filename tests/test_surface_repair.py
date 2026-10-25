import numpy as np
from lib.qc_surface import compute_qc, compute_topology
from lib.surface_repair import (
    fill_holes,
    remove_bad_faces,
    repair_mesh,
    split_nonmanifold_vertices,
)
from test_qc_surface import disk, grid


def holed_disk(half_width=0.2):
    points, faces, u, v = disk()
    centre = (np.abs(u - 0.5) < half_width) & (np.abs(v - 0.5) < half_width)
    keep = ~np.all(np.isin(faces, np.flatnonzero(centre)), axis=1)
    return points, faces[keep]


def topo(points, faces):
    t = compute_topology(points, faces)
    return t["boundary_loops"], t["genus"], len(t["misoriented_faces"])


def test_fill_holes():
    points, faces = holed_disk()
    assert topo(points, faces) == (2, 0, 0)

    filled_points, filled, radii = fill_holes(points, faces, max_radius=5)
    assert len(radii) == 1
    assert len(filled_points) == len(points) + 1  # the fan centre
    assert topo(filled_points, filled) == (1, 0, 0)  # winding consistent
    assert np.array_equal(filled[: len(faces)], faces)  # existing faces kept

    # too large to fill, and the outer boundary is never filled
    _, unfilled, radii = fill_holes(points, faces, max_radius=0.1)
    assert radii == [] and len(unfilled) == len(faces)
    p, f, _ = disk()[:3]
    assert fill_holes(p, f, max_radius=100)[2] == []


def test_repair_misoriented_and_bad_faces():
    points, faces, _, _ = disk()
    bad = faces.copy()
    bad[[3, 20]] = bad[[3, 20], ::-1]
    bad = np.vstack([bad, bad[5], [[0, 0, 1]]])  # duplicate, repeated vertex
    repaired_points, repaired, counts = repair_mesh(points, bad)
    assert counts["misoriented_faces_flipped"] == 2
    assert counts["duplicate_faces"] == 1
    assert counts["repeated_vertex_faces"] == 1
    assert topo(repaired_points, repaired) == (1, 0, 0)
    assert len(repaired) == len(faces)
    # the majority winding is kept
    assert np.array_equal(repaired, faces)


def test_remove_bad_faces_keeps_first():
    faces = np.array([[0, 1, 2], [1, 2, 3], [2, 0, 1], [3, 3, 1]])
    kept, counts = remove_bad_faces(None, faces)
    assert kept.tolist() == [[0, 1, 2], [1, 2, 3]]
    assert counts == {"repeated_vertex_faces": 1, "duplicate_faces": 1}


def test_repair_nonmanifold_edge():
    points, faces, _, _ = disk(4)
    e = faces[5, :2]
    points = np.r_[points, [[5, 5, 0.01]]]
    fin = np.r_[faces, [[e[1], e[0], len(points) - 1]]]  # small fin
    repaired_points, repaired, counts = repair_mesh(points, fin)
    assert counts["nonmanifold_edge_faces"] == 1
    assert counts["unreferenced_vertices"] == 1
    t = compute_topology(repaired_points, repaired)
    assert len(t["nonmanifold_edges"]) == 0 and t["chi"] == 1


def test_split_bowtie_and_small_component():
    # two disks touching at one vertex: split, then keep the larger one
    u, v, faces = grid(4, 4)
    a = np.c_[u, v, np.zeros_like(u)]
    u2, v2, faces2 = grid(2, 2)
    b = np.c_[u2 - 0.5 + 0.0, v2 - 0.5, np.zeros_like(u2)] * 0.5 - [0.25, 0.25, 0]
    points = np.r_[a, b]
    faces2 = faces2 + len(a)
    corner = np.flatnonzero((u2 == 1) & (v2 == 1))[0] + len(a)
    faces2[faces2 == corner] = 0  # share vertex 0 with the first disk
    bowtie = np.r_[faces, faces2]
    assert len(compute_topology(points, bowtie)["nonmanifold_vertices"]) == 1

    _, split, n = split_nonmanifold_vertices(points, bowtie)
    assert (
        n == 1
        and len(
            compute_topology(np.r_[points, points[[0]]], split)["nonmanifold_vertices"]
        )
        == 0
    )

    repaired_points, repaired, counts = repair_mesh(points, bowtie)
    assert counts["nonmanifold_vertices_split"] == 1
    assert counts["small_component_faces"] == len(faces2)
    assert len(repaired) == len(faces)
    assert compute_topology(repaired_points, repaired)["components"] == 1


def test_qc_without_coords():
    points, faces = holed_disk()
    result = compute_qc({"midthickness": (points, faces)}, None, None, {}, {})
    names = {c["name"] for c in result.checks}
    assert "holes" in names and not any(n.startswith("unfold") for n in names)
    hole = next(d for d in result.data["defects"] if d["type"] == "hole")
    assert hole["ap"] is None
