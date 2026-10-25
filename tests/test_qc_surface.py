import numpy as np
from lib.qc import ERROR
from lib.qc_surface import (
    compute_qc,
    compute_topology,
    face_clusters,
    plot_defects,
    plot_distortion,
    plot_histograms,
    self_intersecting_faces,
    unfold_distortion,
)


def grid(n=8, m=8, wrap_u=False, wrap_v=False):
    """Triangulated n x m grid in [0, 1]^2, optionally wrapped (consistent winding)."""
    nu, nv = (n if wrap_u else n + 1), (m if wrap_v else m + 1)
    u, v = np.meshgrid(np.arange(nu), np.arange(nv), indexing="ij")
    idx = lambda i, j: (i % nu) * nv + (j % nv)
    faces = []
    for i in range(n):
        for j in range(m):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            faces += [[a, b, c], [a, c, d]]
    return u.ravel() / n, v.ravel() / m, np.array(faces)


def disk(n=8):
    u, v, faces = grid(n, n)
    points = np.c_[u * 10, v * 10, 0.5 * np.sin(3 * u)]
    return points, faces, u, v


def torus(n=12, m=8):
    u, v, faces = grid(n, m, wrap_u=True, wrap_v=True)
    a, b = 2 * np.pi * u, 2 * np.pi * v
    points = np.c_[(3 + np.cos(b)) * np.cos(a), (3 + np.cos(b)) * np.sin(a), np.sin(b)]
    return points, faces


def checks_by_name(result):
    return {c["name"]: c for c in result.checks}


def qc(points, faces, ap, pd, **thresholds):
    surfaces = {"midthickness": (points, faces)}
    return compute_qc(surfaces, ap, pd, {"hemi": "L", "label": "hipp"}, thresholds)


def test_disk_topology():
    points, faces, u, v = disk()
    t = compute_topology(points, faces)
    assert (t["chi"], t["components"], t["boundary_loops"], t["genus"]) == (1, 1, 1, 0)
    assert t["orientable"] and len(t["misoriented_faces"]) == 0
    assert len(t["nonmanifold_edges"]) == 0 and len(t["nonmanifold_vertices"]) == 0

    # the grid diagonals make ears at two corners: the only defects
    result = qc(points, faces, u, v)
    assert [c["name"] for c in result.failed()] == ["boundary_chords"]
    assert {d["type"] for d in result.data["defects"]} == {"boundary_chord"}
    assert result.metrics["boundary_chords"]["value"] == 2


def test_hole():
    points, faces, u, v = disk()
    centre = np.flatnonzero((np.abs(u - 0.5) < 0.2) & (np.abs(v - 0.5) < 0.2))
    keep = ~np.all(np.isin(faces, centre), axis=1)
    t = compute_topology(points, faces[keep])
    assert (t["boundary_loops"], t["genus"], t["chi"]) == (2, 0, 0)

    result = qc(points, faces[keep], u, v)
    checks = checks_by_name(result)
    assert not checks["holes"]["passed"]
    assert checks["genus"]["passed"]
    holes = [d for d in result.data["defects"] if d["type"] == "hole"]
    assert len(holes) == 1 and not holes[0]["at_boundary"]


def test_handle():
    # torus with one face removed: one boundary loop, genus 1
    points, faces = torus()
    t = compute_topology(points, faces[1:])
    assert (t["boundary_loops"], t["genus"], t["chi"]) == (1, 1, -1)
    # closed torus: no boundary is an error
    result = qc(points, faces, points[:, 0], points[:, 1])
    checks = checks_by_name(result)
    assert not checks["has_boundary"]["passed"]
    assert checks["has_boundary"]["severity"] == ERROR
    assert not checks["genus"]["passed"]


def test_misoriented_face():
    points, faces, u, v = disk()
    faces = faces.copy()
    faces[10] = faces[10, ::-1]
    t = compute_topology(points, faces)
    assert t["orientable"]
    assert list(t["misoriented_faces"]) == [10]

    result = qc(points, faces, u, v)
    check = checks_by_name(result)["orientable"]
    assert not check["passed"] and check["severity"] != ERROR
    types = {d["type"] for d in result.data["defects"]}
    assert {"misoriented", "unfold_flipped"} <= types
    flipped = next(d for d in result.data["defects"] if d["type"] == "unfold_flipped")
    assert "misoriented" in flipped["detail"]


def test_moebius_not_orientable():
    n = 12
    u, v, faces = grid(n, 1, wrap_u=False)
    nv = 2
    # glue the last column to the first with a half twist
    faces = np.where(faces >= n * nv, (1 - (faces - n * nv)), faces)
    a = 2 * np.pi * u
    w = v - 0.5
    points = np.c_[
        (1 + w * np.cos(a / 2)) * np.cos(a),
        (1 + w * np.cos(a / 2)) * np.sin(a),
        w * np.sin(a / 2),
    ]
    t = compute_topology(points, faces)
    assert not t["orientable"]
    assert t["boundary_loops"] == 1


def test_nonmanifold():
    points, faces, _, _ = disk(4)
    # a fin: extra face on an interior edge
    e = faces[5, :2]
    fin = np.r_[points, [[5, 5, 5]]]
    t = compute_topology(fin, np.r_[faces, [[e[0], e[1], len(points)]]])
    assert len(t["nonmanifold_edges"]) == 1

    # bowtie: two triangles meeting at one vertex
    bow = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
    t = compute_topology(bow, np.array([[0, 1, 2], [0, 3, 4]]))
    assert list(t["nonmanifold_vertices"]) == [0]


def test_self_intersection():
    points, faces, _, _ = disk(6)
    points = points.copy()
    points[:, 2] = 0
    n = len(points)
    # a triangle piercing the sheet near its centre
    spike = np.array([[4.6, 4.6, -1], [5.4, 4.6, 1], [5.0, 5.4, 1]])
    hits = self_intersecting_faces(
        np.r_[points, spike], np.r_[faces, [[n, n + 1, n + 2]]]
    )
    assert len(faces) in hits and len(hits) >= 2
    assert len(self_intersecting_faces(points, faces)) == 0


def test_unfold_flip_and_collapse():
    points, faces, u, v = disk()
    ap, pd = v.copy(), u.copy()
    # fold the mapping over at one interior vertex
    k = np.flatnonzero((np.abs(u - 0.5) < 1e-9) & (np.abs(v - 0.5) < 1e-9))[0]
    pd[k] += 0.3
    # collapse a terminal band: first column at PD = 0 and next column too
    pd[u <= 1 / 8 + 1e-9] = 0
    _, collapsed, flipped, log_ratio = unfold_distortion(points, faces, ap, pd)
    assert len(flipped) > 0
    assert len(collapsed) == 16 and np.isnan(log_ratio[collapsed]).all()

    result = qc(points, faces, ap, pd)
    checks = checks_by_name(result)
    assert not checks["unfold_flipped"]["passed"]
    assert not checks["unfold_collapsed"]["passed"]
    band = [d for d in result.data["defects"] if d["type"] == "unfold_collapsed"]
    assert band[0]["detail"] == "in PD src band"


def test_face_clusters():
    _, _, faces = grid(6, 6)
    clusters = face_clusters(faces, [0, 1, 70, 71], 49)
    assert [sorted(c) for c in clusters] == [[0, 1], [70, 71]]


def test_surface_figures():
    points, faces, u, v = disk()
    faces = faces.copy()
    faces[10] = faces[10, ::-1]
    result = qc(points, faces, u, v)
    _, _, flipped, log_ratio = unfold_distortion(points, faces, u, v)
    assert plot_defects(result, points, faces) is not None
    assert len(plot_distortion(points, faces, v, u, flipped, log_ratio)) == 2
    assert plot_histograms({"midthickness": np.full(10, 30.0)}, log_ratio, 5, 2)
