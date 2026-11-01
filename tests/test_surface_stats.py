import numpy as np
import pytest
from lib.surface_stats import summarize_surface_metrics, vertex_areas, weighted_stats


def strip(xs, y=(0.0, 1.0)):
    """Triangulated strip of quads between the given x and y positions."""
    xx, yy = np.meshgrid(xs, y, indexing="ij")
    points = np.c_[xx.ravel(), yy.ravel(), np.zeros(xx.size)]
    ny = len(y)
    faces = []
    for i in range(len(xs) - 1):
        for j in range(ny - 1):
            a, b = i * ny + j, (i + 1) * ny + j
            faces += [[a, b, b + 1], [a, b + 1, a + 1]]
    return points, np.array(faces)


def test_vertex_areas_unit_square():
    # two triangles of area 1/2; the two vertices on the diagonal are in both
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    np.testing.assert_allclose(
        vertex_areas(points, faces), [1 / 3, 1 / 6, 1 / 3, 1 / 6]
    )


def test_vertex_areas_sum_and_rigid_invariance():
    points, faces = strip(np.linspace(0, 2, 5) ** 2)  # rectangle 4 x 1
    areas = vertex_areas(points, faces)
    assert areas.sum() == pytest.approx(4.0)
    c, s = np.cos(0.7), np.sin(0.7)
    rotation = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    moved = points @ rotation.T + [5, -3, 2]
    np.testing.assert_allclose(vertex_areas(moved, faces), areas)


def test_vertex_areas_unused_vertex():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [9, 9, 9]], dtype=float)
    np.testing.assert_allclose(
        vertex_areas(points, np.array([[0, 1, 2]])), [1 / 6, 1 / 6, 1 / 6, 0]
    )


def test_weighted_stats_by_hand():
    # weights 1 and 3: mean (1*2 + 3*6)/4 = 5 (the vertex mean would be 4),
    # variance (1*(2-5)^2 + 3*(6-5)^2)/4 = 3
    mean, std, area, n = weighted_stats([2.0, 6.0], [1.0, 3.0])
    assert mean == pytest.approx(5.0)
    assert std == pytest.approx(np.sqrt(3.0))
    assert (area, n) == (4.0, 2)


def test_weighted_stats_nan_and_empty():
    mean, std, area, n = weighted_stats(
        [2.0, np.nan, 6.0, np.inf], [1.0, 5.0, 3.0, 2.0]
    )
    assert (mean, area, n) == (pytest.approx(5.0), 4.0, 2)
    mean, std, area, n = weighted_stats([np.nan], [1.0])
    assert np.isnan(mean) and np.isnan(std) and (area, n) == (0.0, 0)


def test_mean_is_independent_of_triangulation():
    # a linear metric, f = x, on the rectangle [0, 2] x [0, 1] has mean 1 over
    # the surface. With vertices 10 times denser on the left half, the vertex
    # mean is pulled to the left, but the area-weighted mean is exact.
    xs = np.r_[np.linspace(0, 1, 41), np.linspace(1, 2, 5)[1:]]
    points, faces = strip(xs, np.linspace(0, 1, 5))
    x = points[:, 0]
    df = summarize_surface_metrics({"x": x}, vertex_areas(points, faces), name="all")
    assert df["mean"].item() == pytest.approx(1.0, abs=1e-12)
    assert df["area"].item() == pytest.approx(2.0)
    assert x.mean() < 0.7


def test_subfields_from_labels():
    # three columns of quads, of widths 1, 1 and 2, with labels 1 (left
    # vertices), 2 (right vertices) and 0 (none) on the middle ones
    points, faces = strip([0.0, 1.0, 2.0, 4.0])
    areas = vertex_areas(points, faces)
    # vertices are ordered (x, y): x = 0, 0, 1, 1, 2, 2, 4, 4
    labels = np.array([1, 1, 0, 0, 2, 2, 2, 2])
    thickness = np.array([1.0, 1.0, 9.0, 9.0, 2.0, 2.0, 4.0, 4.0])
    names = {0: "???", 1: "Sub", 2: "CA1", 6: "DG"}

    df = summarize_surface_metrics({"thickness": thickness}, areas, labels, names)

    assert list(df.columns) == [
        "subfield",
        "metric",
        "mean",
        "std",
        "area",
        "nvertices",
    ]
    assert df["subfield"].tolist() == ["Sub", "CA1"]  # no unassigned (0), no DG
    sub, ca1 = df.iloc[0], df.iloc[1]
    assert (sub["mean"], sub["std"], sub["nvertices"]) == (1.0, 0.0, 2)
    assert sub["area"] == pytest.approx(0.5)  # 3 corners of triangles of area 1/2
    # CA1: the two vertices at x = 2 are 3 corners of the triangles of area 1/2
    # to their left and 3 of the triangles of area 1 to their right, so have
    # area (3 * 1/2 + 3 * 1) / 3 = 1.5; those at x = 4 have area 3 * 1 / 3 = 1.
    # So the mean is (1.5 * 2 + 1 * 4) / 2.5 = 2.8, not the vertex mean, 3
    assert ca1["area"] == pytest.approx(2.5)
    assert ca1["mean"] == pytest.approx(2.8)
    assert ca1["nvertices"] == 4


def test_unnamed_label_and_metric_order():
    points, faces = strip([0.0, 1.0])
    areas = vertex_areas(points, faces)
    df = summarize_surface_metrics(
        {"b": np.ones(4), "a": np.zeros(4)}, areas, np.array([7, 7, 7, 7]), {}
    )
    assert df[["subfield", "metric"]].values.tolist() == [["7", "b"], ["7", "a"]]


def test_mismatched_lengths():
    with pytest.raises(ValueError, match="3 values for a surface with 4 vertices"):
        summarize_surface_metrics({"thickness": np.ones(3)}, np.ones(4), name="all")
    with pytest.raises(ValueError, match="parcel labels"):
        summarize_surface_metrics({"thickness": np.ones(4)}, np.ones(4), np.ones(3))
