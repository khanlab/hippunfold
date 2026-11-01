import numpy as np
import pytest
from lib.gyrification import boundary_vertices, smoothing_gyrification
from scipy.spatial import Delaunay

SIZE = 12.0  # mm, side of the test sheet


def structured_sheet(n=97):
    """square sheet, each grid cell split into two triangles"""
    x, y = np.meshgrid(np.linspace(0, SIZE, n), np.linspace(0, SIZE, n), indexing="ij")
    idx = np.arange(n * n).reshape(n, n)
    a, b, c, d = (
        idx[:-1, :-1].ravel(),
        idx[1:, :-1].ravel(),
        idx[:-1, 1:].ravel(),
        idx[1:, 1:].ravel(),
    )
    faces = np.r_[np.c_[a, b, c], np.c_[c, b, d]]
    return np.c_[x.ravel(), y.ravel()], faces


def random_sheet(n_interior=7000, n_edge=96, seed=0):
    """the same square, Delaunay-triangulated from random points"""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, SIZE, n_edge, endpoint=False)
    z = np.zeros_like(t)
    edge = np.r_[
        np.c_[t, z], np.c_[z + SIZE, t], np.c_[SIZE - t, z + SIZE], np.c_[z, SIZE - t]
    ]
    xy = np.r_[edge, rng.uniform(0.1, SIZE - 0.1, (n_interior, 2))]
    return xy, Delaunay(xy).simplices


def corrugate(xy, amplitude=0.4, wavelength=1.5):
    """sinusoidal folds along x, like digitations"""
    z = amplitude * np.sin(2 * np.pi * xy[:, 0] / wavelength)
    return np.c_[xy, z]


def interior(xy, margin=3.0):
    return np.all((xy > margin) & (xy < SIZE - margin), axis=1)


def expected_ratio(amplitude=0.4, wavelength=1.5):
    """arc length of one period of the sine over its wavelength"""
    x = np.linspace(0, wavelength, 20001)
    slope = amplitude * 2 * np.pi / wavelength * np.cos(2 * np.pi * x / wavelength)
    return np.trapezoid(np.sqrt(1 + slope**2), x) / wavelength


def test_flat_sheet_is_one():
    xy, faces = structured_sheet(41)
    sgi = smoothing_gyrification(np.c_[xy, np.zeros(len(xy))], faces)
    np.testing.assert_allclose(sgi, 1.0, atol=1e-9)


def test_corrugated_sheet_matches_area_ratio():
    """folds well below the scale are flattened, so in the interior the index
    is the area the folds add"""
    xy, faces = structured_sheet()
    sgi = smoothing_gyrification(corrugate(xy), faces, scale=3.0)
    assert np.median(sgi[interior(xy)]) == pytest.approx(expected_ratio(), rel=0.03)


def test_folds_much_wider_than_scale_hardly_count():
    xy, faces = structured_sheet()
    narrow = smoothing_gyrification(corrugate(xy, 0.4, 1.5), faces, scale=1.5)
    wide = smoothing_gyrification(corrugate(xy, 0.4, 12.0), faces, scale=1.5)
    assert np.median(wide[interior(xy)]) < 1.01 < np.median(narrow[interior(xy)])


def test_independent_of_triangulation():
    """the same corrugated sheet triangulated two ways gives the same map"""
    xy1, f1 = structured_sheet()
    xy2, f2 = random_sheet()

    def sheet(xy):  # fold amplitude varies across the sheet, so the map is not flat
        z = (0.15 + 0.03 * xy[:, 1]) * np.sin(2 * np.pi * xy[:, 0] / 2.0)
        return np.c_[xy, z]

    s1 = smoothing_gyrification(sheet(xy1), f1)
    s2 = smoothing_gyrification(sheet(xy2), f2)
    # compare at the random sheet's interior vertices, interpolating the grid
    from scipy.interpolate import LinearNDInterpolator

    keep = interior(xy2, margin=2.0)
    s1_at_2 = LinearNDInterpolator(xy1, s1)(xy2[keep])
    assert np.corrcoef(s1_at_2, s2[keep])[0, 1] > 0.99
    np.testing.assert_allclose(s2[keep], s1_at_2, rtol=0.02)


def test_hole_rims_are_pinned():
    """a free rim would retract when smoothed and inflate the index around
    the hole; with every boundary vertex pinned a flat holed sheet stays 1"""
    xy, faces = structured_sheet(61)
    centre = np.all(np.abs(xy - SIZE / 2) < 1.5, axis=1)
    faces = faces[~np.any(centre[faces], axis=1)]
    used = np.unique(faces)
    remap = np.full(len(xy), -1)
    remap[used] = np.arange(len(used))
    xy, faces = xy[used], remap[faces]

    rim = boundary_vertices(faces)
    assert (
        np.sum(np.all(np.abs(xy[rim] - SIZE / 2) < 2.0, axis=1)) > 0
    )  # the hole has a rim
    sgi = smoothing_gyrification(np.c_[xy, np.zeros(len(xy))], faces)
    np.testing.assert_allclose(sgi, 1.0, atol=1e-9)
