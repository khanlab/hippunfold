import numpy as np
import pytest
import pyvista as pv
from lib.surface import orient_faces_up_gradient, vertex_normals


def contour_sphere(affine):
    """isosurface of a radial ramp (increasing outward), contoured in voxel
    space and moved to world space as gen_isosurface does"""
    n = 24
    ijk = np.stack(np.meshgrid(*[np.arange(n)] * 3, indexing="ij"), -1)
    volume = np.linalg.norm(ijk + 0.5 - n / 2, axis=-1) / (n / 2)
    grid = pv.ImageData(dimensions=np.array(volume.shape) + 1)
    grid.cell_data["values"] = volume.flatten(order="F")
    surface = grid.cells_to_points("values").contour([0.6])
    surface = surface.transform(affine, inplace=False)
    return surface.points, surface.faces.reshape((-1, 4))[:, 1:4], volume


@pytest.mark.parametrize("flip_axis", [None, 0, 2])
def test_normals_point_up_gradient_for_either_handedness(flip_axis):
    affine = np.diag([0.3, 0.3, 0.3, 1.0])
    if flip_axis is not None:
        affine[flip_axis, flip_axis] *= -1  # negative determinant (e.g. LAS)
    points, faces, volume = contour_sphere(affine)
    for start in (faces, faces[:, ::-1]):  # whatever the contour's winding
        oriented, _ = orient_faces_up_gradient(points, start, volume, affine)
        outward = points - points.mean(0)
        assert np.all(np.sum(vertex_normals(points, oriented) * outward, 1) > 0)


def test_flip_reported_and_consistent_winding_kept():
    affine = np.diag([0.3, 0.3, 0.3, 1.0])
    points, faces, volume = contour_sphere(affine)
    oriented, _ = orient_faces_up_gradient(points, faces, volume, affine)
    same, flipped = orient_faces_up_gradient(points, oriented, volume, affine)
    again, flipped_again = orient_faces_up_gradient(
        points, oriented[:, ::-1], volume, affine
    )
    assert not flipped and flipped_again
    assert np.array_equal(same, oriented) and np.array_equal(again, oriented)
    # whole-mesh decision: every face reversed or none
    assert np.array_equal(oriented, faces) or np.array_equal(oriented, faces[:, ::-1])
