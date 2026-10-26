"""Smoothing-based gyrification index (sGI) of an open surface.

The surface is compared with a heat-smoothed copy of itself: every fold
narrower than about the smoothing scale is flattened in the copy, so the local
ratio of the surface's area to the copy's area measures how much area those
folds add. The heat flow uses the cotangent Laplacian with a lumped mass
matrix (a finite-element discretisation), so the index does not depend on the
triangulation, and all smoothing stays on the surface, so it never mixes parts
of the sheet that are close in space but far apart along it.

Boundary vertices (the outer loop and the rims of any holes) are pinned while
smoothing the copy: a free rim would retract and inflate the ratio next to it.
A flat sheet gives exactly 1.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import factorized


def face_areas(points, faces):
    a, b, c = (points[faces[:, k]] for k in range(3))
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def vertex_areas(points, faces, areas=None):
    """lumped mass: a third of each incident face's area (or of `areas`)"""
    areas = face_areas(points, faces) if areas is None else areas
    out = np.zeros(len(points))
    np.add.at(out, faces.ravel(), np.repeat(areas / 3, 3))
    return out


def cotangent_laplacian(points, faces):
    """cotangent stiffness matrix, with negative edge weights (opposite obtuse
    angles) clamped to zero so the heat flow stays monotone"""
    n = len(points)
    rows, cols, vals = [], [], []
    for k in range(3):
        i, j, o = faces[:, (k + 1) % 3], faces[:, (k + 2) % 3], faces[:, k]
        u, w = points[i] - points[o], points[j] - points[o]
        cross = np.linalg.norm(np.cross(u, w), axis=1)
        cot = np.zeros(len(faces))
        safe = cross > 1e-12
        cot[safe] = np.einsum("ij,ij->i", u[safe], w[safe]) / cross[safe]
        rows += [i, j]
        cols += [j, i]
        vals += [cot / 2, cot / 2]
    weights = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsr()
    weights.data = np.maximum(weights.data, 0)
    return sp.diags(weights.sum(axis=1).A1) - weights


def boundary_vertices(faces):
    """vertices on any boundary edge (outer loop and hole rims)"""
    edges = np.sort(
        np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1
    )
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return np.unique(unique[counts == 1])


class HeatSmoother:
    """Gaussian-like smoothing on the surface: `steps` implicit heat steps of
    total time sigma^2 / 2 (so sigma is comparable to a Gaussian's standard
    deviation, in mm). Vertices in `fixed` keep their values."""

    def __init__(self, points, faces, sigma, steps=4, fixed=None):
        self.mass = vertex_areas(points, faces)
        self.steps = steps
        system = sp.diags(self.mass) + (sigma**2 / 2 / steps) * cotangent_laplacian(
            points, faces
        )
        self.keep = np.zeros(len(points), bool)
        if fixed is not None:
            self.keep[fixed] = True
            free = sp.diags((~self.keep).astype(float))
            system = free @ system + sp.diags(self.keep.astype(float))
        self.solve = factorized(system.tocsc())

    def __call__(self, values):
        values = np.asarray(values, float)
        columns = values.reshape(len(values), -1)
        out = []
        for c in columns.T:
            for _ in range(self.steps):
                rhs = self.mass * c
                rhs[self.keep] = c[self.keep]
                c = self.solve(rhs)
            out.append(c)
        return np.stack(out, 1).reshape(values.shape)


def smoothing_gyrification(points, faces, scale=1.5, sigma=1.5, steps=4):
    """Smoothing-based gyrification index at each vertex.

    sGI = 1 / S_sigma[ A(X_scale) / A(X) ], where X_scale is the surface
    heat-smoothed at `scale` (mm) with its boundary pinned, A(.) is the lumped
    vertex area of the same triangulation, and S_sigma is heat smoothing on
    the original surface (`sigma`, mm), which pools the area ratio over a
    neighbourhood before inverting it.
    """
    points = np.asarray(points, float)
    faces = np.asarray(faces, int)
    smoothed = HeatSmoother(
        points, faces, scale, steps, fixed=boundary_vertices(faces)
    )(points)
    density = vertex_areas(points, faces, face_areas(smoothed, faces)) / vertex_areas(
        points, faces
    )
    return 1.0 / HeatSmoother(points, faces, sigma, steps)(density)
