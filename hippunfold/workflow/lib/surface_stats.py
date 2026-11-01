"""
Summary statistics of surface metrics over subfields, weighted by surface area.

The vertices of a surface are not evenly spread over it: e.g. the vertices of
the den-8k surfaces are evenly spaced in the unfolded space, so in native space
they are sparser where the surface is more gyrified. A plain mean over the
vertices of a subfield therefore weights parts of it by how densely they are
sampled. Weighting each vertex by its area (a third of the area of each
triangle it is part of) instead estimates the mean over the surface, whatever
the triangulation: over a whole surface, sum(area * value) is exactly the
integral of a metric that is linear over each triangle.
"""

import numpy as np
import pandas as pd

COLUMNS = ["subfield", "metric", "mean", "std", "area", "nvertices"]


def vertex_areas(points, faces):
    """Area of each vertex: a third of the area of each triangle it is part of,
    as in wb_command -surface-vertex-areas (vertices in no triangle get 0)."""
    points = np.asarray(points, dtype=np.float64)
    faces = np.asarray(faces)
    tri = points[faces]
    face_area = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1
    )
    return np.bincount(
        faces.ravel(), weights=np.repeat(face_area / 3, 3), minlength=len(points)
    )


def weighted_stats(values, weights):
    """Weighted mean and standard deviation of the finite values, with the
    total weight and number of those values. The standard deviation is
    normalised by the total weight (no bias correction, as the weights are
    areas, not counts). The mean and standard deviation are NaN if the total
    weight is 0."""
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    finite = np.isfinite(values)
    values, weights = values[finite], weights[finite]
    total = weights.sum()
    if total > 0:
        mean = np.sum(weights * values) / total
        std = np.sqrt(np.sum(weights * (values - mean) ** 2) / total)
    else:
        mean = std = np.nan
    return mean, std, total, int(finite.sum())


def summarize_surface_metrics(metrics, areas, parcels=None, names=None, name=None):
    """
    Area-weighted mean and standard deviation of each metric over each parcel
    (e.g. subfield) of a surface.

    Args:
        metrics (dict): metric name -> value at each vertex.
        areas (array): area of each vertex (see vertex_areas).
        parcels (array, optional): integer label of each vertex. If None, the
            whole surface is summarized as a single parcel called name.
        names (dict, optional): label -> parcel name (e.g. from the label
            table of a label.gii); labels missing from it are named by their
            number.
        name (str, optional): name of the parcel if parcels is None.

    Only the labels found on the surface are summarized, except label 0
    (unassigned), and vertices with a non-finite value (NaN or inf) of a
    metric are left out for that metric. So area and nvertices are of the
    vertices each mean was computed over.

    Returns:
        pandas.DataFrame: one row per metric and parcel, with columns
        subfield, metric, mean, std, area and nvertices.
    """
    areas = np.asarray(areas, dtype=np.float64)
    if parcels is None:
        parcel_masks = [(name, np.ones(len(areas), dtype=bool))]
    else:
        parcels = np.asarray(parcels)
        if parcels.shape != areas.shape:
            raise ValueError(
                f"{len(parcels)} parcel labels for a surface with {len(areas)} vertices"
            )
        names = names or {}
        parcel_masks = [
            (names.get(label, str(label)), parcels == label)
            for label in np.unique(parcels)
            if label != 0
        ]

    rows = []
    for metric, values in metrics.items():
        values = np.asarray(values, dtype=np.float64)
        if values.shape != areas.shape:
            raise ValueError(
                f"{metric} has {len(values)} values for a surface with {len(areas)} vertices"
            )
        for parcel, mask in parcel_masks:
            mean, std, area, nvertices = weighted_stats(values[mask], areas[mask])
            rows.append([parcel, metric, mean, std, area, nvertices])

    return pd.DataFrame(rows, columns=COLUMNS)
