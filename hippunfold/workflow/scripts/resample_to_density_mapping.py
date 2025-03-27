import re

import nibabel as nib
import numpy as np
import pandas as pd
import pyvista as pv
from nibabel.gifti.gifti import intent_codes


def read_surface_from_gifti(surf_gii):
    """Load a surface mesh from a GIFTI file."""
    surf = nib.load(surf_gii)
    vertices = surf.agg_data("NIFTI_INTENT_POINTSET")
    faces = surf.agg_data("NIFTI_INTENT_TRIANGLE")
    faces_pv = np.hstack([np.full((faces.shape[0], 1), 3), faces])  # PyVista format

    # Find the first darray that represents vertices (NIFTI_INTENT_POINTSET)
    vertices_darray = next(
        (
            darray
            for darray in surf.darrays
            if darray.intent == intent_codes["NIFTI_INTENT_POINTSET"]
        ),
        None,
    )

    # Extract metadata as a dictionary (return empty dict if no metadata)
    metadata = dict(vertices_darray.meta) if vertices_darray else {}

    return pv.PolyData(vertices, faces_pv), metadata


def compute_edge_lengths(surface):

    # Extract edges
    edges = surface.extract_all_edges()

    # Extract individual edge segments
    edge_lengths = []
    lines = edges.lines.reshape(-1, 3)  # Each row: [2, point1, point2]

    for line in lines:
        _, p1, p2 = line  # First entry is always "2" (pairs of points)
        length = np.linalg.norm(edges.points[p1] - edges.points[p2])
        edge_lengths.append(length)

    edge_lengths = np.array(edge_lengths)

    return edge_lengths


density_list = []


rows = []

for surf_gii in snakemake.input.surf_giis:
    # get the resample value from the filename
    resample = re.search(r"_resample-(\d+)_", surf_gii).group(1)
    hemi = re.search(r"_hemi-([LR])_", surf_gii).group(1)
    label = re.search(r"_label-(hipp|dentate)_", surf_gii).group(1)

    surface, metadata = read_surface_from_gifti(surf_gii)

    n_vertices = surface.points.shape[0]
    if n_vertices > 850:
        density = "{n}k".format(n=int(np.round(n_vertices / 1000)))
    else:
        density = str(n_vertices)

    # calculate vertex spacing (edge lengths)
    edge_lengths = compute_edge_lengths(surface)

    # Calculate statistics
    mean_length = np.mean(edge_lengths)
    min_length = np.min(edge_lengths)
    max_length = np.max(edge_lengths)
    median_length = np.median(edge_lengths)

    row = {
        "surf_gii": surf_gii,
        "hemi": hemi,
        "label": label,
        "resample": resample,
        "n_vertices": n_vertices,
        "min_spacing": min_length,
        "max_spacing": max_length,
        "mean_spacing": mean_length,
        "median_spacing": median_length,
        "density": density,
    }

    rows.append(row)


pd.DataFrame(rows).astype({"density": str, "resample": str}).to_csv(
    snakemake.output.mapping_csv, index=False
)
