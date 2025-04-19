import re
import numpy as np
import pandas as pd
import pyvista as pv
from lib.surface import compute_edge_lengths, read_surface_from_gifti


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
