import json
import os
import numpy as np
from lib.surface import read_surface_from_gifti, write_label_gii

# --- Snakemake Inputs ---
vertex_hits_json = snakemake.input.hits
surface_paths = snakemake.input.surfs  # just midthickness
label_gii_paths = snakemake.output.labelgii  # same length/order

# --- Load vertex hits ---
with open(vertex_hits_json, "r") as f:
    vertex_hits = json.load(f)

# --- Loop over each (surface, label output) pair ---
for surf_path, out_path in zip(surface_paths, label_gii_paths):
    basename = os.path.basename(surf_path)
    mesh, metadata = read_surface_from_gifti(surf_path)
    n_verts = mesh.n_points
    label_scalars = np.zeros(n_verts, dtype=np.int32)

    # Always start with background label
    label_table = {
        "Background": {"key": 0, "red": 1.0, "green": 1.0, "blue": 1.0, "alpha": 0.0}
    }

    if basename in vertex_hits:
        for integer_id, (vert_str, label) in enumerate(vertex_hits[basename].items(), start=1):
            label_table[f'sub-{snakemake.wildcards.subject}_{label}'] = {
                "key": integer_id,
                "red": 1.0,
                "green": 0.0,
                "blue": 0.0,
                "alpha": 1.0,
            }
            label_scalars[int(vert_str)] = integer_id
    else:
        print(f"Note: {basename} has no vertex hits — writing background-only label GIFTI.")

    write_label_gii(label_scalars, out_path, label_dict=label_table, metadata=metadata)

