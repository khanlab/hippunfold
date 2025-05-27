import json
import re
import os
import numpy as np
from lib.surface import read_surface_from_gifti, write_label_gii

# --- Snakemake Inputs ---
vertex_hits_json = snakemake.input.hits
surface_paths = snakemake.input.surfs  # all 3 per group
label_gii_paths = snakemake.output.labelgii  # now 1 per group (not per surface)

# --- Load vertex hits ---
with open(vertex_hits_json, "r") as f:
    vertex_hits = json.load(f)  # keys are e.g. sub-XX_ses-YY_hemi-L_..._midthickness.surf.gii


def extract_group_key(filename):
    fname = os.path.basename(filename)
    hemi_match = re.search(r"hemi-[LR]", fname)
    label_match = re.search(r"label-[^_]+", fname)

    if not hemi_match or not label_match:
        raise ValueError(f"Could not extract group key from filename: {filename}")

    return f"{hemi_match.group()}_{label_match.group()}"


grouped_hits = {}
for surface_name, hits in vertex_hits.items():
    group_key = extract_group_key(surface_name)
    if group_key not in grouped_hits:
        grouped_hits[group_key] = {}
    grouped_hits[group_key].update(hits)  # merge across midthickness/inner/outer

# --- Map group key → output label path ---
output_map = {
    extract_group_key(os.path.basename(p)): p for p in label_gii_paths
}

# --- Write one label GIFTI per group ---
for group_key, vertex_dict in grouped_hits.items():
    out_path = output_map[group_key]

    # Pick any surface file from the group (e.g., first midthickness)
    surface_path = next(
        s for s in surface_paths if extract_group_key(s) == group_key and "midthickness" in s
    )
    mesh, metadata = read_surface_from_gifti(surface_path)
    n_verts = mesh.n_points
    label_scalars = np.zeros(n_verts, dtype=np.int32)

    # Assign unique label IDs
    label_to_id = {}
    label_table = {}
    next_id = 1

    for vert_str, label in vertex_dict.items():
        if label not in label_to_id:
            label_to_id[label] = next_id
            label_table[next_id] = {
                "Label": label,
                "Red": 255,
                "Green": 0,
                "Blue": 0,
                "Alpha": 1.0,
                "Key": next_id,
            }
            next_id += 1
        label_scalars[int(vert_str)] = label_to_id[label]

    write_label_gii(label_scalars, out_path, label_dict=label_table, metadata=metadata)

