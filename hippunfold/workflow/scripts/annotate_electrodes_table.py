import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import os
import nibabel as nib
from lib.surface import read_surface_from_gifti

# --- Inputs from Snakemake ---
electrode_tsv = snakemake.input.tsv
xfm_path = snakemake.input.xfm
surface_paths = snakemake.input.surfs
label_left_path = snakemake.input.label_left
label_right_path = snakemake.input.label_right
output_path = snakemake.output.annotated
max_dist_mm = 5.0

# --- Helper functions ---
def apply_transform(pts, affine):
    pts_h = np.c_[pts, np.ones(len(pts))]
    return (affine @ pts_h.T).T[:, :3]

def load_label_scalars_and_names(label_gii_path):
    gii = nib.load(label_gii_path)
    scalars = gii.darrays[0].data.astype(int)
    labeltable = gii.labeltable
    label_map = {lbl.key: lbl.label for lbl in labeltable.labels}
    label_names = [label_map.get(idx, "n/a") for idx in scalars]
    return label_names

def determine_side(path):
    name = os.path.basename(path).lower()
    if "hemi-l" in name or "left" in name or "lh" in name:
        return "left"
    elif "hemi-r" in name or "right" in name or "rh" in name:
        return "right"
    return None

# --- Load and transform electrode coordinates ---
df = pd.read_csv(electrode_tsv, sep="\t")
coords = df[["x", "y", "z"]].values
xfm = np.loadtxt(xfm_path)
coords = apply_transform(coords, xfm)

# --- Load label GIFTIs ---
labels_by_side = {
    "left": load_label_scalars_and_names(label_left_path),
    "right": load_label_scalars_and_names(label_right_path),
}

# --- Load all surfaces and collect vertex + label data ---
all_points = []
all_labels = []
all_sources = []
all_local_indices = []  # track original vertex index
all_surface_keys = []   # track surface (filename) of each point

for surf_path in surface_paths:
    mesh, _ = read_surface_from_gifti(surf_path)
    verts = mesh.points
    source = os.path.basename(surf_path)
    side = determine_side(source)
    is_dentate = "dentate" in source

    if is_dentate:
        labels = ["Dentate"] * verts.shape[0]
    elif side in labels_by_side:
        labels = labels_by_side[side]
    else:
        raise ValueError(f"Cannot determine label side for surface: {source}")

    if len(labels) != verts.shape[0]:
        raise ValueError(f"Label count mismatch for {source}: {len(labels)} vs {verts.shape[0]}")

    all_points.append(verts)
    all_labels.extend(labels)
    all_sources.extend([source] * verts.shape[0])
    all_local_indices.extend(list(range(verts.shape[0])))
    all_surface_keys.extend([source] * verts.shape[0])

# --- Build KDTree and query nearest vertex ---
all_points = np.vstack(all_points)
tree = cKDTree(all_points)
dists, global_indices = tree.query(coords)

# --- Recover per-surface local indices and metadata ---
matched_labels = []
matched_sources = []
matched_local_indices = []

for d, gi in zip(dists, global_indices):
    if d <= max_dist_mm:
        matched_labels.append(all_labels[gi])
        matched_sources.append(all_sources[gi])
        matched_local_indices.append(all_local_indices[gi])
    else:
        matched_labels.append("n/a")
        matched_sources.append("n/a")
        matched_local_indices.append(-1)

# --- Add annotations to dataframe ---
df["closest_label"] = matched_labels
df["distance_mm"] = dists
df["surface_source"] = matched_sources
df["vertex_index"] = matched_local_indices

# --- Write output TSV ---
df.to_csv(output_path, sep="\t", index=False)


