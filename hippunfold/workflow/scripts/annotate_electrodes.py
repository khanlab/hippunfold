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

# --- Load and transform electrode coordinates ---
df = pd.read_csv(electrode_tsv, sep="\t")
coords = df[["x", "y", "z"]].values
xfm = np.loadtxt(xfm_path)

def apply_transform(pts, affine):
    pts_h = np.c_[pts, np.ones(len(pts))]
    return (affine @ pts_h.T).T[:, :3]

coords = apply_transform(coords, xfm)

# --- Load label GIFTIs using your helper conventions ---
def load_label_scalars_and_names(label_gii_path):
    gii = nib.load(label_gii_path)
    scalars = gii.darrays[0].data.astype(int)
    labeltable = gii.labeltable
    label_map = {lbl.key: lbl.label for lbl in labeltable.labels}
    label_names = [label_map.get(idx, "n/a") for idx in scalars]
    return label_names

labels_by_side = {
    "left": load_label_scalars_and_names(label_left_path),
    "right": load_label_scalars_and_names(label_right_path),
}

def determine_side(path):
    name = os.path.basename(path).lower()
    if "hemi-l" in name or "left" in name or "lh" in name:
        return "left"
    elif "hemi-r" in name or "right" in name or "rh" in name:
        return "right"
    return None

# --- Load surfaces, build KDTree ---
all_points = []
all_labels = []
all_sources = []

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

all_points = np.vstack(all_points)
tree = cKDTree(all_points)

# --- Query nearest surface vertex for each electrode ---
dists, indices = tree.query(coords)

closest_labels = [
    all_labels[i] if d <= max_dist_mm else "n/a" for i, d in zip(indices, dists)
]
closest_sources = [
    all_sources[i] if d <= max_dist_mm else "n/a" for i, d in zip(indices, dists)
]

# --- Write output ---
df["closest_label"] = closest_labels
df["distance_mm"] = dists
df["surface_source"] = closest_sources
df.to_csv(output_path, sep="\t", index=False)

