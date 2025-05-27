import pandas as pd
import json
import os

# --- Snakemake Inputs ---
annotated_tsv = snakemake.input.annotated
output_json = snakemake.output.hits

# --- Load annotated electrode data ---
df = pd.read_csv(annotated_tsv, sep="\t")

# Filter to matched rows only
df = df[(df["surface_source"] != "n/a") & (df["vertex_index"] >= 0)]


# --- Build surface-to-vertex-electrode mapping ---
surface_hits = {}

for i, row in df.iterrows():
    surface = row["surface_source"]
    if surface == "n/a" or pd.isna(surface):
        continue

    vert_idx = int(row["vertex_index"]) if "vertex_index" in row else None
    if pd.isna(vert_idx):
        continue

    label = str(row["label"]) if "label" in row else f"electrode_{i}"

    if surface not in surface_hits:
        surface_hits[surface] = {}
    surface_hits[surface][str(vert_idx)] = label

# --- Write to JSON ---
with open(output_json, "w") as f:
    json.dump(surface_hits, f, indent=2)

