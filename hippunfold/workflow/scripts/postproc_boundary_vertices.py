import pyvista as pv
import numpy as np
import nibabel as nib
import nibabel.gifti as gifti
from collections import Counter
from scipy.signal import argrelextrema
from lib.utils import setup_logger
from lib.surface import write_label_gii, read_surface_from_gifti

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


nmin = snakemake.params.min_terminal_vertices
max_iterations = snakemake.params.max_iterations  # Prevent infinite loops
shifting_epsilon = snakemake.params.shifting_epsilon

logger.info("Loading surface from GIFTI...")


edges = nib.load(snakemake.input.edges).agg_data()
ap_src = nib.load(snakemake.input.ap_src).agg_data()
ap_sink = nib.load(snakemake.input.ap_sink).agg_data()
pd_src = nib.load(snakemake.input.pd_src).agg_data()
pd_sink = nib.load(snakemake.input.pd_sink).agg_data()

# get structure metadata from edges
metadata = dict(nib.load(snakemake.input.edges).meta)


logger.info("Assigning labels apsrc, apsink, pdsrc, pdsink")

distances = np.vstack(
    (ap_src[edges == 1], ap_sink[edges == 1], pd_src[edges == 1], pd_sink[edges == 1])
).T
shifting_factors = np.zeros((4))
num_labels = 4  # starts at 0

for _ in range(max_iterations):
    # shift distances
    shifted_distances = distances - shifting_factors

    # Assign labels based on min scaled distance
    labels = np.argmin(shifted_distances, axis=1)

    # Count occurrences per label
    unique, counts = np.unique(labels, return_counts=True)
    label_counts = dict(zip(unique, counts))
    # Ensure all labels in range(num_labels) are present, setting missing ones to 0
    label_counts = {k: label_counts.get(k, 0) for k in range(num_labels)}

    logger.info(f"Iteration {_}, have label_counts={label_counts}")

    # Check if all labels meet nmin
    if all(count >= nmin for count in label_counts.values()):
        break  # Stop if all labels are sufficiently represented

    # Update scaling factors for underrepresented labels
    for k in range(num_labels):
        if label_counts[k] < nmin:
            shifting_factors[
                k
            ] += shifting_epsilon  # Increase competitiveness of the label
            logger.info(
                f"Shifting distances of label {k} to bump up competitiveness, factor = {shifting_factors[k]}"
            )

# Ensure all labels are represented
for k in range(num_labels):
    if label_counts[k] < nmin:
        raise ValueError(
            f"Label {k} has less than minimum number of vertices, {nmin}, label_counts={label_counts}"
        )

logger.info(["Label counts before enforcing label contiguity:", label_counts])


def enforce_label_contiguity(labels, boundary_vertices, surface, logger=None):
    """
    Enforce spatial contiguity of each label on boundary vertices.
    Keeps the largest component and reassigns strays using surface topology.
    """
    import networkx as nx

    label_map = dict(zip(boundary_vertices, labels))
    new_labels = labels.copy()

    # Build connectivity graph from mesh topology
    connectivity_graph = nx.Graph()
    faces = surface.faces.reshape((-1, 4))[:, 1:4]

    for face in faces:
        for i in range(3):
            for j in range(i + 1, 3):
                u, v = face[i], face[j]
                if u in label_map and v in label_map:
                    connectivity_graph.add_edge(u, v)

    for label in np.unique(labels):
        nodes = [v for v in boundary_vertices if label_map[v] == label]
        subgraph = connectivity_graph.subgraph(nodes)
        components = list(nx.connected_components(subgraph))

        if logger:
            logger.info(f"Label {label}: found {len(components)} connected components")

        if len(components) <= 1:
            continue

        components = sorted(components, key=len, reverse=True)
        keep = components[0]
        reassigned_count = 0

        for comp in components[1:]:
            for v in comp:
                neighbors = list(connectivity_graph.neighbors(v))
                neighbor_labels = [
                    label_map[n] for n in neighbors if label_map[n] != label
                ]
                if neighbor_labels:
                    new_label = max(set(neighbor_labels), key=neighbor_labels.count)
                else:
                    new_label = label  # fallback
                new_labels[np.where(boundary_vertices == v)[0][0]] = new_label
                reassigned_count += 1

        if logger:
            logger.info(f"Label {label}: reassigned {reassigned_count} stray vertices")

    return new_labels


logger.info("Loading surface from GIFTI...")
surface, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
logger.info(f"Surface loaded: {surface.n_points} vertices, {surface.n_faces} faces.")

boundary_vertices = np.where(edges == 1)[0]
labels = enforce_label_contiguity(labels, boundary_vertices, surface, logger=logger)


ap_srcsink = np.zeros((len(edges)), dtype=np.int32)
idx_edges = np.where(edges == 1)[0]
ap_srcsink[idx_edges[labels == 0]] = 1
ap_srcsink[idx_edges[labels == 1]] = 2

pd_srcsink = np.zeros((len(edges)), dtype=np.int32)
pd_srcsink[idx_edges[labels == 2]] = 1
pd_srcsink[idx_edges[labels == 3]] = 2


# save to gifti:
label_dict = {
    "Background": {"key": 0, "red": 1.0, "green": 1.0, "blue": 1.0, "alpha": 0.0},
    "Source": {"key": 1, "red": 1.0, "green": 0.0, "blue": 0.0, "alpha": 1.0},
    "Sink": {"key": 2, "red": 0.0, "green": 0.0, "blue": 1.0, "alpha": 1.0},
}

write_label_gii(
    ap_srcsink, snakemake.output.ap, label_dict=label_dict, metadata=metadata
)
write_label_gii(
    pd_srcsink, snakemake.output.pd, label_dict=label_dict, metadata=metadata
)
