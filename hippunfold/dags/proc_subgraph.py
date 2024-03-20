#!/usr/bin/env python3
import pygraphviz as pgv
import sys
from pathlib import Path

if len(sys.argv) < 2:
    print("DAG should be provided in stdin (i.e. piped in from snakemake or cat)")
    print(f"Usage: {sys.argv[0]}  <output_prefix> ")
    print(f"e.g.: snakemake --dag | {sys.argv[0]} out_dags/mydag")
    sys.exit(0)


import re
import sys
from glob import glob

# first, parse the smk files
pattern = re.compile(r"rule (\w+):")
subgraphs = dict()

for smk in glob("../workflow/rules/*.smk"):
    smk_name = re.findall(re.compile(r"/(\w+).smk"), smk)[0]

    with open(smk, "r") as f:
        text = f.read()

    subgraphs[smk_name] = re.findall(pattern, text)


# ok now we have the rule names
# want to get the nodes with those labels


in_dag = sys.stdin.read()
out_prefix = Path(sys.argv[1])

# make the folder
out_prefix.parents[0].mkdir(parents=True, exist_ok=True)


G = pgv.AGraph(in_dag)

Gwithclusters = G.copy()


for i, cluster_name in enumerate(subgraphs.keys()):
    cluster_nodes = list()
    for n in G.nodes():
        label = n.attr["label"]
        label = label.split("\\")[0]  # strip off wildcards in the label

        if label in subgraphs[cluster_name]:
            cluster_nodes.append(n)

    # create a set of nodes in the subgraph, along with neighbor
    nodes_to_retain = set(cluster_nodes)

    for n in cluster_nodes:
        # add it's input neighbours
        for m in G.in_neighbors(n):
            nodes_to_retain.add(m)

    # print(f'nodes to retain: {nodes_to_retain}')

    # now, create a new graph and remove everything but these nodes
    Gsubset = G.copy()

    for n in Gsubset.nodes():
        # print(n)
        if n not in nodes_to_retain:
            # print(f'node {n} is not in the set, removing it')
            Gsubset.remove_node(n)
    #        else:
    #            print(f'node {n} is in the set, keeping it')

    final_nodes = Gsubset.nodes()
    final_edges = Gsubset.edges()
    # print(f'final_nodes: {final_nodes}')
    # print(f'final_edges: {final_edges}')

    # add the cluster back
    Gsubset.add_subgraph(
        cluster_nodes, name="cluster_0", label=cluster_name, color="blue"
    )

    Gwithclusters.add_subgraph(
        cluster_nodes, name=f"cluster_{i}", label=cluster_name, color="blue"
    )

    if len(Gsubset.nodes()) > 1:
        Gsubset.write(f"{out_prefix}.{cluster_name}.txt")
        #        Gsubset.draw(f"{out_prefix}.{cluster_name}.svg", prog="dot")
        #        Gsubset.draw(f'{out_prefix}.{cluster_name}.pdf',prog='dot')
        Gsubset.draw(f"{out_prefix}.{cluster_name}.png", prog="dot")


Gwithclusters.write(f"{out_prefix}.txt")
# Gwithclusters.draw(f"{out_prefix}.svg", prog="dot")
Gwithclusters.draw(f"{out_prefix}.png", prog="dot")
# Gwithclusters.draw(f'{out_prefix}.pdf',prog='dot')
