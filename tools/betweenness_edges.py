# Author: Kanak Varma
"""
Edge-painted betweenness map: every street segment colored by the mean
betweenness centrality of its two endpoint nodes (no node markers).

Usage:  python betweenness_edges.py [path/to/graph.graphml]
Node betweenness is cached in results/<tag>_node_bc.npz (reused if present).
"""
import os
import sys
import time

import networkx as nx
import numpy as np
import osmnx as ox

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
K_NODE_BC = 300
SEED = 42
tag = os.path.splitext(os.path.basename(GRAPH))[0]
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
Gs = ox.truncate.largest_component(G, strongly=True)
log(f"SCC: {Gs.number_of_nodes():,} nodes, {Gs.number_of_edges():,} edges")

bc_file = f"results/{tag}_node_bc.npz"
if os.path.exists(bc_file):
    z = np.load(bc_file, allow_pickle=True)
    bc = dict(zip(z["nodes"].tolist(), z["bc"]))
    log(f"loaded cached node betweenness from {bc_file}")
else:
    D = ox.convert.to_digraph(Gs, weight="length")
    log(f"computing node betweenness (k={K_NODE_BC}, length-weighted) ... ~7 min")
    bc = nx.betweenness_centrality(D, k=min(K_NODE_BC, D.number_of_nodes()),
                                   weight="length", normalized=True, seed=SEED)
    np.savez_compressed(bc_file, nodes=np.array(list(D.nodes), dtype=object),
                        bc=np.array([bc[u] for u in D.nodes]))
    log(f"cached -> {bc_file}")

import matplotlib as mpl
import matplotlib.pyplot as plt

BLUES = mpl.colors.LinearSegmentedColormap.from_list("blues1", ["#dbe9f6", "#2171b5", "#08306b"])

log("rendering edge-painted betweenness map ...")
ev = np.array([(bc.get(u, 0.0) + bc.get(v, 0.0)) / 2 for u, v, _ in Gs.edges(keys=True)])
s = np.sqrt(ev / max(ev.max(), 1e-12))        # sqrt scale so mid-range roads stay visible
ec = [BLUES(x) for x in s]
ew = list(0.25 + 2.2 * s)
fig, ax = ox.plot_graph(Gs, bgcolor="white", node_size=0, edge_color=ec, edge_linewidth=ew,
                        show=False, close=False, figsize=(12, 12))
sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, ev.max()), cmap=BLUES)
cb = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.01)
cb.set_label("betweenness centrality (mean of edge endpoints)", color="#444")
cb.ax.tick_params(colors="#444", labelsize=7)
ax.set_title(f"Betweenness centrality, edge-painted — {tag}", color="#222", fontsize=12)
fig.savefig(f"figures/{tag}_betweenness_edges.png", dpi=220, bbox_inches="tight", facecolor="white")
plt.close(fig)
log(f"done: figures/{tag}_betweenness_edges.png")
