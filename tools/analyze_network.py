# Author: Kanak Varma
"""
Basic structural analysis of a city street network (GraphML from fetch_namma.py):

  1. Betweenness centrality  — which intersections carry the most shortest paths
     (exact for small graphs, k-source sampled approximation for big ones).
  2. Cyclic structure        — how much of the network lies on cycles (loops/grid)
     vs. tree-like parts (dead ends, bridges): cyclomatic number, alpha/beta/gamma
     indices, bridge share, 2-core, biconnected components.

Usage:
    python analyze_network.py [path/to/graph.graphml] [k_samples]

Writes figures to ./figures/ and a top-nodes CSV to ./results/.
"""
import csv
import os
import sys
import time

import networkx as nx
import numpy as np
import osmnx as ox

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
K_SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 300   # sources for approximate betweenness
EXACT_LIMIT = 4000                                           # exact betweenness below this many nodes
SEED = 42

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
tag = os.path.splitext(os.path.basename(GRAPH))[0]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def edge_names(G, u):
    """Street names on edges incident to node u."""
    names = set()
    for _, _, d in list(G.edges(u, data=True)) + list(G.in_edges(u, data=True)):
        nm = d.get("name")
        if isinstance(nm, list):
            names.update(nm)
        elif nm:
            names.add(nm)
    return "; ".join(sorted(names)) if names else "(unnamed)"


# ---------------------------------------------------------------- load
log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
n, m = G.number_of_nodes(), G.number_of_edges()
log(f"full graph: {n:,} nodes, {m:,} directed edges")

# keep the largest strongly connected component (betweenness needs reachability)
Gs = ox.truncate.largest_component(G, strongly=True)
log(f"largest strongly connected component: {Gs.number_of_nodes():,} nodes "
    f"({Gs.number_of_nodes()/n:.1%}), {Gs.number_of_edges():,} edges")

# ---------------------------------------------------------------- betweenness
D = ox.convert.to_digraph(Gs, weight="length")   # collapse parallel edges, keep shortest
nD = D.number_of_nodes()
if nD <= EXACT_LIMIT:
    log(f"computing EXACT betweenness centrality on {nD:,} nodes ...")
    bc = nx.betweenness_centrality(D, weight="length", normalized=True)
    bc_kind = "exact"
else:
    k = min(K_SAMPLES, nD)
    log(f"computing APPROXIMATE betweenness centrality ({k} sampled sources, "
        f"length-weighted) on {nD:,} nodes ... this is the slow step")
    bc = nx.betweenness_centrality(D, k=k, weight="length", normalized=True, seed=SEED)
    bc_kind = f"approx, k={k}"
log("betweenness done")

vals = np.array(list(bc.values()))
top = sorted(bc.items(), key=lambda kv: -kv[1])[:20]

with open(f"results/{tag}_top_betweenness.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "osmid", "betweenness", "lat", "lon", "streets"])
    for r, (node, b) in enumerate(top, 1):
        w.writerow([r, node, f"{b:.6f}", f"{Gs.nodes[node]['y']:.6f}",
                    f"{Gs.nodes[node]['x']:.6f}", edge_names(Gs, node)])

# ---------------------------------------------------------------- cyclic structure
log("analysing cyclic structure ...")
U = ox.convert.to_undirected(G)          # undirected MultiGraph of the full network
Un, Um = U.number_of_nodes(), U.number_of_edges()
c = nx.number_connected_components(U)

mu = Um - Un + c                          # cyclomatic number = dimension of cycle space
alpha = mu / (2 * Un - 5)                 # fraction of the max possible independent cycles (planar)
beta = Um / Un                            # edges per node
gamma = Um / (3 * (Un - 2))               # fraction of max possible edges (planar)

S = nx.Graph(U)                           # simple graph for bridge/core analysis
S.remove_edges_from(nx.selfloop_edges(S))  # self-loops break k_core; they're cycles by definition
multi_uv = {frozenset((u, v)) for u, v, k in U.edges(keys=True) if k > 0}
bridges = {frozenset(e) for e in nx.bridges(S)} - multi_uv   # parallel edges are never bridges

bridge_len = cycle_len = 0.0
for u, v, d in U.edges(data=True):
    L = float(d.get("length", 0.0))
    if frozenset((u, v)) in bridges:
        bridge_len += L
    else:
        cycle_len += L
total_len = bridge_len + cycle_len

core2 = nx.k_core(S, 2)                   # the "cyclic skeleton": strip dead-end trees
bicomp = [b for b in nx.biconnected_components(S) if len(b) > 2]
largest_bi = max((len(b) for b in bicomp), default=0)

report = f"""
================= {tag} =================
Nodes (intersections/dead-ends) : {Un:,}
Edges (street segments, undir.) : {Um:,}
Connected components            : {c}
Total street length             : {total_len/1000:,.1f} km

--- Betweenness centrality ({bc_kind}, length-weighted, on largest SCC) ---
mean = {vals.mean():.5f}   max = {vals.max():.5f}   (normalized)
Gini-like concentration: top 1% of nodes hold {np.sort(vals)[::-1][:max(1,len(vals)//100)].sum()/vals.sum():.1%} of total betweenness
Top 10 intersections:
""" + "\n".join(
    f"  {r:2d}. bc={b:.5f}  ({Gs.nodes[nd]['y']:.5f}, {Gs.nodes[nd]['x']:.5f})  {edge_names(Gs, nd)[:70]}"
    for r, (nd, b) in enumerate(top[:10], 1)
) + f"""

--- Cyclic structure (undirected) ---
Cyclomatic number mu = m - n + c       : {mu:,}   (independent cycles)
alpha index  mu/(2n-5)                 : {alpha:.3f}
beta  index  m/n                       : {beta:.3f}
gamma index  m/3(n-2)                  : {gamma:.3f}
Bridge edges (on NO cycle)             : {len(bridges):,} ({len(bridges)/Um:.1%} of edges, {bridge_len/1000:,.1f} km = {bridge_len/total_len:.1%} of length)
Edges on >=1 cycle                     : {Um-len(bridges):,} ({cycle_len/1000:,.1f} km = {cycle_len/total_len:.1%} of length)
2-core (cyclic skeleton) nodes         : {core2.number_of_nodes():,} ({core2.number_of_nodes()/Un:.1%} of nodes)
Biconnected blocks (>2 nodes)          : {len(bicomp):,}, largest spans {largest_bi:,} nodes ({largest_bi/Un:.1%})
"""
print(report)
with open(f"results/{tag}_report.txt", "w") as f:
    f.write(report)

# ---------------------------------------------------------------- figures
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

log("plotting betweenness map ...")
# single-hue sequential ramp (light -> dark blue) on sqrt scale; faint gray base net
bvals = np.array([bc[nd] for nd in Gs.nodes])
norm = mpl.colors.Normalize(vmin=0, vmax=np.sqrt(bvals.max()) or 1)
cmap = mpl.colors.LinearSegmentedColormap.from_list("blues1", ["#c6dbef", "#2171b5", "#08306b"])
nc = [cmap(norm(np.sqrt(b))) for b in bvals]
ns = list(0.15 + 28 * np.sqrt(bvals / (bvals.max() or 1)))
fig, ax = ox.plot_graph(Gs, bgcolor="white", node_color=nc, node_size=ns, node_zorder=2,
                        edge_color="#d9d9d9", edge_linewidth=0.3, show=False, close=False)
sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, bvals.max()), cmap=cmap)
cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.01)
cb.set_label("node betweenness centrality", color="#444444")
cb.ax.tick_params(colors="#444444", labelsize=7)
ax.set_title(f"Betweenness centrality — {tag} ({bc_kind})", color="#222222", fontsize=11)
fig.savefig(f"figures/{tag}_betweenness.png", dpi=250, bbox_inches="tight", facecolor="white")
plt.close(fig)

log("plotting cyclic-vs-tree map ...")
GRAY, ORANGE = "#b5b5b5", "#E8710A"      # cycle edges recessive gray, bridges pop orange
ec, ew = [], []
for u, v, k in U.edges(keys=True):
    is_br = frozenset((u, v)) in bridges
    ec.append(ORANGE if is_br else GRAY)
    ew.append(0.9 if is_br else 0.35)
fig, ax = ox.plot_graph(U, bgcolor="white", node_size=0, edge_color=ec, edge_linewidth=ew,
                        show=False, close=False)
ax.legend(handles=[Line2D([], [], color=GRAY, lw=2, label="on a cycle (loops/grid)"),
                   Line2D([], [], color=ORANGE, lw=2, label="bridge / tree-like (no cycle)")],
          loc="lower left", fontsize=8, frameon=False, labelcolor="#333333")
ax.set_title(f"Cyclic vs tree-like street segments — {tag}", color="#222222", fontsize=11)
fig.savefig(f"figures/{tag}_cycles.png", dpi=250, bbox_inches="tight", facecolor="white")
plt.close(fig)

log(f"done. figures/{tag}_*.png, results/{tag}_report.txt, results/{tag}_top_betweenness.csv")
