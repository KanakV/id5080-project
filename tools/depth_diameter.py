# Author: Kanak Varma
"""
Eccentricity ("depth") map, diameter, and degree-exponent check for the street network.

  1. Eccentricity of every node (hops): tight lower bound from 64 farthest-point
     landmark BFS sweeps — for near-planar graphs the landmarks hit the extremal
     "corners", so the bound is essentially exact for mapping purposes.
  2. Diameter: lower bound from the sweeps, then an exact computation with
     NetworkX's extrema-bounding (may take a while; runs last).
  3. gamma check vs the AB (Barabasi-Albert) model: discrete power-law MLE
     (Clauset et al. approximation) on the degree distribution. AB predicts
     gamma = 3; a street network has a narrow physical degree range, so the fit
     mostly demonstrates that a power law does not apply.

Usage:  python depth_diameter.py [path/to/graph.graphml]
"""
import gc
import os
import sys
import time

import networkx as nx
import numpy as np
import osmnx as ox

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
N_LANDMARKS = 64
SEED = 42
tag = os.path.splitext(os.path.basename(GRAPH))[0]
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
xs = {u: d["x"] for u, d in G.nodes(data=True)}
ys = {u: d["y"] for u, d in G.nodes(data=True)}
S = nx.Graph()
S.add_nodes_from(G.nodes)
S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
del G; gc.collect()
S = S.subgraph(max(nx.connected_components(S), key=len)).copy()
nodes = list(S.nodes)
n = len(nodes)
idx = {u: i for i, u in enumerate(nodes)}
log(f"GCC: n={n:,}, m={S.number_of_edges():,}")

# ---------------------------------------------------------------- gamma (AB check)
degs = np.array([d for _, d in S.degree()])
for dmin in (1, 2, 3):
    sel = degs[degs >= dmin]
    gamma = 1.0 + len(sel) / np.log(sel / (dmin - 0.5)).sum()
    log(f"power-law MLE with d_min={dmin}: gamma_hat = {gamma:.2f}  (n_tail={len(sel):,})")
gamma_report = (
    f"Degree exponent check vs AB model:\n"
    f"  degrees span only {degs.min()}..{degs.max()} (physical limit of road junctions)\n"
    f"  discrete MLE gamma_hat = {1.0 + len(degs)/np.log(degs/0.5).sum():.2f} (d_min=1), "
    f"{1.0 + (degs>=3).sum()/np.log(degs[degs>=3]/2.5).sum():.2f} (d_min=3)\n"
    f"  AB predicts gamma = 3 over an unbounded tail -> NOT applicable: the empirical\n"
    f"  distribution is narrow/exponential-bounded, not scale-free.\n"
)

# ---------------------------------------------------------------- landmark BFS
log(f"eccentricity bounds: {N_LANDMARKS} farthest-point landmark BFS sweeps ...")
rng = np.random.default_rng(SEED)
ecc_lb = np.zeros(n, dtype=np.int32)
mindist = np.full(n, np.iinfo(np.int32).max, dtype=np.int32)
lm = int(rng.integers(n))
diam_pair = (None, None, 0)
for it in range(N_LANDMARKS):
    dist = nx.single_source_shortest_path_length(S, nodes[lm])
    dv = np.zeros(n, dtype=np.int32)
    for u, d in dist.items():
        dv[idx[u]] = d
    np.maximum(ecc_lb, dv, out=ecc_lb)
    far = int(dv.argmax())
    if dv[far] > diam_pair[2]:
        diam_pair = (nodes[lm], nodes[far], int(dv[far]))
    np.minimum(mindist, dv, out=mindist)
    nxt = int(mindist.argmax())
    if it % 16 == 15:
        log(f"  {it+1}/{N_LANDMARKS} sweeps, diameter lower bound = {ecc_lb.max()}")
    lm = nxt

diam_lb = int(ecc_lb.max())
log(f"landmark phase done: diameter >= {diam_lb} hops "
    f"(between {diam_pair[0]} and {diam_pair[1]})")
np.savez_compressed(f"results/{tag}_eccentricity.npz",
                    nodes=np.array(nodes, dtype=object), ecc_lb=ecc_lb)

# ---------------------------------------------------------------- figure
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

x = np.array([xs[u] for u in nodes]); y = np.array([ys[u] for u in nodes])
segs = np.array([[[xs[u], ys[u]], [xs[v], ys[v]]] for u, v in S.edges])
CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "blues1", ["#dbe9f6", "#2a78d6", "#08306b"])

fig, ax = plt.subplots(figsize=(11, 11), dpi=220)
ax.add_collection(LineCollection(segs, colors="#e3e2da", linewidths=0.25, zorder=1))
sc = ax.scatter(x, y, c=ecc_lb, cmap=CMAP, s=0.5, lw=0, zorder=2)
ax.set_aspect(1 / np.cos(np.deg2rad(y.mean())))
ax.autoscale(); ax.set_axis_off()
cb = fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.01)
cb.set_label("eccentricity (hops) — max BFS depth from the node", color="#444")
cb.ax.tick_params(colors="#444", labelsize=8)
ax.set_title(f"Node eccentricity (graph depth) — {tag}\n"
             f"{N_LANDMARKS}-landmark bound; diameter ≥ {diam_lb} hops",
             color="#1a1a19", fontsize=12)
fig.savefig(f"figures/{tag}_eccentricity.png", bbox_inches="tight", facecolor="white")
plt.close(fig)
log(f"figure saved: figures/{tag}_eccentricity.png")

report = (f"================ {tag} — depth / diameter / gamma ================\n"
          f"{gamma_report}"
          f"Eccentricity (hops): min={ecc_lb.min()} (graph center), "
          f"median={int(np.median(ecc_lb))}, max={diam_lb}\n"
          f"Diameter lower bound: {diam_lb} hops\n")
print(report)
with open(f"results/{tag}_depth_report.txt", "w") as f:
    f.write(report)

# ---------------------------------------------------------------- exact diameter
log("computing EXACT diameter with extrema bounding (this can take a while) ...")
t0 = time.time()
D = nx.diameter(S, usebounds=True)
line = f"EXACT diameter = {D} hops  [{time.time()-t0:.0f}s]"
log(line)
with open(f"results/{tag}_depth_report.txt", "a") as f:
    f.write(line + "\n")
