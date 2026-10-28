# Author: Kanak Varma
"""
Expansion & bottleneck analysis of a street network (GraphML from fetch_namma.py):

  1. Spectral expansion — algebraic connectivity lambda_2 of the normalized
     Laplacian, Cheeger inequality bounds on the conductance/expansion constant
     h(G):  lambda_2/2  <=  h(G)  <=  sqrt(2*lambda_2),
     and the Ramanujan comparison: a d-regular Ramanujan graph has all nontrivial
     adjacency eigenvalues <= 2*sqrt(d-1); for our irregular graph we compare the
     normalized second eigenvalue mu_2 = 1 - lambda_2 against 2*sqrt(d-1)/d at d = d_avg.
  2. Sweep cut on the Fiedler vector — the constructive half of Cheeger's theorem:
     exhibits the (approximately) minimum-conductance cut = the global bottleneck.
  3. Max-flow / min-cut (Menger): fewest unit-capacity street segments whose removal
     disconnects the north band from the south band (and east from west) — the
     literal bottleneck roads.

Usage:  python expander_bottleneck.py [path/to/graph.graphml]
Writes figures/ and results/ next to the cwd.
"""
import gc
import os
import sys
import time

import networkx as nx
import numpy as np
import osmnx as ox
import scipy.sparse.linalg as spla
from networkx.algorithms.flow import preflow_push
from scipy.sparse import csgraph

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
tag = os.path.splitext(os.path.basename(GRAPH))[0]
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------- load light
log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
xs = {u: d["x"] for u, d in G.nodes(data=True)}
ys = {u: d["y"] for u, d in G.nodes(data=True)}
S = nx.Graph()                                # lightweight simple undirected copy
S.add_nodes_from(G.nodes)
S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
del G
gc.collect()
S = S.subgraph(max(nx.connected_components(S), key=len)).copy()
nodes = list(S.nodes)
n, m = S.number_of_nodes(), S.number_of_edges()
x = np.array([xs[u] for u in nodes])
y = np.array([ys[u] for u in nodes])
d_avg = 2 * m / n
log(f"largest component (simple, undirected): n={n:,}  m={m:,}  d_avg={d_avg:.3f}")

# ---------------------------------------------------------------- spectrum
A = nx.to_scipy_sparse_array(S, nodelist=nodes, weight=None, format="csr", dtype=np.float64)
deg = np.asarray(A.sum(axis=1)).ravel()
L = csgraph.laplacian(A, normed=True).tocsc()

log("computing lambda_2 of normalized Laplacian (shift-invert Lanczos) ...")
vals, vecs = spla.eigsh(L, k=2, sigma=-1e-2, which="LM")
lam2 = float(vals[1])
log(f"lambda_2 = {lam2:.3e}")
log("computing lambda_max ...")
lam_max = float(spla.eigsh(L, k=1, which="LA", return_eigenvectors=False)[0])

# sweep vector for conductance: v / sqrt(deg)
fiedler = vecs[:, 1] / np.sqrt(np.maximum(deg, 1))

# ---------------------------------------------------------------- sweep cut
log("sweep cut over Fiedler vector ...")
order = np.argsort(fiedler)
rank = np.empty(n, dtype=np.int64)
rank[order] = np.arange(n)
indptr, indices = A.indptr, A.indices
vol_tot = float(deg.sum())
cut = 0.0
vol = 0.0
phi = np.full(n - 1, np.inf)
for step, i in enumerate(order[:-1]):
    nbrs = indices[indptr[i]:indptr[i + 1]]
    inside = np.count_nonzero(rank[nbrs] <= step)      # neighbours already swept in
    cut += len(nbrs) - 2 * inside
    vol += deg[i]
    phi[step] = cut / min(vol, vol_tot - vol)

best = int(np.argmin(phi))
phi_best = float(phi[best])
# best cut with at least 10% of volume on the small side (a *regional* bottleneck)
volcum = np.cumsum(deg[order])
balanced = (volcum[:-1] >= 0.10 * vol_tot) & (volcum[:-1] <= 0.90 * vol_tot)
best_bal = int(np.flatnonzero(balanced)[np.argmin(phi[balanced])])
phi_bal = float(phi[best_bal])

side_bal = np.zeros(n, dtype=bool)
side_bal[order[:best_bal + 1]] = True                   # True = side A
side_nodes = {nodes[i] for i in np.flatnonzero(side_bal)}
sweep_cut_edges = [(u, v) for u, v in nx.edge_boundary(S, side_nodes)]
log(f"best conductance phi = {phi_best:.3e}; balanced (>=10% vol) phi = {phi_bal:.3e} "
    f"with {len(sweep_cut_edges)} cut edges")

# ---------------------------------------------------------------- min cuts
nx.set_edge_attributes(S, 1, "capacity")

def regional_min_cut(coord, label):
    """Max-flow min-cut between the top and bottom 10% bands of `coord`."""
    hi, lo = np.quantile(coord, 0.9), np.quantile(coord, 0.1)
    src_set = {nodes[i] for i in np.flatnonzero(coord >= hi)}
    snk_set = {nodes[i] for i in np.flatnonzero(coord <= lo)}
    F = S.copy()
    F.add_edges_from((("__SRC__", u) for u in src_set), capacity=np.inf)
    F.add_edges_from((("__SNK__", u) for u in snk_set), capacity=np.inf)
    t0 = time.time()
    cutval, (P1, P2) = nx.minimum_cut(F, "__SRC__", "__SNK__", flow_func=preflow_push)
    P1 = set(P1) - {"__SRC__"}
    cut_edges = [(u, v) for u, v in nx.edge_boundary(S, P1)]
    log(f"{label}: min cut = {int(cutval)} road segments  [{time.time()-t0:.0f}s]")
    return int(cutval), cut_edges, src_set, snk_set

ns_val, ns_edges, north, south = regional_min_cut(y, "north<->south")
ew_val, ew_edges, east, west = regional_min_cut(x, "east<->west")

# ---------------------------------------------------------------- report
mu2 = 1.0 - lam2                               # 2nd eigenvalue of normalized adjacency
ram = 2 * np.sqrt(d_avg - 1) / d_avg           # Ramanujan threshold, normalized
report = f"""
================= {tag} — expansion & bottlenecks =================
Largest component: n = {n:,} nodes, m = {m:,} undirected edges, d_avg = {d_avg:.3f}

--- Spectral expansion ---
lambda_2 (normalized Laplacian)        : {lam2:.3e}
lambda_max                             : {lam_max:.4f}   (2.0 would mean bipartite)
Cheeger bounds:  {lam2/2:.3e}  <=  h(G)  <=  {np.sqrt(2*lam2):.3e}
Best sweep-cut conductance phi         : {phi_best:.3e}  (certified h(G) <= this)
Balanced bottleneck (>=10% vol side)   : phi = {phi_bal:.3e}, {len(sweep_cut_edges)} edges cut,
                                         small side = {min(best_bal+1, n-best_bal-1):,} nodes
--- Expander / Ramanujan verdict ---
An expander family needs lambda_2 bounded away from 0; a (d_avg-)regular Ramanujan
graph would need mu_2 = 1 - lambda_2 <= 2*sqrt(d-1)/d.
  mu_2            = {mu2:.6f}
  Ramanujan bound = {ram:.6f}
  -> mu_2 exceeds the bound by a factor {(mu2/ram):.4f}; lambda_2 ~ {lam2:.1e} vs the
     ~{1-ram:.3f} an expander of this degree would need: a gap of ~{(1-ram)/lam2:,.0f}x.
VERDICT: decisively NOT an expander and far from Ramanujan — as expected for a
near-planar graph (O(sqrt n) separators), random-walk mixing is ~1/lambda_2 ~ {1/lam2:,.0f} steps.

--- Max-flow / min-cut bottlenecks (unit capacity per street segment) ---
North <-> South: cutting just {ns_val} road segments disconnects the top-10% latitude
                 band from the bottom-10% band.
East  <-> West : cutting {ew_val} segments disconnects east from west.
(sqrt(n) = {np.sqrt(n):.0f} — Lipton–Tarjan-scale separators, tiny vs m = {m:,}.)
"""
print(report)
with open(f"results/{tag}_expander_report.txt", "w") as f:
    f.write(report)

# ---------------------------------------------------------------- figures
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

BLUE, ORANGE, INK, GRID = "#3987e5", "#ff8a3d", "#f0efec", "#3a3a37"
SURFACE = "#1a1a19"
pos_i = {u: i for i, u in enumerate(nodes)}

def seg_array(edges):
    return np.array([[[x[pos_i[u]], y[pos_i[u]]], [x[pos_i[v]], y[pos_i[v]]]] for u, v in edges])

all_segs = seg_array(list(S.edges))

def base_ax():
    fig, ax = plt.subplots(figsize=(10, 10), dpi=220)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.add_collection(LineCollection(all_segs, colors=GRID, linewidths=0.25, zorder=1))
    ax.set_aspect(1 / np.cos(np.deg2rad(y.mean())))     # ~equal metres in lon/lat
    ax.autoscale()
    ax.set_axis_off()
    return fig, ax

log("plotting spectral partition ...")
fig, ax = base_ax()
ax.scatter(x[side_bal], y[side_bal], s=0.25, c=BLUE, lw=0, zorder=2)
ax.scatter(x[~side_bal], y[~side_bal], s=0.25, c=ORANGE, lw=0, zorder=2)
ax.add_collection(LineCollection(seg_array(sweep_cut_edges), colors=INK, linewidths=1.6, zorder=3))
ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", lw=0, label="side A"),
                   Line2D([], [], color=ORANGE, marker="o", lw=0, label="side B"),
                   Line2D([], [], color=INK, lw=2, label=f"bottleneck cut ({len(sweep_cut_edges)} edges)")],
          loc="lower left", fontsize=9, frameon=False, labelcolor=INK)
ax.set_title(f"Spectral (Fiedler) bottleneck — {tag}\nconductance φ = {phi_bal:.2e}", fontsize=11, color=INK)
fig.savefig(f"figures/{tag}_spectral_bottleneck.png", bbox_inches="tight", facecolor=SURFACE)
plt.close(fig)

for label, (cv, ce, sa, sb), fname in [
        ("North ↔ South", (ns_val, ns_edges, north, south), f"figures/{tag}_mincut_NS.png"),
        ("East ↔ West", (ew_val, ew_edges, east, west), f"figures/{tag}_mincut_EW.png")]:
    log(f"plotting {fname} ...")
    ia = [pos_i[u] for u in sa]; ib = [pos_i[u] for u in sb]
    fig, ax = base_ax()
    ax.scatter(x[ia], y[ia], s=0.3, c=BLUE, lw=0, zorder=2)
    ax.scatter(x[ib], y[ib], s=0.3, c=ORANGE, lw=0, zorder=2)
    ax.add_collection(LineCollection(seg_array(ce), colors=INK, linewidths=1.6, zorder=3))
    ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", lw=0, label="source band"),
                       Line2D([], [], color=ORANGE, marker="o", lw=0, label="sink band"),
                       Line2D([], [], color=INK, lw=2, label=f"min cut = {cv} segments")],
              loc="lower left", fontsize=9, frameon=False, labelcolor=INK)
    ax.set_title(f"Max-flow min-cut bottleneck, {label} — {tag}", fontsize=11, color=INK)
    fig.savefig(fname, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

log(f"done. results/{tag}_expander_report.txt + 3 figures")
