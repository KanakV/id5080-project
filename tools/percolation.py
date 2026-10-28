# Author: Kanak Varma
"""
Robustness / site-percolation analysis of the street network.

Removes nodes (intersections) one at a time under four strategies and tracks the
relative size S(q) of the giant connected component as a function of the removed
fraction q — computed exactly for every q with a reverse Newman–Ziff sweep
(union–find while re-adding nodes in reverse order), so each full curve costs
O((n+m) alpha).

Strategies:
  random            — uniformly random order (averaged over 5 seeds)
  degree (static)   — descending initial degree
  degree (adaptive) — recompute degrees as nodes are removed (lazy-heap peeling)
  betweenness       — descending approximate betweenness (cached, k=300, static)

Usage:  python percolation.py [path/to/graph.graphml]
Writes results/<tag>_percolation.npz, results/<tag>_percolation_report.txt,
figures/<tag>_percolation.png
"""
import gc
import heapq
import os
import sys
import time

import networkx as nx
import numpy as np
import osmnx as ox

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
N_RANDOM_SEEDS = 5
SEED = 42
tag = os.path.splitext(os.path.basename(GRAPH))[0]
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
S = nx.Graph()
S.add_nodes_from(G.nodes)
S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
del G; gc.collect()
S = S.subgraph(max(nx.connected_components(S), key=len)).copy()
nodes = list(S.nodes)
n = len(nodes)
idx = {u: i for i, u in enumerate(nodes)}
A = nx.to_scipy_sparse_array(S, nodelist=nodes, weight=None, format="csr")
indptr, indices = A.indptr, A.indices
deg0 = np.diff(indptr)
log(f"GCC: n={n:,}, m={S.number_of_edges():,}")


def giant_curve(removal_order):
    """removal_order: array of node indices, first-removed first.
    Returns S_rel[k] = giant component fraction after removing the first k nodes."""
    parent = np.arange(n)
    size = np.ones(n, dtype=np.int64)
    present = np.zeros(n, dtype=bool)

    def find(a):
        root = a
        while parent[root] != root:
            root = parent[root]
        while parent[a] != root:
            parent[a], a = root, parent[a]
        return root

    best = 0
    S_rel = np.empty(n + 1)
    S_rel[n] = 0.0
    for pos in range(n - 1, -1, -1):          # re-add in reverse removal order
        v = removal_order[pos]
        present[v] = True
        for w in indices[indptr[v]:indptr[v + 1]]:
            if present[w]:
                ra, rb = find(v), find(w)
                if ra != rb:
                    if size[ra] < size[rb]:
                        ra, rb = rb, ra
                    parent[rb] = ra
                    size[ra] += size[rb]
        best = max(best, size[find(v)])
        S_rel[pos] = best / n
    return S_rel


orders = {}

rng = np.random.default_rng(SEED)
log("orders: random x5, degree-static, betweenness-static ...")
rand_curves = []
for s in range(N_RANDOM_SEEDS):
    o = np.random.default_rng(SEED + s).permutation(n)
    rand_curves.append(giant_curve(o))
orders["random"] = None  # averaged below

orders["degree (static)"] = np.argsort(-deg0, kind="stable")

bc_file = f"results/{tag}_node_bc.npz"
if os.path.exists(bc_file):
    z = np.load(bc_file, allow_pickle=True)
    bcmap = dict(zip(z["nodes"].tolist(), z["bc"]))
    bcv = np.array([bcmap.get(u, 0.0) for u in nodes])
    orders["betweenness (static)"] = np.argsort(-bcv, kind="stable")
else:
    log(f"  WARNING: {bc_file} missing — skipping betweenness strategy")

log("order: degree-adaptive (lazy-heap peeling) ...")
livedeg = deg0.copy()
alive = np.ones(n, dtype=bool)
heap = [(-int(d), i) for i, d in enumerate(livedeg)]
heapq.heapify(heap)
adaptive = np.empty(n, dtype=np.int64)
k = 0
while heap:
    d, v = heapq.heappop(heap)
    if not alive[v] or -d != livedeg[v]:
        continue
    alive[v] = False
    adaptive[k] = v; k += 1
    for w in indices[indptr[v]:indptr[v + 1]]:
        if alive[w]:
            livedeg[w] -= 1
            heapq.heappush(heap, (-int(livedeg[w]), w))
orders["degree (adaptive)"] = adaptive

log("computing percolation curves ...")
curves = {"random": np.mean(rand_curves, axis=0)}
for name, o in orders.items():
    if o is not None and name != "random":
        curves[name] = giant_curve(o)

q = np.arange(n + 1) / n
qc = {}
for name, c in curves.items():
    below = np.flatnonzero(c < 0.01)
    qc[name] = q[below[0]] if len(below) else float("nan")

np.savez_compressed(f"results/{tag}_percolation.npz", q=q,
                    **{k.replace(" ", "_"): v for k, v in curves.items()})

report = [f"================ {tag} — robustness / site percolation ================",
          f"Giant-component collapse threshold q_c (first q with S < 1% of n):"]
for name, v in qc.items():
    report.append(f"  {name:22s}: q_c = {v:.3f}")
report.append("(random vs targeted gap = classic robust-yet-fragile signature)")
report = "\n".join(report) + "\n"
print(report)
with open(f"results/{tag}_percolation_report.txt", "w") as f:
    f.write(report)

# ---------------------------------------------------------------- figure (chart: light theme)
import matplotlib.pyplot as plt

COLORS = {"random": "#2a78d6", "degree (static)": "#eb6834",
          "degree (adaptive)": "#1baf7a", "betweenness (static)": "#eda100"}
MARKERS = {"random": "o", "degree (static)": "s",
           "degree (adaptive)": "^", "betweenness (static)": "D"}
INK, SUB = "#1a1a19", "#5f5e56"

fig, ax = plt.subplots(figsize=(8, 6), dpi=200)
step = max(1, n // 400)
for name in COLORS:
    if name not in curves:
        continue
    ax.plot(q[::step], curves[name][::step], color=COLORS[name], lw=1.6,
            marker=MARKERS[name], ms=4, markevery=25, label=name)
ax.set_xlabel("fraction of nodes removed  q", color=INK)
ax.set_ylabel("giant component share  S(q)", color=INK)
ax.set_xlim(0, max(min(1.0, max(v for v in qc.values() if v == v) + 0.1), 0.5))
ax.set_ylim(0, 1.0)
ax.set_title("Robustness under node removal (site percolation) — " + tag,
             color=INK, fontsize=12)
ax.legend(frameon=False, labelcolor=INK)
ax.grid(True, color="#e8e7e0", lw=0.6)
ax.set_axisbelow(True)
ax.tick_params(colors=SUB)
for s_ in ax.spines.values():
    s_.set_color("#d4d3c9")
fig.savefig(f"figures/{tag}_percolation.png", bbox_inches="tight", facecolor="white")
plt.close(fig)
log(f"done: figures/{tag}_percolation.png, results/{tag}_percolation_report.txt")
