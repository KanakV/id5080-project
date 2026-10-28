# Author: Kanak Varma
"""
Percolation comparison: the city vs ER and WS (dark-theme chart).

City curves are reused from results/<tag>_percolation.npz (percolation.py);
n and m come from results/<tag>_model_comparison.csv (model_comparison.py),
so run percolation.py and model_comparison.py first.
ER G(n,m) and WS(k=4, p=0.1) are generated with the same n (and m for ER) and
seed as model_comparison.py, then run through the same reverse Newman-Ziff sweep
for two strategies: random failure (5-seed average) and adaptive degree attack.
NOTE: WS is ~52% denser (k must be an even integer), which inflates its robustness.

Usage:  python percolation_models.py [path/to/graph.graphml]
Writes figures/<tag>_percolation_vs_models_dark.png and appends q_c values to
results/<tag>_percolation_report.txt.
"""
import csv
import gc
import heapq
import os
import time

import sys

import networkx as nx
import numpy as np

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
TAG = os.path.splitext(os.path.basename(GRAPH))[0]
CITY = TAG.split("_")[0].capitalize()
with open(f"results/{TAG}_model_comparison.csv") as f:   # city GCC size from step 3
    row = next(r for r in csv.DictReader(f) if r["graph"] == CITY)
N_TARGET, M_TARGET = int(row["n_gcc"]), int(row["m_gcc"])
SEED = 42
N_RANDOM_SEEDS = 5
os.makedirs("figures", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def to_csr(G):
    import scipy.sparse  # noqa: F401  (networkx needs scipy for this)
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    A = nx.to_scipy_sparse_array(G, weight=None, format="csr")
    return A.indptr, A.indices, A.shape[0]


def giant_curve(indptr, indices, n, removal_order):
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
    S = np.empty(n + 1)
    S[n] = 0.0
    for pos in range(n - 1, -1, -1):
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
        S[pos] = best / n
    return S


def adaptive_degree_order(indptr, indices, n):
    livedeg = np.diff(indptr).copy()
    alive = np.ones(n, dtype=bool)
    heap = [(-int(d), i) for i, d in enumerate(livedeg)]
    heapq.heapify(heap)
    order = np.empty(n, dtype=np.int64)
    k = 0
    while heap:
        d, v = heapq.heappop(heap)
        if not alive[v] or -d != livedeg[v]:
            continue
        alive[v] = False
        order[k] = v; k += 1
        for w in indices[indptr[v]:indptr[v + 1]]:
            if alive[w]:
                livedeg[w] -= 1
                heapq.heappush(heap, (-int(livedeg[w]), w))
    return order


def curves_for(indptr, indices, n, label):
    log(f"[{label}] random x{N_RANDOM_SEEDS} + adaptive degree (n={n:,}) ...")
    rand = np.mean([giant_curve(indptr, indices, n,
                                np.random.default_rng(SEED + s).permutation(n))
                    for s in range(N_RANDOM_SEEDS)], axis=0)
    adapt = giant_curve(indptr, indices, n, adaptive_degree_order(indptr, indices, n))
    return {"random": rand, "degree (adaptive)": adapt}


def qc(q, c):
    below = np.flatnonzero(c < 0.01)
    return q[below[0]] if len(below) else float("nan")


# ---------------------------------------------------------------- data
log(f"loading cached {CITY} percolation curves (n={N_TARGET:,}, m={M_TARGET:,}) ...")
z = np.load(f"results/{TAG}_percolation.npz")
data = {CITY: {"q": z["q"],
                      "random": z["random"],
                      "degree (adaptive)": z["degree_(adaptive)"]}}

log("generating ER G(n, m) ...")
ip, ix, n = to_csr(nx.gnm_random_graph(N_TARGET, M_TARGET, seed=SEED))
data["ER"] = {"q": np.arange(n + 1) / n, **curves_for(ip, ix, n, "ER")}
del ip, ix; gc.collect()

log("generating WS (k=4, p=0.1) ...")
ip, ix, n = to_csr(nx.watts_strogatz_graph(N_TARGET, 4, 0.1, seed=SEED))
data["WS"] = {"q": np.arange(n + 1) / n, **curves_for(ip, ix, n, "WS")}
del ip, ix; gc.collect()

lines = ["--- percolation comparison vs models (q_c = first q with S < 1%) ---"]
for g in (CITY, "ER", "WS"):
    lines.append(f"  {g:9s}: random q_c = {qc(data[g]['q'], data[g]['random']):.3f}   "
                 f"adaptive-degree q_c = {qc(data[g]['q'], data[g]['degree (adaptive)']):.3f}")
report = "\n".join(lines) + "\n"
print(report)
with open(f"results/{TAG}_percolation_report.txt", "a") as f:
    f.write(report)

# ---------------------------------------------------------------- dark chart
import matplotlib.pyplot as plt

SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
COLORS = {CITY: "#3987e5", "ER": "#d95926", "WS": "#199e70"}   # dark-mode steps
MARKERS = {CITY: "o", "ER": "s", "WS": "^"}

fig, ax = plt.subplots(figsize=(8, 6), dpi=200)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
for g in (CITY, "ER", "WS"):
    q, c = data[g]["q"], data[g]["random"]
    step = max(1, len(q) // 500)
    ax.plot(q[::step], c[::step], color=COLORS[g], lw=1.8,
            marker=MARKERS[g], ms=4.5, markevery=30, label=g)
ax.set_title(f"Site percolation, random failure: {CITY} vs ER and WS\n(WS is ~52% denser: k=4)",
             color=INKD, fontsize=12)
ax.set_xlabel("fraction of nodes removed  q", color=SUBD)
ax.set_ylabel("giant component share  S(q)", color=SUBD)
ax.set_xlim(0, 0.75)
ax.set_ylim(0, 1.0)
ax.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.tick_params(colors=SUBD)
for s_ in ax.spines.values():
    s_.set_color(SPINE)
ax.legend(frameon=False, labelcolor=INKD, fontsize=10)
fig.tight_layout()
fig.savefig(f"figures/{TAG}_percolation_vs_models_dark.png",
            bbox_inches="tight", facecolor=SURFACE)
plt.close(fig)
log(f"done: figures/{TAG}_percolation_vs_models_dark.png")
