"""
Compare the Bengaluru street network against four canonical random-graph models,
each generated at the same number of nodes (matching edges where the model's
integer parameters allow):

  ER : Erdos-Renyi G(n, m)            — same n AND m exactly
  WS : Watts-Strogatz (k=4, p=0.1)    — k must be an even integer, so m ~ 2n (52% denser)

Metrics: degree distribution (CCDF), max degree, assortativity, clustering
(approx), mean shortest path + diameter lower bound (sampled BFS), spectral gap
lambda_2 (expansion), bridges share, cyclomatic number / alpha index.

Usage:  python model_comparison.py [path/to/graph.graphml]
Writes results/model_comparison.csv, results/model_degrees.npz, 2 figures.
"""
import csv
import gc
import os
import sys
import time

import networkx as nx
import numpy as np
import osmnx as ox
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from networkx.algorithms.approximation import average_clustering as approx_clustering
from scipy.sparse import csgraph

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
SEED = 42
N_BFS = 48            # sampled BFS sources for path lengths
CL_TRIALS = 20000     # triangle-sampling trials for approximate clustering
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

COLORS = {"Bengaluru": "#2a78d6", "ER": "#eb6834", "WS": "#1baf7a",
          "AB": "#eda100", "HK": "#e87ba4"}
MARKERS = {"Bengaluru": "o", "ER": "s", "WS": "^", "AB": "D", "HK": "v"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def lambda2(G, planar_like):
    """Second-smallest eigenvalue of the normalized Laplacian."""
    A = nx.to_scipy_sparse_array(G, weight=None, format="csr", dtype=np.float64)
    L = csgraph.laplacian(A, normed=True)
    try:
        if planar_like:   # good separators -> cheap LU factorization, shift-invert
            vals = spla.eigsh(L.tocsc(), k=2, sigma=-1e-2, which="LM",
                              return_eigenvectors=False)
            return float(np.sort(vals)[1])
        # expanders: Lanczos on 2I - L for the top of the spectrum (no factorization)
        M = 2.0 * sp.identity(L.shape[0], format="csr") - L.tocsr()
        vals = spla.eigsh(M, k=2, which="LA", return_eigenvectors=False,
                          maxiter=8000, tol=1e-7)
        return float(2.0 - np.sort(vals)[0])
    except Exception as e:
        log(f"  lambda2 failed ({type(e).__name__}); retrying with shift-invert")
        try:
            vals = spla.eigsh(L.tocsc(), k=2, sigma=-1e-2, which="LM",
                              return_eigenvectors=False)
            return float(np.sort(vals)[1])
        except Exception:
            return float("nan")


def measure(G, name, planar_like=False):
    """G: simple undirected Graph. Returns metric dict + degree histogram."""
    t0 = time.time()
    G.remove_edges_from(nx.selfloop_edges(G))
    n_all = G.number_of_nodes()
    if not nx.is_connected(G):
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    n, m = G.number_of_nodes(), G.number_of_edges()
    degs = np.array([d for _, d in G.degree()])
    rng = np.random.default_rng(SEED)

    log(f"[{name}] GCC n={n:,} ({n/n_all:.1%}), m={m:,}; paths ...")
    nodes = list(G.nodes)
    tot = cnt = diam = 0
    for s in rng.choice(len(nodes), size=min(N_BFS, n), replace=False):
        dist = nx.single_source_shortest_path_length(G, nodes[s])
        v = np.fromiter(dist.values(), dtype=np.int64)
        tot += v.sum(); cnt += len(v) - 1; diam = max(diam, v.max())
    log(f"[{name}] clustering / assortativity / bridges ...")
    C = approx_clustering(G, trials=CL_TRIALS, seed=SEED)
    r = nx.degree_assortativity_coefficient(G)
    n_bridges = sum(1 for _ in nx.bridges(G))
    log(f"[{name}] lambda_2 ...")
    l2 = lambda2(G, planar_like)

    mu = m - n + 1
    row = {
        "graph": name, "n_gcc": n, "m_gcc": m, "gcc_frac": round(n / n_all, 4),
        "d_avg": round(2 * m / n, 3), "d_max": int(degs.max()),
        "assortativity": round(r, 4), "clustering_approx": round(C, 4),
        "mean_path_hops": round(tot / cnt, 2), "diameter_lb": int(diam),
        "lambda2": f"{l2:.3e}", "bridge_frac": round(n_bridges / m, 4),
        "cyclomatic_mu": mu, "alpha_index": round(mu / (2 * n - 5), 4),
    }
    log(f"[{name}] done in {time.time()-t0:.0f}s: {row}")
    return row, np.bincount(degs)


rows, deg_hists = [], {}

def record(row, hist, name):
    rows.append(row)
    deg_hists[name] = hist
    with open("results/model_comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    np.savez_compressed("results/model_degrees.npz",
                        **{k: v for k, v in deg_hists.items()})


# ---------------------------------------------------------------- Bengaluru
log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
S = nx.Graph()
S.add_nodes_from(G.nodes)
S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
del G; gc.collect()
n_target, m_target = S.number_of_nodes(), S.number_of_edges()
row, hist = measure(S, "Bengaluru", planar_like=True)
record(row, hist, "Bengaluru")
del S; gc.collect()

# ---------------------------------------------------------------- models
log(f"targets: n={n_target:,}, m={m_target:,} (d_avg={2*m_target/n_target:.3f})")

log("generating ER: G(n, m) ...")
M = nx.gnm_random_graph(n_target, m_target, seed=SEED)
record(*measure(M, "ER"), "ER")
del M; gc.collect()

log("generating WS: watts_strogatz(k=4, p=0.1)  [k integer -> ~52% denser] ...")
M = nx.watts_strogatz_graph(n_target, 4, 0.1, seed=SEED)
record(*measure(M, "WS"), "WS")
del M; gc.collect()

# ---------------------------------------------------------------- figures
import matplotlib.pyplot as plt

ORDER = ["Bengaluru", "ER", "WS"]
INK, SUB = "#1a1a19", "#5f5e56"

# degree distribution (grouped bars, counts + normalized) is rendered by
# degree_distribution.py from the cached results/model_degrees.npz

log("figure: metric panels ...")
PANELS = [
    ("clustering_approx", "avg clustering C", None),
    ("mean_path_hops", "mean shortest path (hops)", None),
    ("diameter_lb", "diameter (lower bound, hops)", None),
    ("lambda2", "spectral gap λ₂  (log scale)", "log"),
    ("assortativity", "degree assortativity r", None),
    ("bridge_frac", "fraction of edges that are bridges", None),
]
by = {r["graph"]: r for r in rows}
fig, axes = plt.subplots(2, 3, figsize=(14, 8), dpi=200)
for ax, (key, title, scale) in zip(axes.ravel(), PANELS):
    vals = [float(by[name][key]) for name in ORDER]
    ax.bar(range(len(ORDER)), vals, color=[COLORS[o] for o in ORDER], width=0.62)
    for i, v in enumerate(vals):
        txt = f"{v:.1e}" if scale == "log" else (f"{v:,.2f}" if abs(v) < 100 else f"{v:,.0f}")
        ax.text(i, v, f" {txt}", ha="center",
                va="bottom" if v >= 0 else "top", fontsize=8, color=INK)
    if scale == "log":
        ax.set_yscale("log")
    if min(vals) < 0:
        ax.axhline(0, color="#d4d3c9", lw=0.8)
    ax.set_xticks(range(len(ORDER)), ORDER, color=INK)
    ax.set_title(title, color=INK, fontsize=11)
    ax.grid(True, axis="y", color="#e8e7e0", lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=SUB)
    for s_ in ax.spines.values():
        s_.set_color("#d4d3c9")
fig.suptitle("Bengaluru street network vs canonical models (same n; m matched where integer params allow)",
             color=INK, fontsize=13)
fig.tight_layout()
fig.savefig("figures/model_metric_panels.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

log("done: results/model_comparison.csv, results/model_degrees.npz, "
    "figures/model_degree_ccdf.png, figures/model_metric_panels.png")
