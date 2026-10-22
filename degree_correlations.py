"""
Degree correlations of the Bengaluru drive network (Newman, Networks, sec. 7.13):

  1. e_jk       — joint distribution of the degrees at the two ends of an edge
  2. k_nn(k)    — mean degree of a node's neighbours vs its own degree k
  3. map        — nodes coloured by local assortativity (Piraveenan et al. 2010,
                  contributions that sum exactly to Newman's r, eq. 7.82)

Important nodes = top betweenness intersections from
results/<tag>_top_betweenness.csv (analyze_bengaluru.py), numbered on every panel.

Usage:  python degree_correlations.py [path/to/graph.graphml]
Writes figures/bengaluru_degree_correlations_dark.png and
results/bengaluru_degree_correlations_report.txt
"""
import csv
import gc
import os
import sys
import time

import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
from matplotlib.collections import LineCollection

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "bengaluru_city_bengaluru_drive (10).graphml"
TAG = "bengaluru_city_bengaluru_drive"
N_LABEL = 8
SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
BLUE, ORANGE, BASE = "#3987e5", "#ff8a3d", "#3a3a37"

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------- graph
log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
S = nx.Graph()                                   # same simplification as model_comparison.py
S.add_nodes_from(G.nodes)
S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
xy_all = {n: (d["x"], d["y"]) for n, d in G.nodes(data=True)}
del G
gc.collect()
S = S.subgraph(max(nx.connected_components(S), key=len)).copy()
nodes = list(S.nodes)
idx = {n: i for i, n in enumerate(nodes)}
xy = np.array([xy_all[n] for n in nodes])
deg = np.array([S.degree(n) for n in nodes])
edges = np.array([(idx[u], idx[v]) for u, v in S.edges])
M = len(edges)
log(f"giant component: {len(nodes):,} nodes, {M:,} edges")

# ---------------------------------------------------------------- assortativity (eq. 7.82)
# r = sum_ij (A_ij - k_i k_j/2m) k_i k_j / sum_ij (k_i d_ij - k_i k_j/2m) k_i k_j
ki, kj = deg[edges[:, 0]].astype(float), deg[edges[:, 1]].astype(float)
two_m = 2.0 * M
num = 2 * (ki * kj).sum() - (deg.astype(float) ** 2).sum() ** 2 / two_m
den = (deg.astype(float) ** 3).sum() - (deg.astype(float) ** 2).sum() ** 2 / two_m
r = num / den
r_nx = nx.degree_assortativity_coefficient(S)
assert abs(r - r_nx) < 1e-9, (r, r_nx)
log(f"assortativity r = {r:.4f} (networkx agrees)")

# ---------------------------------------------------------------- e_jk
kmax = deg.max()
ejk = np.zeros((kmax + 1, kmax + 1))
np.add.at(ejk, (deg[edges[:, 0]], deg[edges[:, 1]]), 1)
np.add.at(ejk, (deg[edges[:, 1]], deg[edges[:, 0]]), 1)
ejk /= ejk.sum()

# ---------------------------------------------------------------- k_nn(k)
nbr_mean = np.array([np.mean([deg[idx[w]] for w in S[n]]) for n in nodes])
ks = np.arange(1, kmax + 1)
knn = np.array([nbr_mean[deg == k].mean() if (deg == k).any() else np.nan for k in ks])
knn_sd = np.array([nbr_mean[deg == k].std() if (deg == k).any() else np.nan for k in ks])
cnt = np.array([(deg == k).sum() for k in ks])
k1, k2 = deg.mean(), (deg.astype(float) ** 2).mean()
baseline = k2 / k1                                # k_nn(k) if degrees were uncorrelated
slope = np.polyfit(deg, nbr_mean, 1)[0]

# ---------------------------------------------------------------- local assortativity
# rho_v = q_v * sum_{w~v}(q_w - mu_q) / (2M sigma_q^2),  q = k-1 (excess degree); sum_v rho_v = r
q = deg - 1.0
qdist = q[edges].ravel()                           # excess degree at both ends of every edge
mu_q, var_q = qdist.mean(), qdist.var()
nbr_sum = np.zeros(len(nodes))
np.add.at(nbr_sum, edges[:, 0], q[edges[:, 1]] - mu_q)
np.add.at(nbr_sum, edges[:, 1], q[edges[:, 0]] - mu_q)
rho = q * nbr_sum / (two_m * var_q)
assert abs(rho.sum() - r) < 1e-6, (rho.sum(), r)

# ---------------------------------------------------------------- important nodes
important = []
seen = set()
with open(f"results/{TAG}_top_betweenness.csv", newline="") as f:
    for row in csv.DictReader(f):
        nd = int(row["osmid"])
        if nd not in idx or row["streets"] in seen:    # skip twin nodes of one junction
            continue
        seen.add(row["streets"])
        important.append((nd, row["streets"], float(row["betweenness"])))
        if len(important) == N_LABEL:
            break
imp_idx = [idx[nd] for nd, _, _ in important]

# ---------------------------------------------------------------- report
report = f"""================= Degree correlations: {TAG} =================
Nodes / edges (simple, undirected, giant component) : {len(nodes):,} / {M:,}
<k> = {k1:.3f}   <k^2>/<k> = {baseline:.3f}   k_max = {kmax}
Assortativity r (Newman eq. 7.82)    : {r:+.4f}
Slope of neighbour degree vs degree  : {slope:+.4f}
k_nn(k):  """ + "  ".join(f"k={k}: {v:.2f}" for k, v in zip(ks, knn) if not np.isnan(v)) + """

Important nodes (top betweenness) — degree, mean neighbour degree, local assortativity:
""" + "\n".join(
    f"  {i}. {name[:50]:50s} k={deg[ix]}  k_nn={nbr_mean[ix]:.2f}  rho={rho[ix]:+.2e}  bc={bc:.4f}"
    for i, ((nd, name, bc), ix) in enumerate(zip(important, imp_idx), 1)) + "\n"
print(report)
with open(f"results/bengaluru_degree_correlations_report.txt", "w") as f:
    f.write(report)

# ---------------------------------------------------------------- figure
log("plotting ...")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)


def tag_pts(ax, xs, ys):
    """Ring the important nodes; nodes sharing a spot get one combined label (e.g. '1,4,6')."""
    groups = {}
    for i, (x, y) in enumerate(zip(xs, ys), 1):
        groups.setdefault((round(float(x), 2), round(float(y), 2)), []).append(i)
    for (x, y), ids in groups.items():
        ax.scatter([x], [y], s=90, facecolors="none", edgecolors=ORANGE, linewidths=1.5, zorder=5)
        ax.annotate(",".join(map(str, ids)), (x, y), xytext=(8, 8), textcoords="offset points",
                    color=ORANGE, fontsize=9, fontweight="bold", zorder=6,
                    path_effects=[mpl.patheffects.withStroke(linewidth=2.5, foreground=SURFACE)])


fig = plt.figure(figsize=(19, 6.5), dpi=200)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], wspace=0.28)

# panel 1: e_jk
ax = fig.add_subplot(gs[0])
style(ax)
disp = np.ma.masked_less_equal(ejk[1:, 1:], 0)
im = ax.imshow(disp, origin="lower", extent=(0.5, kmax + 0.5, 0.5, kmax + 0.5),
               cmap=mpl.colors.LinearSegmentedColormap.from_list("dk", ["#232b36", BLUE, "#bcd9ff"]),
               norm=mpl.colors.LogNorm(vmin=disp.min(), vmax=disp.max()))
ax.set_xlabel("degree j at one end of an edge", color=SUBD)
ax.set_ylabel("degree k at the other end", color=SUBD)
ax.set_title("Joint degree distribution  $e_{jk}$", color=INKD, fontsize=12)
ax.set_xticks(ks)
ax.set_yticks(ks)
cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cb.set_label("fraction of edge ends (log)", color=SUBD)
cb.ax.tick_params(colors=SUBD)
tag_pts(ax, deg[imp_idx], np.round(nbr_mean[imp_idx]).clip(1, kmax))

# panel 2: k_nn(k)
ax = fig.add_subplot(gs[1])
style(ax)
ok = cnt >= 5                                      # drop sparsely populated degrees
ax.errorbar(ks[ok], knn[ok], yerr=knn_sd[ok], fmt="o-", color=BLUE, ecolor="#2e5f98",
            capsize=3, ms=6, label="Bengaluru  $k_{nn}(k)$  (± 1 sd)")
ax.axhline(baseline, color=SUBD, ls="--", lw=1.2, label=f"no correlation:  <k²>/<k> = {baseline:.2f}")
ax.set_xlabel("node degree  k", color=SUBD)
ax.set_ylabel("mean degree of neighbours", color=SUBD)
ax.set_title(f"Average neighbour degree      r = {r:+.3f}", color=INKD, fontsize=12)
ax.set_xticks(ks)
ax.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9, loc="lower right")
tag_pts(ax, deg[imp_idx], nbr_mean[imp_idx])

# panel 3: map coloured by local assortativity
ax = fig.add_subplot(gs[2])
style(ax)
segs = np.stack([xy[edges[:, 0]], xy[edges[:, 1]]], axis=1)
ax.add_collection(LineCollection(segs, colors=BASE, linewidths=0.25, zorder=1))
lim = np.percentile(np.abs(rho), 99.5)
order = np.argsort(np.abs(rho))                    # strongest contributors drawn on top
sc = ax.scatter(xy[order, 0], xy[order, 1], c=rho[order], s=1.2, zorder=2,
                cmap=mpl.colors.LinearSegmentedColormap.from_list(
                    "div", [BLUE, "#2b2b29", "#ff5a5a"]),
                norm=mpl.colors.TwoSlopeNorm(0, -lim, lim), linewidths=0)
ax.set_aspect(1 / np.cos(np.radians(xy[:, 1].mean())))
ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("Local assortativity  (blue: disassortative, red: assortative)", color=INKD, fontsize=12)
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.01)
cb.set_label("node contribution to r", color=SUBD)
cb.ax.tick_params(colors=SUBD, labelsize=7)
px, py = xy[imp_idx, 0], xy[imp_idx, 1]
ax.scatter(px, py, s=60, facecolors="none", edgecolors=ORANGE, linewidths=1.2, zorder=5)
padx, pady = 0.012, 0.010
x0, x1, y0, y1 = px.min() - padx, px.max() + padx, py.min() - pady, py.max() + pady
ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=ORANGE, lw=0.8, zorder=4))
ins = ax.inset_axes([0.0, 0.0, 0.52, 0.52])          # zoom on the central business district
ins.set_facecolor(SURFACE)
ins.add_collection(LineCollection(segs, colors="#55554f", linewidths=0.5, zorder=1))
inb = (xy[:, 0] > x0) & (xy[:, 0] < x1) & (xy[:, 1] > y0) & (xy[:, 1] < y1)
ins.scatter(xy[inb, 0], xy[inb, 1], c=rho[inb], s=6, zorder=2, linewidths=0,
            cmap=sc.get_cmap(), norm=sc.norm)
ins.scatter(px, py, s=110, facecolors="none", edgecolors=ORANGE, linewidths=1.6, zorder=5)
for i, (x, y) in enumerate(zip(px, py), 1):
    ins.annotate(str(i), (x, y), xytext=(7, 7), textcoords="offset points", color=ORANGE,
                 fontsize=10, fontweight="bold", zorder=6,
                 path_effects=[mpl.patheffects.withStroke(linewidth=2.5, foreground=SURFACE)])
ins.set_xlim(x0, x1)
ins.set_ylim(y0, y1)
ins.set_aspect(1 / np.cos(np.radians(xy[:, 1].mean())))
ins.set_xticks([])
ins.set_yticks([])
for s_ in ins.spines.values():
    s_.set_color(ORANGE)

legend = "Important nodes (top betweenness)\n" + "\n".join(
    f"{i}. {name[:34]}   (k={deg[ix]})" for i, ((_, name, _), ix) in enumerate(zip(important, imp_idx), 1))
fig.text(0.015, 0.015, legend, color=INKD, fontsize=9, va="bottom", family="monospace")
fig.suptitle("Degree correlations — Bengaluru road network", color=INKD, fontsize=14, y=0.99)
fig.subplots_adjust(bottom=0.30, top=0.9, left=0.05, right=0.98)
fig.savefig("figures/bengaluru_degree_correlations_dark.png", facecolor=SURFACE)
log("done: figures/bengaluru_degree_correlations_dark.png, "
    "results/bengaluru_degree_correlations_report.txt")
