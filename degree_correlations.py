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
from collections import Counter
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

# road attributes per node, from all incident edges (OSM highway tag, bridge, name, ref)
CLASS_RANK = {"motorway": 5, "trunk": 5, "primary": 4, "secondary": 3, "tertiary": 2}
CLASS_NAMES = {5: "Expressway / trunk", 4: "Primary", 3: "Secondary", 2: "Tertiary", 1: "Local street"}
node_rank, node_flags, node_names = {}, {}, {}


def as_list(v):
    return v if isinstance(v, list) else [v]


for u, v, d in G.edges(data=True):
    rk = max(CLASS_RANK.get(h.replace("_link", ""), 1) for h in as_list(d["highway"]))
    names = " ".join(str(x) for x in as_list(d.get("name", "")))
    fl = set()
    if d.get("bridge") not in (None, "no"):
        fl.add("Flyover / viaduct")
    if "Ring Road" in names:
        fl.add("Ring road")
    if "Outer Ring Road" in names:
        fl.add("Outer Ring Road")
    if any(str(x).startswith(("NH", "SH")) for x in as_list(d.get("ref", ""))):
        fl.add("National / state highway")
    for n in (u, v):
        node_rank[n] = max(node_rank.get(n, 1), rk)
        for nm in as_list(d.get('name', [])):
            if nm:
                node_names.setdefault(n, Counter())[nm] += 1
        node_flags.setdefault(n, set()).update(fl)
del G
gc.collect()
S = S.subgraph(max(nx.connected_components(S), key=len)).copy()
nodes = list(S.nodes)
idx = {n: i for i, n in enumerate(nodes)}
xy = np.array([xy_all[n] for n in nodes])
deg = np.array([S.degree(n) for n in nodes])
rank = np.array([node_rank[n] for n in nodes])
flags = [node_flags.get(n, set()) for n in nodes]
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


fig = plt.figure(figsize=(13, 6.5), dpi=200)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(1, 2, wspace=0.3)

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

legend = "Important nodes (top betweenness)\n" + "\n".join(
    f"{i}. {name[:34]}   (k={deg[ix]})" for i, ((_, name, _), ix) in enumerate(zip(important, imp_idx), 1))
fig.text(0.04, 0.015, legend, color=INKD, fontsize=9, va="bottom", family="monospace")
fig.suptitle("Degree correlations — Bengaluru road network (important nodes numbered)", color=INKD, fontsize=14, y=0.99)
fig.subplots_adjust(bottom=0.30, top=0.9, left=0.06, right=0.98)
fig.savefig("figures/bengaluru_degree_correlations_dark.png", facecolor=SURFACE)
plt.close(fig)


# ================================================================ local assortativity maps
DIV = mpl.colors.LinearSegmentedColormap.from_list("div", [BLUE, "#2b2b29", "#ff5a5a"])
lim = np.percentile(np.abs(rho), 99.5)
DNORM = mpl.colors.TwoSlopeNorm(0, -lim, lim)
segs = np.stack([xy[edges[:, 0]], xy[edges[:, 1]]], axis=1)
ASPECT = 1 / np.cos(np.radians(xy[:, 1].mean()))
px, py = xy[imp_idx, 0], xy[imp_idx, 1]
halo = [mpl.patheffects.withStroke(linewidth=2.5, foreground=SURFACE)]


def rho_map(ax, mask, size, base_lw, base_color, number_size):
    ax.set_facecolor(SURFACE)
    ax.add_collection(LineCollection(segs, colors=base_color, linewidths=base_lw, zorder=1))
    sel = np.where(mask)[0]
    sel = sel[np.argsort(np.abs(rho[sel]))]            # strongest contributors on top
    sc = ax.scatter(xy[sel, 0], xy[sel, 1], c=rho[sel], s=size, zorder=2, cmap=DIV, norm=DNORM,
                    linewidths=0)
    ax.set_aspect(ASPECT)
    ax.set_xticks([])
    ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    for i, (x, y) in enumerate(zip(px, py), 1):
        ax.scatter([x], [y], s=number_size * 8, facecolors="none", edgecolors=ORANGE,
                   linewidths=1.3, zorder=5)
        ax.annotate(str(i), (x, y), xytext=(7, 7), textcoords="offset points", color=ORANGE,
                    fontsize=number_size, fontweight="bold", zorder=6, path_effects=halo)
    return sc


def add_cb(fig, ax, sc):
    cb = fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.01)
    cb.set_label("node contribution to r  (sums to r)", color=SUBD)
    cb.ax.tick_params(colors=SUBD, labelsize=8)
    cb.formatter.set_powerlimits((0, 0))
    cb.update_ticks()
    cb.ax.yaxis.offsetText.set_color(SUBD)


# full-city map, nothing covering it
log("plotting local assortativity map ...")
fig, ax = plt.subplots(figsize=(12, 11), dpi=200)
fig.patch.set_facecolor(SURFACE)
sc = rho_map(ax, np.ones(len(nodes), bool), 1.6, 0.2, "#33332f", 10)
ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
add_cb(fig, ax, sc)
ax.set_title("Local assortativity of Bengaluru intersections\n"
             "red: links to higher-degree nodes than average (assortative) · blue: disassortative",
             color=INKD, fontsize=12)
fig.savefig("figures/bengaluru_local_assortativity_map_dark.png", bbox_inches="tight", facecolor=SURFACE)
plt.close(fig)

# zoom on the central business district, separate figure
fig, ax = plt.subplots(figsize=(10, 9), dpi=200)
fig.patch.set_facecolor(SURFACE)
x0, x1 = px.min() - 0.012, px.max() + 0.012
y0, y1 = py.min() - 0.010, py.max() + 0.010
inb = (xy[:, 0] > x0) & (xy[:, 0] < x1) & (xy[:, 1] > y0) & (xy[:, 1] < y1)
sc = rho_map(ax, inb, 9, 0.5, "#55554f", 12)
ax.set_xlim(x0, x1)
ax.set_ylim(y0, y1)
add_cb(fig, ax, sc)
ax.set_title("Local assortativity — central business district", color=INKD, fontsize=12)
fig.savefig("figures/bengaluru_local_assortativity_cbd_dark.png", bbox_inches="tight", facecolor=SURFACE)
plt.close(fig)

# ================================================================ by road category
log("categorising nodes by road type ...")
cbd = np.array([px.mean(), py.mean()])                  # centroid of the top-betweenness junctions
dist_km = np.hypot((xy[:, 0] - cbd[0]) * np.cos(np.radians(cbd[1])), xy[:, 1] - cbd[1]) * 111.2
cats = {}                                               # name -> boolean mask
for r_ in (5, 4, 3, 2, 1):
    cats[CLASS_NAMES[r_]] = rank == r_                  # exclusive: highest class at the node
for fl in ("Flyover / viaduct", "Ring road", "Outer Ring Road", "National / state highway"):
    cats[fl] = np.array([fl in f_ for f_ in flags])     # overlapping flags
cats["Central (< 3 km of CBD)"] = dist_km < 3
cats["Outskirts (> 12 km)"] = dist_km > 12
cats["Important nodes (top betweenness)"] = np.isin(np.arange(len(nodes)), imp_idx)
rho_n = rho * len(nodes)                                # rescale so the mean over all nodes = r
rng = np.random.default_rng(0)
rows = []
for name, m_ in cats.items():
    v = rho_n[m_]
    boot = [rng.choice(v, len(v)).mean() for _ in range(300)] if len(v) > 1 else [v.mean()]
    rows.append((name, m_.sum(), v.mean(), np.percentile(boot, 2.5), np.percentile(boot, 97.5),
                 deg[m_].mean(), (rho_n < v.mean()).mean()))

cat_report = f"\nLocal assortativity by road category (rho x N, so the mean over all nodes = r = {r:+.4f}):\n"
cat_report += f"{'category':36s} {'nodes':>8s} {'mean':>9s} {'95% CI':>20s} {'<k>':>6s}\n"
for name, n_, mu, lo, hi, kmean, _ in rows:
    cat_report += f"{name:36s} {n_:8,d} {mu:+9.4f} [{lo:+8.4f},{hi:+8.4f}] {kmean:6.2f}\n"
print(cat_report)
with open("results/bengaluru_degree_correlations_report.txt", "a", encoding="utf-8") as f:
    f.write(cat_report)

fig = plt.figure(figsize=(19, 8.5), dpi=200)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.08)
ax = fig.add_subplot(gs[0])
ax.set_facecolor(SURFACE)
ypos = np.arange(len(rows))[::-1]
for y_, (name, n_, mu, lo, hi, kmean, _) in zip(ypos, rows):
    ax.errorbar(mu, y_, xerr=[[mu - lo], [hi - mu]], fmt="o", color=ORANGE if mu > r else BLUE,
                ecolor=SUBD, capsize=4, ms=8)
ax.axvline(r, color=SUBD, ls="--", lw=1.2)
ax.text(r, -0.75, f" r = {r:+.3f} (whole city)", color=SUBD, fontsize=9, va="bottom")
ax.set_yticks(ypos)
ax.set_yticklabels([f"{r_[0]}\n(n={r_[1]:,}, <k>={r_[5]:.2f})" for r_ in rows], color=INKD, fontsize=9)
ax.axhline(ypos[4] - 0.5, color=SPINE, lw=0.8)         # class | flags
ax.axhline(ypos[8] - 0.5, color=SPINE, lw=0.8)         # flags | location
ax.set_xlabel("mean local assortativity  ρ·N   (95% bootstrap CI)", color=SUBD)
ax.set_title("Where each road type sits  (orange: above city average, blue: below)", color=INKD, fontsize=12)
ax.tick_params(colors=SUBD)
ax.grid(True, axis="x", color=GRID, lw=0.6)
ax.set_axisbelow(True)
for s_ in ax.spines.values():
    s_.set_color(SPINE)

ax = fig.add_subplot(gs[1])
ax.set_facecolor(SURFACE)
ax.add_collection(LineCollection(segs, colors="#2a2a28", linewidths=0.15, zorder=1))
LAYERS = [("Local street", "#6e6e68", 0.3, "local streets"),
          ("Primary", "#f2c14e", 1.2, "primary"),
          ("Expressway / trunk", "#ff5a5a", 3.0, "expressway / trunk"),
          ("Flyover / viaduct", "#3ddc97", 3.0, "flyovers"),
          ("Outer Ring Road", "#c77dff", 3.0, "Outer Ring Road")]
handles = []
for name, col, size, lab in LAYERS:
    m_ = cats[name]
    ax.scatter(xy[m_, 0], xy[m_, 1], s=size, color=col, linewidths=0, zorder=2, alpha=0.9 if size > 1 else 0.35)
    handles.append(plt.Line2D([], [], marker="o", ls="", color=col, ms=7, label=f"{lab}  ({m_.sum():,})"))
ax.legend(handles=handles, loc="upper right", frameon=False, labelcolor=INKD, fontsize=9)
ax.set_aspect(ASPECT)
ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("Where those road types are in the city", color=INKD, fontsize=12)
for s_ in ax.spines.values():
    s_.set_color(SPINE)
fig.savefig("figures/bengaluru_assortativity_by_road_type_dark.png", bbox_inches="tight", facecolor=SURFACE)
log("done: figures/bengaluru_degree_correlations_dark.png, _local_assortativity_map_dark.png, "
    "_local_assortativity_cbd_dark.png, _assortativity_by_road_type_dark.png")


# ================================================================ extremes: disassortative + outlier clusters
from collections import Counter

from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

log("extremes: disassortative nodes/edges, boundary clusters ...")


def top_name(i, n=2):
    c = node_names.get(nodes[i])
    return "; ".join(nm for nm, _ in c.most_common(n)) if c else "(unnamed)"


def compass(dx, dy):
    return ["E", "NE", "N", "NW", "W", "SW", "S", "SE"][int(((np.degrees(np.arctan2(dy, dx)) + 22.5) % 360) // 45)]


ext = "\n================ Most DISassortative nodes (lowest rho*N) ================\n"
ext += "(rho_v = q_v * sum_nbrs(q_w - mu) / (2M sigma^2); negative = hub-like node whose neighbours are low-degree)\n"
for j in np.argsort(rho_n)[:15]:
    nb = sorted(int(deg[idx[w]]) for w in S[nodes[j]])
    ext += (f"  rho*N={rho_n[j]:+.3f}  k={deg[j]}  nbr degrees={nb}  "
            f"({xy[j, 1]:.4f}, {xy[j, 0]:.4f})  {dist_km[j]:4.1f} km from CBD  {top_name(j)[:55]}\n")
low = rho_n < -0.5
ext += (f"\nNodes with rho*N < -0.5: {low.sum():,} ({low.mean():.1%});  of these, "
        f"{np.mean([any(deg[idx[w]] == 1 for w in S[nodes[j]]) for j in np.where(low)[0]]):.0%} touch a dead end (k=1), "
        f"degree mix {dict(sorted(Counter(deg[low].tolist()).items()))}\n")
ext += f"Mean distance from CBD: disassortative {dist_km[low].mean():.1f} km vs all nodes {dist_km.mean():.1f} km\n"

# edge-level contribution: c_e = (q_i - mu)(q_j - mu) / (M sigma^2); sum_e c_e = r
ce = (q[edges[:, 0]] - mu_q) * (q[edges[:, 1]] - mu_q) / (M * var_q)
assert abs(ce.sum() - r) < 1e-6, (ce.sum(), r)
ce_n = ce * M                                      # mean over edges = r
ext += "\n================ Most DISassortative edges (hub linked to low-degree node) ================\n"
for e in np.argsort(ce_n)[:15]:
    a_, b_ = edges[e]
    mid = (xy[a_] + xy[b_]) / 2
    ext += (f"  c*M={ce_n[e]:+.3f}  degrees {deg[a_]}-{deg[b_]}  ({mid[1]:.4f}, {mid[0]:.4f})  "
            f"{np.hypot((mid[0]-cbd[0])*np.cos(np.radians(cbd[1])), mid[1]-cbd[1])*111.2:4.1f} km  "
            f"{top_name(a_, 1)[:30]} | {top_name(b_, 1)[:30]}\n")
pair = Counter()
for e in range(M):
    pair[tuple(sorted((int(deg[edges[e, 0]]), int(deg[edges[e, 1]]))))] += ce[e]
ext += "\nTotal contribution to r by edge type (deg_a, deg_b), most negative / most positive:\n"
srt = sorted(pair.items(), key=lambda kv: kv[1])
ext += "  neg: " + "  ".join(f"{k}:{v:+.4f}" for k, v in srt[:6]) + "\n"
ext += "  pos: " + "  ".join(f"{k}:{v:+.4f}" for k, v in srt[::-1][:6]) + "\n"

# radial trend (the "centre is higher" pattern)
bins = [0, 2, 4, 6, 8, 10, 12, 15, 20, 40]
radial = []
ext += "\nRadial profile (distance from CBD): mean rho*N, share of nodes with k>=3, mean k\n"
for lo_, hi_ in zip(bins[:-1], bins[1:]):
    m_ = (dist_km >= lo_) & (dist_km < hi_)
    if m_.sum():
        radial.append(((lo_ + hi_) / 2, rho_n[m_].mean(), rho_n[m_].std() / np.sqrt(m_.sum())))
        ext += f"  {lo_:2d}-{hi_:2d} km: n={m_.sum():6,d}  rho*N={rho_n[m_].mean():+.3f}  k>=3: {(deg[m_] >= 3).mean():.0%}  <k>={deg[m_].mean():.2f}\n"

# boundary clusters of high assortativity
OUTER_KM, LINK_KM = 10.0, 0.15
thr = np.percentile(rho_n, 90)
hi_nodes = np.where((rho_n >= thr) & (dist_km > OUTER_KM))[0]
X = (xy - cbd) * np.array([111.2 * np.cos(np.radians(cbd[1])), 111.2])
tree = cKDTree(X[hi_nodes])
pairs = tree.query_pairs(LINK_KM, output_type="ndarray")
g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(hi_nodes),) * 2)
_, lab = connected_components(g, directed=False)
clusters = [hi_nodes[lab == c] for c in range(lab.max() + 1) if (lab == c).sum() >= 40]
clusters.sort(key=len, reverse=True)
outer = dist_km > OUTER_KM
ext += (f"\n================ Boundary clusters of high assortativity ================\n"
        f"(nodes > {OUTER_KM:.0f} km from CBD with rho*N >= {thr:.2f} (city 90th pct), linked within {LINK_KM*1000:.0f} m, >= 40 nodes)\n"
        f"Outskirts baseline: k>=3 share {(deg[outer] >= 3).mean():.0%}, k=4 share {(deg[outer] == 4).mean():.0%}, "
        f"mean rho*N {rho_n[outer].mean():+.3f}\n")
cl_rows = []
for ci, cl in enumerate(clusters[:12], 1):
    c_xy = X[cl].mean(0)
    hw = Counter(rank[cl])
    names = Counter()
    for j in cl:
        for nm, cnt_ in node_names.get(nodes[j], Counter()).items():
            names[nm] += cnt_
    ll = xy[cl].mean(0)
    # grid-likeness: share of cluster nodes whose neighbours are all junctions (k>=3)
    allj = np.mean([all(deg[idx[w]] >= 3 for w in S[nodes[j]]) for j in cl])
    cl_rows.append((ci, ll, len(cl)))
    ext += (f"  #{ci:2d} {len(cl):5,d} nodes  {np.hypot(*c_xy):4.1f} km {compass(*c_xy)}  ({ll[1]:.4f}, {ll[0]:.4f})  "
            f"<k>={deg[cl].mean():.2f}  k=4: {(deg[cl] == 4).mean():.0%}  all-nbrs-junctions: {allj:.0%}  "
            f"local-class: {hw.get(1, 0) / len(cl):.0%}  streets: {', '.join(n for n, _ in names.most_common(3)) or '(unnamed)'}\n")
print(ext)
with open("results/bengaluru_degree_correlations_report.txt", "a", encoding="utf-8") as f:
    f.write(ext)

# figure
fig = plt.figure(figsize=(20, 9), dpi=200)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(1, 3, width_ratios=[0.7, 1.15, 1.15], wspace=0.12)

ax = fig.add_subplot(gs[0])
style(ax)
rx, ry, re = zip(*radial)
ax.errorbar(rx, ry, yerr=re, fmt="o-", color=BLUE, ecolor=SUBD, capsize=3, ms=6)
ax.axhline(r, color=SUBD, ls="--", lw=1)
ax.set_xlabel("distance from centre of CBD (km)", color=SUBD)
ax.set_ylabel("mean local assortativity  ρ·N", color=SUBD)
ax.set_title("Centre vs boundary", color=INKD, fontsize=12)
ax.grid(True, color=GRID, lw=0.6)

ax = fig.add_subplot(gs[1])
ax.set_facecolor(SURFACE)
ax.add_collection(LineCollection(segs, colors="#2f2f2c", linewidths=0.2, zorder=1))
ax.scatter(xy[:, 0], xy[:, 1], s=0.5, c="#4a4a45", linewidths=0, zorder=2)
for ci, cl in enumerate(clusters[:12], 1):
    ax.scatter(xy[cl, 0], xy[cl, 1], s=4, color="#ff5a5a", linewidths=0, zorder=3)
    cx, cy = xy[cl].mean(0)
    ax.annotate(str(ci), (cx, cy), xytext=(6, 6), textcoords="offset points", color="#ffd0d0",
                fontsize=11, fontweight="bold", zorder=6, path_effects=halo)
th = plt.Circle((cbd[0], cbd[1]), OUTER_KM / (111.2 * np.cos(np.radians(cbd[1]))), fill=False,
                ec=ORANGE, ls="--", lw=1)
ax.add_patch(th)
ax.set_aspect(ASPECT)
ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
ax.set_xticks([])
ax.set_yticks([])
ax.set_title(f"High-assortativity clusters beyond {OUTER_KM:.0f} km (dashed ring)", color=INKD, fontsize=12)
for s_ in ax.spines.values():
    s_.set_color(SPINE)

ax = fig.add_subplot(gs[2])
ax.set_facecolor(SURFACE)
ax.add_collection(LineCollection(segs, colors="#2f2f2c", linewidths=0.2, zorder=1))
ax.scatter(xy[:, 0], xy[:, 1], s=0.5, c="#4a4a45", linewidths=0, zorder=2)
worst_e = np.argsort(ce_n)[:300]
ax.add_collection(LineCollection(segs[worst_e], colors=BLUE, linewidths=2.2, zorder=3))
worst_n = np.argsort(rho_n)[:300]
ax.scatter(xy[worst_n, 0], xy[worst_n, 1], s=10, color="#6cb0ff", linewidths=0, zorder=4)
ax.set_aspect(ASPECT)
ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("Most disassortative: 300 edges (lines) and 300 nodes (dots)", color=INKD, fontsize=12)
for s_ in ax.spines.values():
    s_.set_color(SPINE)
fig.savefig("figures/bengaluru_assortativity_extremes_dark.png", bbox_inches="tight", facecolor=SURFACE)
log("done: figures/bengaluru_assortativity_extremes_dark.png")
