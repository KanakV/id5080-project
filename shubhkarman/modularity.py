"""
Python source export of `modularity.ipynb`.
The notebook is the primary, sectioned version of this analysis.
"""

# %%
import osmnx as ox
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
from collections import defaultdict
from tqdm.auto import tqdm

# %%
cities = ["chennai", "bengaluru"]
city = cities[1]
"""
# Download drivable road network of Chennai
G = ox.graph_from_place(
    "Chennai, Tamil Nadu, India",
    network_type="drive",
    simplify=True
)"""

""" 
# # 1. Load your pre-saved Chennai GraphML file (if not already loaded)
filepath = "chennai_roads.graphml"

# Load the graph back into OSMnx
G_whole = ox.load_graphml(filepath)

# 2. Define Tambaram coordinates
tambaram_lat = 12.9249
tambaram_lon = 80.1274

# 3. Find the nearest graph node to Tambaram
center_node = ox.distance.nearest_nodes(G_whole, X=tambaram_lon, Y=tambaram_lat)

# 4. Extract subgraph within a 10 km network distance (or straight-line distance)
#    Depth to stop the search. Only target nodes where the shortest path to
#    this node from the source node contains <= ``cutoff + 1`` nodes will be
#    included in the returned results.
radius = 50

# Get all nodes within 10 km network distance from Tambaram node
subgraph_nodes = nx.single_source_shortest_path_length(
    G_whole, center_node, cutoff=radius
).keys()

# Create the subgraph
G = G_whole.subgraph(subgraph_nodes).copy()

print("Tamabaram")
print(f"Nodes: {len(G.nodes):,}")
print(f"Edges: {len(G.edges):,}\n")

print("Whole Chenni")
print(f"Nodes: {len(G_whole.nodes):,}")
print(f"Edges: {len(G_whole.edges):,}")
# Optional: save for reuse
ox.save_graphml(G, "chennai_tambaram_roads.graphml")"""

filepath = f"{city}_roads.graphml"#"chennai_tambaram_roads.graphml"

G = ox.load_graphml(filepath)
# Plot
print(f"Nodes: {len(G.nodes):,}")
print(f"Edges: {len(G.edges):,}")
fig, ax = ox.plot_graph(
    G,
    node_size=0,
    edge_color="black",
    edge_linewidth=0.4,
    bgcolor="white"
)

# %%
# ---------- CELL 1: modularity (Leiden) with a max module size ----------
import numpy as np, osmnx as ox, igraph as ig
import matplotlib.pyplot as plt

# ---- config ----
MAX_SIZE = 5000          # hard cap on nodes per module
RESOLUTION = 1.0         # higher -> more, smaller modules on the first pass
EDGE_WEIGHT = "unit"     # "unit" = every road link equal, "inv_time" = faster roads bind tighter
SEED = 42
np.random.seed(SEED)

# ---- prep graph (same edge order as before) ----
G = ox.truncate.largest_component(G, strongly=True)
G = ox.add_edge_speeds(G)
G = ox.add_edge_travel_times(G)

nodes = list(G.nodes)
idx = {n: i for i, n in enumerate(nodes)}
N = len(nodes)
edges = list(G.edges(keys=True, data=True))       # list position == edge id
es = np.array([idx[u] for u, v, k, d in edges])
ed = np.array([idx[v] for u, v, k, d in edges])
w = np.array([d["travel_time"] for *_, d in edges], dtype=float)
xs = np.array([G.nodes[n]["x"] for n in nodes])
ys = np.array([G.nodes[n]["y"] for n in nodes])

# ---- undirected weighted graph for community detection ----
gu = ig.Graph(n=N, edges=list(zip(es.tolist(), ed.tolist())), directed=False)
gu.es["cw"] = np.ones(len(es)) if EDGE_WEIGHT == "unit" else (1.0 / w).tolist()
gu.simplify(multiple=True, loops=True, combine_edges="max")

def leiden(sub, res):
    return np.array(sub.community_leiden(
        objective_function="modularity", weights="cw",
        resolution=res, n_iterations=-1).membership)

def split_module(vids, res):
    """Recursively split a module until every piece is <= MAX_SIZE."""
    vids = np.sort(np.asarray(vids))
    if len(vids) <= MAX_SIZE:
        return [vids]
    sub = gu.induced_subgraph(vids.tolist())
    r, memb = res, None
    for _ in range(10):                            # raise resolution until it actually splits
        memb = leiden(sub, r)
        if memb.max() > 0:
            break
        r *= 2
    if memb.max() == 0:                            # fallback: spatial bisection
        a = xs if np.ptp(xs[vids]) >= np.ptp(ys[vids]) else ys
        memb = (a[vids] > np.median(a[vids])).astype(int)
        if memb.min() == memb.max():
            memb[: len(memb) // 2] = 1
    out = []
    for c in np.unique(memb):
        out += split_module(vids[memb == c], res * 1.5)
    return out

# ---- first pass on the whole city, then enforce the size cap ----
memb0 = leiden(gu, RESOLUTION)
print(f"First pass: {memb0.max()+1} modules, modularity Q = "
      f"{gu.modularity(memb0.tolist(), weights='cw'):.3f}")

pieces = []
for c in np.unique(memb0):
    pieces += split_module(np.where(memb0 == c)[0], RESOLUTION)

pieces.sort(key=len, reverse=True)
module = np.empty(N, dtype=int)
for m, vids in enumerate(pieces):
    module[vids] = m

sizes = np.array([len(p) for p in pieces])
intra = (module[es] == module[ed])
print(f"Final: {len(pieces)} modules (max size {sizes.max()}, median {int(np.median(sizes))}, min {sizes.min()})")
print(f"Final modularity Q = {gu.modularity(module.tolist(), weights='cw'):.3f}")
print(f"Edges inside modules: {intra.mean():.1%}   between modules: {(~intra).mean():.1%}")

# ---- plot the modules ----
gdf = ox.convert.graph_to_gdfs(G, nodes=False)
cols = plt.cm.tab20(module[es] % 20)
cols[~intra] = (0, 0, 0, 1)                        # inter-module links in black
wid = np.where(intra, 0.4, 0.9)
order = np.argsort(intra)[::-1]                    # draw intra first, black links on top
fig, ax = plt.subplots(figsize=(12, 12))
gdf.iloc[order].plot(ax=ax, color=cols[order].tolist(), linewidth=wid[order].tolist())
ax.set_axis_off()
ax.set_title(f"{city} road modules (max size {MAX_SIZE}); black = links between modules")
plt.show()

# %%
# ---------- CELL 2: exact detour analysis inside each module ----------
import random, numpy as np, igraph as ig, pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm

# ---- config ----
MIN_MODULE_SIZE = 30        # skip tiny modules
TOP_K_PER_MODULE = None       # highest-betweenness edges tested per module (None = all edges)
SOURCES_PER_MODULE = None    # None = every node in the module is an origin; or e.g. 500
BET_CUTOFF = None            # exact betweenness is fine at module scale
MIN_PAIRS = 30
BIG = 1e12
N_JOBS = -1
SEED = 42

def analyze_module(mid, eids, src, dst, wt, top_k, n_src, bet_cutoff, big, seed):
    # relabel nodes locally
    uniq, inv = np.unique(np.concatenate([src, dst]), return_inverse=True)
    ls, ld = inv[:len(src)], inv[len(src):]
    n = len(uniq)
    g0 = ig.Graph(n=n, edges=list(zip(ls.tolist(), ld.tolist())), directed=True)

    # keep largest strongly connected component only
    comp = np.array(g0.connected_components(mode="strong").membership)
    big_c = np.bincount(comp).argmax()
    keep = (comp[ls] == big_c) & (comp[ld] == big_c)
    eids, ls, ld, wt = eids[keep], ls[keep], ld[keep], wt[keep]
    uniq2, inv2 = np.unique(np.concatenate([ls, ld]), return_inverse=True)
    ls, ld = inv2[:len(ls)], inv2[len(ls):]
    n = len(uniq2)
    if n < 20 or len(eids) < 10:
        return mid, None

    g = ig.Graph(n=n, edges=list(zip(ls.tolist(), ld.tolist())), directed=True)
    wl = wt.tolist()

    bet = np.array(g.edge_betweenness(directed=True, cutoff=bet_cutoff, weights=wl))
    K = len(eids) if top_k is None else min(top_k, len(eids))
    top = np.argsort(bet)[::-1][:K]

    A = np.empty((K, n)); B = np.empty((K, n))
    for k, e in enumerate(top):
        u, v = g.es[int(e)].tuple
        A[k] = g.distances(source=u, weights=wl, mode="in")[0]
        B[k] = g.distances(source=v, weights=wl, mode="out")[0]
    w_top = wt[top]

    rng = random.Random(seed + mid)
    srcs = range(n) if n_src is None or n_src >= n else rng.sample(range(n), n_src)

    s_ratio = np.zeros(K); n_pairs = np.zeros(K, dtype=np.int64); n_cut = np.zeros(K, dtype=np.int64)
    wl2 = list(wl)
    for s in srcs:
        d0 = np.array(g.distances(source=s, weights=wl2, mode="out")[0])
        for k in range(K):
            on = np.isclose(A[k, s] + w_top[k] + B[k], d0, rtol=1e-9, atol=1e-6)
            on[s] = False
            if not on.any():
                continue
            e = int(top[k])
            orig = wl2[e]; wl2[e] = big
            d1 = np.array(g.distances(source=s, weights=wl2, mode="out")[0])
            wl2[e] = orig
            d1t, d0t = d1[on], d0[on]
            cut = d1t >= big / 2
            n_cut[k] += cut.sum()
            ok = ~cut
            s_ratio[k] += (d1t[ok] / d0t[ok]).sum()
            n_pairs[k] += ok.sum()

    return mid, dict(eids=eids[top], bet=bet[top], s_ratio=s_ratio,
                     n_pairs=n_pairs, n_cut=n_cut, n_nodes=n)

# ---- group intra-module edges by module ----
intra_ids = np.where(module[es] == module[ed])[0]
by_mod = {}
for e in intra_ids:
    by_mod.setdefault(int(module[es[e]]), []).append(e)

jobs = [(m, np.array(eids)) for m, eids in by_mod.items() if len(pieces[m]) >= MIN_MODULE_SIZE]
jobs.sort(key=lambda j: len(j[1]), reverse=True)        # biggest first -> better load balance
print(f"Analysing {len(jobs)} modules")

out = Parallel(n_jobs=N_JOBS)(
    delayed(analyze_module)(m, eids, es[eids], ed[eids], w[eids],
                            TOP_K_PER_MODULE, SOURCES_PER_MODULE, BET_CUTOFF, BIG, SEED)
    for m, eids in tqdm(jobs, desc="Modules")
)

# ---- attach results to the graph ----
for e, (u, v, k, d) in enumerate(edges):
    d["mean_detour"] = np.nan; d["n_pairs"] = 0; d["n_cut"] = 0
    d["module"] = int(module[es[e]]) if module[es[e]] == module[ed[e]] else -1

rows = []
for mid, r in out:
    if r is None:
        continue
    vals_m = []
    for j, e in enumerate(r["eids"]):
        d = edges[int(e)][3]
        d["n_pairs"] = int(r["n_pairs"][j]); d["n_cut"] = int(r["n_cut"][j])
        if r["n_pairs"][j] >= MIN_PAIRS:
            d["mean_detour"] = float(r["s_ratio"][j] / r["n_pairs"][j])
            vals_m.append(d["mean_detour"])
    rows.append(dict(module=mid, nodes=r["n_nodes"], tested=len(r["eids"]),
                     with_detour=len(vals_m),
                     mean_detour=np.mean(vals_m) if vals_m else np.nan,
                     median_detour=np.median(vals_m) if vals_m else np.nan,
                     max_detour=np.max(vals_m) if vals_m else np.nan,
                     cut_pairs=int(r["n_cut"].sum())))
# ---- summary (safe when nothing was analysed) ----
cols = ["module", "nodes", "tested", "with_detour",
        "mean_detour", "median_detour", "max_detour", "cut_pairs"]
summary = pd.DataFrame(rows, columns=cols)

if summary.empty:
    print("No module produced results. Check module sizes vs MIN_MODULE_SIZE, "
          "and whether modules are too fragmented after the strongly-connected filter "
          "(n < 20 nodes or < 10 edges returns None).")
else:
    summary = summary.sort_values("mean_detour", ascending=False)
    display(summary.head(20))

n_ok = sum(1 for *_, d in edges if not np.isnan(d["mean_detour"]))
print(f"Edges with a detour value: {n_ok:,}")

# %%
# ---------- CELL 3: Gephi-style dark plot (full, corrected order) ----------
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
import osmnx as ox

# ---- config ----
MODE = "detour"             # "detour" or "module"
BG = "#07090f"
BASE_COLOR = "#4a5d92"      # untested roads / inter-module links
BASE_LW = 0.55
BASE_ALPHA = 0.85
GLOW = True
GLOW_MIN = 0.35             # only edges above this strength get halos (much faster)
SHOW_NODES = False
P_LOW, P_HIGH = 5, 97
BBOX = None                 # (west, south, east, north) in lon/lat
FIGSIZE = (16, 16)
SAVE = "chennai_detour_dark.png"
DPI = 300

# ---- data (same order as G.edges(keys=True)) ----
gdf = ox.convert.graph_to_gdfs(G, nodes=False)
segs = [np.asarray(geom.coords) for geom in gdf.geometry]
E = list(G.edges(keys=True, data=True))
vals = np.array([d.get("mean_detour", np.nan) for *_, d in E], dtype=float)
mods = np.array([d.get("module", -1) for *_, d in E])
has = ~np.isnan(vals)

if BBOX is not None:
    w_, s_, e_, n_ = BBOX
    keep = np.array([(sg[:, 0].max() >= w_) and (sg[:, 0].min() <= e_) and
                     (sg[:, 1].max() >= s_) and (sg[:, 1].min() <= n_) for sg in segs])
else:
    keep = np.ones(len(segs), bool)

# ---- colours ----
gephi = LinearSegmentedColormap.from_list(
    "gephi_glow",
    ["#22b8ff", "#4d6bff", "#9b4dff", "#e0409a", "#ff7a3d", "#ffd23f", "#fffbe6"])

if MODE == "detour":
    lo, hi = np.percentile(vals[has], [P_LOW, P_HIGH])
    lo = max(lo, 1.0)
    norm = Normalize(vmin=lo, vmax=hi, clip=True)
    strength = np.zeros(len(vals)); strength[has] = norm(vals[has])
    edge_cols = np.tile(to_rgba(BASE_COLOR), (len(vals), 1))
    edge_cols[has] = gephi(strength[has])
    emph = has
else:
    rng = np.random.default_rng(7)
    n_mod = max(mods.max() + 1, 1)
    hues = rng.permutation(n_mod) / n_mod
    pal = mpl.colors.hsv_to_rgb(np.c_[hues, np.full(n_mod, 0.65), np.full(n_mod, 1.0)])
    edge_cols = np.tile(to_rgba("#ffffff", 0.9), (len(vals), 1))
    inside = mods >= 0
    edge_cols[inside] = np.c_[pal[mods[inside]], np.full(inside.sum(), 0.9)]
    emph = inside
    strength = np.where(inside, 0.5, 1.0)
    norm = None

# ---- figure ----
fig, ax = plt.subplots(figsize=FIGSIZE, facecolor=BG)
ax.set_facecolor(BG)

def add(mask, lw, alpha=None, z=1):
    idx = np.where(mask & keep)[0]
    if len(idx) == 0:
        return
    c = edge_cols[idx].copy()
    if alpha is not None:
        c[:, 3] = alpha
    ax.add_collection(LineCollection([segs[i] for i in idx], colors=c,
                                     linewidths=lw, capstyle="round", zorder=z))

# 1) base layer (now drawn AFTER add() exists)
if MODE == "detour":
    add(~emph, BASE_LW * 5, alpha=0.05, z=0)       # soft haze
    add(~emph, BASE_LW, alpha=BASE_ALPHA, z=1)
else:
    add(~emph, 1.0, alpha=0.95, z=1)

# 2) tested edges, coolest first so the hottest end up on top
order = np.array([i for i in np.argsort(np.where(emph, strength, -1)) if emph[i] and keep[i]])

def add_sorted(sel, lw, alpha, z):
    if len(sel) == 0:
        return
    c = edge_cols[sel].copy(); c[:, 3] = alpha
    widths = lw * (0.6 + 1.4 * strength[sel])
    ax.add_collection(LineCollection([segs[i] for i in sel], colors=c,
                                     linewidths=widths, capstyle="round", zorder=z))

if GLOW and len(order):
    hot_sel = order[strength[order] > GLOW_MIN]
    add_sorted(hot_sel, 9.0, 0.035, 2)
    add_sorted(hot_sel, 5.0, 0.07, 3)
    add_sorted(hot_sel, 2.8, 0.16, 4)
add_sorted(order, 1.1, 0.95, 5)                    # bright core for ALL tested edges

# 3) glowing nodes on the hottest edges
if SHOW_NODES and MODE == "detour" and len(order):
    hot = order[strength[order] > 0.6]
    if len(hot):
        pts = np.array([segs[i][0] for i in hot] + [segs[i][-1] for i in hot])
        st = np.concatenate([strength[hot], strength[hot]])
        pc = gephi(st)
        for s_mul, a in [(10, 0.05), (5, 0.12), (2, 0.35)]:
            ax.scatter(pts[:, 0], pts[:, 1], s=(8 + 40 * st) * s_mul / 2,
                       c=pc, alpha=a, linewidths=0, zorder=6)
        ax.scatter(pts[:, 0], pts[:, 1], s=4 + 14 * st, c=pc,
                   edgecolors="white", linewidths=0.3, zorder=7)

# sanity check
print(f"Edges in graph: {len(segs):,} | drawn: {int(keep.sum()):,} | "
      f"with detour: {int(has.sum()):,} | untested: {int((~has).sum()):,}")

# ---- view ----
allx = np.concatenate([sg[:, 0] for sg, k in zip(segs, keep) if k])
ally = np.concatenate([sg[:, 1] for sg, k in zip(segs, keep) if k])
ax.set_xlim(allx.min(), allx.max()); ax.set_ylim(ally.min(), ally.max())
ax.set_aspect(1 / np.cos(np.radians(ally.mean())))
ax.set_axis_off()

# ---- title and colourbar ----
title = (f"{city} road network: detour when a critical road is removed"
         if MODE == "detour" else f"{city} road network: modules")
fig.text(0.04, 0.955, title, color="white", fontsize=22, fontweight="bold", ha="left")
fig.text(0.04, 0.928,
         "Colour and thickness = mean shortest-path detour for trips that used the edge"
         if MODE == "detour" else "Colour = module; white = links between modules",
         color="#8b95b3", fontsize=12, ha="left")

if MODE == "detour":
    cax = fig.add_axes([0.04, 0.06, 0.22, 0.012])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=gephi), cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors="#c8d0e8", labelsize=10, length=0)
    cb.set_label("Mean detour ratio (new time / old time)", color="#c8d0e8", fontsize=11, labelpad=8)

if SAVE:
    fig.savefig(SAVE, dpi=DPI, facecolor=BG, bbox_inches="tight", pad_inches=0.3)
plt.show()

