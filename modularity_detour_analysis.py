"""
Python source export of `modularity_detour_analysis.ipynb`.
The notebook is the primary, sectioned version of this analysis.
"""

# %% [markdown]
# # Where does a city's road network break? Detour analysis by community#
# **Question:** if a single road segment is closed, how much longer do the trips that used it take?#
# **Pipeline** (run for Chennai and Bengaluru):#
# 1. **Load** the drivable road network (OSMnx graph) and attach a travel time to every road segment.# 2. **Find communities** with the Leiden algorithm (modularity maximisation), recursively splitting any community larger than a size cap. Roads inside a community are tightly interlinked; roads between communities are the sparse "bridges".# 3. **Stress-test every road inside each community.** For one edge at a time:#    - find all origin–destination pairs whose shortest path uses that edge,#    - remove the edge and recompute the shortest path,#    - record the **detour ratio** = new travel time / old travel time.# 4. **Summarise** each edge by its mean detour ratio over all the trips it carried.# 5. **Plot** the results: one map per city, then both cities side by side with a shared colour scale.#
# > Running the analysis inside communities (instead of on the whole city) keeps each shortest-path problem small enough to be computed *exactly*, and the work is spread across CPU cores.
# %% [markdown]
# ## 1. Imports
# %%
import random

import igraph as ig
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
from IPython.display import display
from joblib import Parallel, delayed
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from tqdm.auto import tqdm

# %% [markdown]
# ## 2. Configuration#
# All tunable parameters live in this one cell.#
# | Group | Parameter | Meaning |# |---|---|---|# | Community detection | `MAX_SIZE` | Hard cap on nodes per community; larger ones are split again. |# | | `RESOLUTION` | Leiden resolution. Higher gives more, smaller communities. |# | | `EDGE_WEIGHT` | `"unit"`: every road link counts equally. `"inv_time"`: faster roads bind nodes more tightly. |# | Detour analysis | `MIN_MODULE_SIZE` | Communities smaller than this are skipped. |# | | `TOP_K_PER_MODULE` | Test only the *k* highest-betweenness edges per community (`None` = every edge). |# | | `SOURCES_PER_MODULE` | Number of origin nodes per community (`None` = every node is an origin). Set e.g. `500` for a faster, sampled run. |# | | `MIN_PAIRS` | An edge needs at least this many affected trips for its detour value to be trusted. |# | Plotting | colours, line widths, glow, `P_LOW`/`P_HIGH` | Visual styling. The colour scale is clipped to these percentiles so a few extreme edges don't wash out the map. |
# %%
# ── Which cities to analyse (each needs a "<city>_roads.graphml" file next to the notebook) ──
CITIES = ["chennai", "bengaluru"]

# ── Community detection ──
MAX_SIZE    = 5000      # hard cap on nodes per community
RESOLUTION  = 1.0       # higher -> more, smaller communities on the first pass
EDGE_WEIGHT = "unit"    # "unit" or "inv_time"
SEED        = 42

# ── Detour analysis ──
MIN_MODULE_SIZE    = 30     # skip tiny communities
TOP_K_PER_MODULE   = None   # None = test every edge; or e.g. 200
SOURCES_PER_MODULE = None   # None = every node is an origin; or e.g. 500
BET_CUTOFF         = None   # exact betweenness is fine at community scale
MIN_PAIRS          = 30     # minimum affected trips for a reliable edge value
BIG                = 1e12   # "infinite" weight used to remove an edge
N_JOBS             = -1     # -1 = use all CPU cores

# ── Caching (the detour analysis is the slow step) ──
USE_CACHE = True            # reuse "<city>_detour_cache.npz" if it exists

# ── Plot style (dark, Gephi-like) ──
BG          = "#07090f"     # background
BASE_COLOR  = "#4a5d92"     # untested roads / links between communities
BASE_LW     = 0.55
BASE_ALPHA  = 0.85
GLOW        = True          # soft halos around strong edges
GLOW_MIN    = 0.35          # only edges above this strength get halos (much faster)
P_LOW, P_HIGH = 5, 97       # colour scale is clipped to these percentiles
DPI         = 300
SHOW_INDIVIDUAL_MAPS = True # also draw one map per city before the combined figure

np.random.seed(SEED)

# Colour map: cool blue (small detour) -> hot white (large detour)
GEPHI = LinearSegmentedColormap.from_list(
    "gephi_glow",
    ["#22b8ff", "#4d6bff", "#9b4dff", "#e0409a", "#ff7a3d", "#ffd23f", "#fffbe6"],
)

# %% [markdown]
# ## 3. Load and prepare a city's road network#
# `load_city` reads the saved GraphML file and keeps only the **largest strongly connected component**, so every node can reach every other node (otherwise "shortest path" is undefined for some pairs). It then adds speeds and travel times (seconds) to every edge.#
# It also builds flat NumPy arrays (edge endpoints, travel times, node coordinates). The position of an edge in `edges` is its **edge id** and is used consistently in the rest of the notebook.#
# `download_city` is an optional helper for fetching and saving a city the first time.
# %%
def download_city(place, city):
    '''One-off helper: download a drivable network from OpenStreetMap and save it as GraphML.'''
    G = ox.graph_from_place(place, network_type="drive", simplify=True)
    ox.save_graphml(G, f"{city}_roads.graphml")
    return G


def load_city(city):
    '''Load, clean and flatten the road network of one city into a dict ("net").'''
    G = ox.load_graphml(f"{city}_roads.graphml")
    G = ox.truncate.largest_component(G, strongly=True)
    G = ox.add_edge_speeds(G)
    G = ox.add_edge_travel_times(G)          # adds "travel_time" in seconds

    nodes = list(G.nodes)
    idx   = {n: i for i, n in enumerate(nodes)}
    edges = list(G.edges(keys=True, data=True))    # list position == edge id

    gdf = ox.convert.graph_to_gdfs(G, nodes=False)  # same order as G.edges(keys=True)
    assert len(gdf) == len(edges), "edge order mismatch between G and its GeoDataFrame"

    return dict(
        city=city, G=G, nodes=nodes, edges=edges, N=len(nodes),
        es=np.array([idx[u] for u, v, k, d in edges]),                 # edge start (node index)
        ed=np.array([idx[v] for u, v, k, d in edges]),                 # edge end   (node index)
        w=np.array([d["travel_time"] for *_, d in edges], dtype=float),
        xs=np.array([G.nodes[n]["x"] for n in nodes]),                 # longitude
        ys=np.array([G.nodes[n]["y"] for n in nodes]),                 # latitude
        gdf=gdf,
        segs=[np.asarray(geom.coords) for geom in gdf.geometry],       # polylines for plotting
    )

# %% [markdown]
# ## 4. Community detection (Leiden + size cap)#
# **Modularity** `Q` measures how much denser the connections *inside* groups are than expected by chance. A value near 1 means very clear communities.#
# Steps in `detect_modules`:#
# 1. Build an **undirected, simplified** copy of the network (community detection ignores one-way streets and parallel edges).# 2. Run **Leiden** on the whole city.# 3. Any community still larger than `MAX_SIZE` is split **recursively**: Leiden is re-run on just that community, doubling the resolution until it actually splits. As a last resort it is cut in half spatially along its longer axis.# 4. Communities are renumbered from largest to smallest, and `Q` plus the share of edges inside communities are reported.
# %%
def _leiden(graph, resolution):
    '''Run Leiden (modularity objective) and return the membership vector.'''
    return np.array(graph.community_leiden(
        objective_function="modularity", weights="cw",
        resolution=resolution, n_iterations=-1).membership)


def _split_module(gu, vids, resolution, xs, ys, max_size):
    '''Recursively split a set of node ids until every piece has <= max_size nodes.'''
    vids = np.sort(np.asarray(vids))
    if len(vids) <= max_size:
        return [vids]

    sub = gu.induced_subgraph(vids.tolist())
    r, memb = resolution, None
    for _ in range(10):                       # raise the resolution until it really splits
        memb = _leiden(sub, r)
        if memb.max() > 0:
            break
        r *= 2

    if memb.max() == 0:                       # fallback: cut in half spatially
        axis = xs if np.ptp(xs[vids]) >= np.ptp(ys[vids]) else ys
        memb = (axis[vids] > np.median(axis[vids])).astype(int)
        if memb.min() == memb.max():
            memb[: len(memb) // 2] = 1

    pieces = []
    for c in np.unique(memb):
        pieces += _split_module(gu, vids[memb == c], resolution * 1.5, xs, ys, max_size)
    return pieces


def detect_modules(net):
    '''Assign every node to a community. Returns (module_of_node, list_of_node_id_arrays).'''
    N, es, ed, w = net["N"], net["es"], net["ed"], net["w"]

    # Undirected weighted graph for community detection
    gu = ig.Graph(n=N, edges=list(zip(es.tolist(), ed.tolist())), directed=False)
    gu.es["cw"] = np.ones(len(es)) if EDGE_WEIGHT == "unit" else (1.0 / w).tolist()
    gu.simplify(multiple=True, loops=True, combine_edges="max")

    # First pass on the whole city, then enforce the size cap
    memb0 = _leiden(gu, RESOLUTION)
    print(f"[{net['city']}] first pass: {memb0.max() + 1} modules, "
          f"Q = {gu.modularity(memb0.tolist(), weights='cw'):.3f}")

    pieces = []
    for c in np.unique(memb0):
        pieces += _split_module(gu, np.where(memb0 == c)[0], RESOLUTION,
                                net["xs"], net["ys"], MAX_SIZE)

    pieces.sort(key=len, reverse=True)       # module 0 = largest
    module = np.empty(N, dtype=int)
    for m, vids in enumerate(pieces):
        module[vids] = m

    sizes = np.array([len(p) for p in pieces])
    intra = module[es] == module[ed]
    print(f"[{net['city']}] final: {len(pieces)} modules "
          f"(max {sizes.max()}, median {int(np.median(sizes))}, min {sizes.min()}), "
          f"Q = {gu.modularity(module.tolist(), weights='cw'):.3f}")
    print(f"[{net['city']}] edges inside modules: {intra.mean():.1%}, "
          f"between modules: {(~intra).mean():.1%}")
    return module, pieces

# %% [markdown]
# ## 5. Exact detour analysis inside one community#
# `analyze_module` is the core of the project. For one community it:#
# 1. **Re-indexes** the nodes locally and keeps only the largest strongly connected part.# 2. Computes **edge betweenness** (how many shortest paths use each edge) and picks the edges to test (all of them by default).# 3. **Pre-computes** two distance tables for each tested edge `u→v`: the distance *from every node to `u`* (`A`) and *from `v` to every node* (`B`).# 4. For every **origin** `s`:#    - Compute shortest-path distances `d0` from `s` to all nodes.#    - An edge lies on a shortest path from `s` to `t` exactly when `dist(s,u) + w(u→v) + dist(v,t) == d0[t]`. This test selects the **affected destinations** without enumerating paths.#    - **Remove** the edge (weight set to `BIG`), recompute the distances, and compare. Destinations that become unreachable are counted as *cut*; the others contribute their **detour ratio** `new_time / old_time`.# 5. Returns, per edge, the sum of ratios, the number of affected pairs and the number of cut pairs. The mean detour is `sum / pairs`.#
# The function is deliberately self-contained (plain arrays in, dict out) so it can be shipped to worker processes by `joblib`.
# %%
def analyze_module(mid, eids, src, dst, wt, top_k, n_src, bet_cutoff, big, seed):
    '''Detour statistics for every (tested) edge inside one community.

    mid  : community id            eids : global edge ids of this community's edges
    src, dst, wt : edge endpoints (global node idx) and travel times
    Returns (mid, results_dict) or (mid, None) if the community is too small.
    '''
    # ── 1. Relabel nodes locally (0..n-1) ──
    uniq, inv = np.unique(np.concatenate([src, dst]), return_inverse=True)
    ls, ld = inv[:len(src)], inv[len(src):]
    n = len(uniq)
    g0 = ig.Graph(n=n, edges=list(zip(ls.tolist(), ld.tolist())), directed=True)

    # Keep only the largest strongly connected component
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

    # ── 2. Choose which edges to test (highest betweenness first) ──
    bet = np.array(g.edge_betweenness(directed=True, cutoff=bet_cutoff, weights=wl))
    K = len(eids) if top_k is None else min(top_k, len(eids))
    top = np.argsort(bet)[::-1][:K]

    # ── 3. Distance tables: A[k] = dist(all -> u_k),  B[k] = dist(v_k -> all) ──
    A = np.empty((K, n)); B = np.empty((K, n))
    for k, e in enumerate(top):
        u, v = g.es[int(e)].tuple
        A[k] = g.distances(source=u, weights=wl, mode="in")[0]
        B[k] = g.distances(source=v, weights=wl, mode="out")[0]
    w_top = wt[top]

    # ── 4. Loop over origins ──
    rng = random.Random(seed + mid)
    srcs = range(n) if n_src is None or n_src >= n else rng.sample(range(n), n_src)

    s_ratio = np.zeros(K)                       # sum of detour ratios per edge
    n_pairs = np.zeros(K, dtype=np.int64)       # affected (origin, destination) pairs
    n_cut   = np.zeros(K, dtype=np.int64)       # pairs disconnected by the removal
    wl2 = list(wl)                              # working copy of the weights

    for s in srcs:
        d0 = np.array(g.distances(source=s, weights=wl2, mode="out")[0])
        for k in range(K):
            # Destinations whose shortest path from s passes through edge k
            on = np.isclose(A[k, s] + w_top[k] + B[k], d0, rtol=1e-9, atol=1e-6)
            on[s] = False
            if not on.any():
                continue

            # Remove the edge, recompute, restore
            e = int(top[k])
            orig = wl2[e]; wl2[e] = big
            d1 = np.array(g.distances(source=s, weights=wl2, mode="out")[0])
            wl2[e] = orig

            d1t, d0t = d1[on], d0[on]
            cut = d1t >= big / 2                # unreachable without the edge
            n_cut[k] += cut.sum()
            ok = ~cut
            s_ratio[k] += (d1t[ok] / d0t[ok]).sum()
            n_pairs[k] += ok.sum()

    return mid, dict(eids=eids[top], bet=bet[top], s_ratio=s_ratio,
                     n_pairs=n_pairs, n_cut=n_cut, n_nodes=n)

# %% [markdown]
# ## 6. Run the analysis on all communities of a city#
# `run_detour_analysis`:#
# 1. Groups the **intra-community** edges by community. Links *between* communities are not tested, so they remain "untested" (grey-blue) on the maps.# 2. Runs `analyze_module` for every community **in parallel** (largest first, for better load balancing).# 3. Converts the sums into a **mean detour ratio per edge**, keeping only edges with at least `MIN_PAIRS` affected trips.# 4. Builds a per-community **summary table**.#
# Results are returned as two arrays aligned with the edge ids: `detour` (NaN = untested) and `edge_module` (community id, or `-1` for inter-community edges). They are also written onto the graph edges as `mean_detour`, `n_pairs`, `n_cut` and `module`.
# %%
def run_detour_analysis(net, module, pieces):
    '''Run analyze_module on every community. Returns (detour, edge_module, summary_df).'''
    edges, es, ed, w = net["edges"], net["es"], net["ed"], net["w"]
    E = len(edges)

    # Group intra-community edges by community
    by_mod = {}
    for e in np.where(module[es] == module[ed])[0]:
        by_mod.setdefault(int(module[es[e]]), []).append(e)

    jobs = [(m, np.array(eids)) for m, eids in by_mod.items()
            if len(pieces[m]) >= MIN_MODULE_SIZE]
    jobs.sort(key=lambda j: len(j[1]), reverse=True)     # biggest first
    print(f"[{net['city']}] analysing {len(jobs)} modules")

    out = Parallel(n_jobs=N_JOBS)(
        delayed(analyze_module)(m, eids, es[eids], ed[eids], w[eids],
                                TOP_K_PER_MODULE, SOURCES_PER_MODULE,
                                BET_CUTOFF, BIG, SEED)
        for m, eids in tqdm(jobs, desc=f"{net['city']} modules")
    )

    # Defaults for every edge
    detour = np.full(E, np.nan)
    edge_module = np.where(module[es] == module[ed], module[es], -1)
    for e, (u, v, k, d) in enumerate(edges):
        d["mean_detour"], d["n_pairs"], d["n_cut"] = np.nan, 0, 0
        d["module"] = int(edge_module[e])

    # Fill in the tested edges and collect per-community statistics
    rows = []
    for mid, r in out:
        if r is None:
            continue
        vals = []
        for j, e in enumerate(r["eids"]):
            e = int(e)
            edges[e][3]["n_pairs"] = int(r["n_pairs"][j])
            edges[e][3]["n_cut"]   = int(r["n_cut"][j])
            if r["n_pairs"][j] >= MIN_PAIRS:
                detour[e] = r["s_ratio"][j] / r["n_pairs"][j]
                edges[e][3]["mean_detour"] = float(detour[e])
                vals.append(detour[e])
        rows.append(dict(
            module=mid, nodes=r["n_nodes"], tested=len(r["eids"]), with_detour=len(vals),
            mean_detour=np.mean(vals) if vals else np.nan,
            median_detour=np.median(vals) if vals else np.nan,
            max_detour=np.max(vals) if vals else np.nan,
            cut_pairs=int(r["n_cut"].sum())))

    columns = ["module", "nodes", "tested", "with_detour",
               "mean_detour", "median_detour", "max_detour", "cut_pairs"]
    summary = pd.DataFrame(rows, columns=columns)
    if not summary.empty:
        summary = summary.sort_values("mean_detour", ascending=False)

    print(f"[{net['city']}] edges with a detour value: {int(np.isfinite(detour).sum()):,} of {E:,}")
    return detour, edge_module, summary

# %% [markdown]
# ## 7. Run both cities#
# For each city this cell loads the network, detects communities, runs the detour analysis and stores everything in the `results` dictionary. If `USE_CACHE` is on and a cache file exists, the slow analysis step is skipped and the saved detour values are reused.#
# > ⚠️ With `TOP_K_PER_MODULE = None` and `SOURCES_PER_MODULE = None` the analysis is exact but can take a long time. For a quick trial, set e.g. `TOP_K_PER_MODULE = 200` and `SOURCES_PER_MODULE = 300` in the configuration cell.
# %%
from pathlib import Path

results = {}

for city in CITIES:
    print(f"\n{'=' * 20}  {city.upper()}  {'=' * 20}")
    net = load_city(city)
    print(f"[{city}] nodes: {net['N']:,}   edges: {len(net['edges']):,}")

    module, pieces = detect_modules(net)

    cache = Path(f"{city}_detour_cache.npz")
    if USE_CACHE and cache.exists():
        saved = np.load(cache)
        detour, edge_module = saved["detour"], saved["edge_module"]
        assert len(detour) == len(net["edges"]), "cache does not match this graph - delete it"
        summary = pd.DataFrame()
        print(f"[{city}] loaded cached detour values from {cache}")
    else:
        detour, edge_module, summary = run_detour_analysis(net, module, pieces)
        if USE_CACHE:
            np.savez_compressed(cache, detour=detour, edge_module=edge_module)

    if not summary.empty:
        print(f"[{city}] 20 communities with the largest mean detour:")
        display(summary.head(20))

    results[city] = dict(net=net, module=module, pieces=pieces,
                         detour=detour, edge_module=edge_module, summary=summary)

# %% [markdown]
# ## 8. Plotting helpers#
# * `plot_modules` colours each road by its community; black roads are the links *between* communities.# * `make_norm` builds the colour scale from detour values. Passing several cities at once gives **one shared scale**, so colours mean the same thing on every map.# * `draw_detour_map` paints a detour map onto a given axis in layers: a faint haze for untested roads, soft glow halos for the strongest edges, and a bright core for every tested edge. **Colour and thickness both grow with the detour ratio**, and the hottest edges are drawn last so they stay on top.
# %%
def plot_modules(res, max_size=None):
    '''Colour every road by its community; black = links between communities.'''
    net, module = res["net"], res["module"]
    es, ed = net["es"], net["ed"]
    intra = module[es] == module[ed]

    colours = plt.cm.tab20(module[es] % 20)
    colours[~intra] = (0, 0, 0, 1)
    widths = np.where(intra, 0.4, 0.9)
    order = np.argsort(intra)[::-1]              # draw intra links first, black links on top

    fig, ax = plt.subplots(figsize=(12, 12))
    net["gdf"].iloc[order].plot(ax=ax, color=colours[order].tolist(),
                                linewidth=widths[order].tolist())
    ax.set_axis_off()
    ax.set_title(f"{net['city'].title()} road communities "
                 f"(max size {max_size or MAX_SIZE}); black = links between communities")
    plt.show()


def make_norm(*detours):
    '''Shared colour normalisation from one or more detour arrays (NaN = untested).'''
    vals = np.concatenate([d[np.isfinite(d)] for d in detours])
    lo, hi = np.percentile(vals, [P_LOW, P_HIGH])
    return Normalize(vmin=max(lo, 1.0), vmax=hi, clip=True)


def draw_detour_map(ax, net, detour, norm, lw_scale=1.0):
    '''Draw one city's detour map on `ax` using the given colour normalisation.'''
    segs = net["segs"]
    has = np.isfinite(detour)

    # Colour + "strength" (0..1) of every edge
    strength = np.zeros(len(detour))
    strength[has] = norm(detour[has])
    colours = np.tile(to_rgba(BASE_COLOR), (len(detour), 1))
    colours[has] = GEPHI(strength[has])

    def add_layer(idx, linewidth, alpha, z):
        if len(idx) == 0:
            return
        c = colours[idx].copy()
        c[:, 3] = alpha
        ax.add_collection(LineCollection([segs[i] for i in idx], colors=c,
                                         linewidths=linewidth, capstyle="round", zorder=z))

    # 1) Untested roads: faint haze + thin line
    untested = np.where(~has)[0]
    add_layer(untested, BASE_LW * 5 * lw_scale, 0.05, 0)
    add_layer(untested, BASE_LW * lw_scale, BASE_ALPHA, 1)

    # 2) Tested roads, coolest first so the hottest are drawn on top
    tested = np.where(has)[0]
    tested = tested[np.argsort(strength[tested])]

    def width(idx, base):                        # thickness grows with strength
        return base * lw_scale * (0.6 + 1.4 * strength[idx])

    if GLOW:
        hot = tested[strength[tested] > GLOW_MIN]
        add_layer(hot, width(hot, 9.0), 0.035, 2)
        add_layer(hot, width(hot, 5.0), 0.07, 3)
        add_layer(hot, width(hot, 2.8), 0.16, 4)
    add_layer(tested, width(tested, 1.1), 0.95, 5)   # bright core for every tested edge

    # 3) View: tight bounds, correct lon/lat aspect ratio
    allx = np.concatenate([s[:, 0] for s in segs])
    ally = np.concatenate([s[:, 1] for s in segs])
    ax.set_xlim(allx.min(), allx.max())
    ax.set_ylim(ally.min(), ally.max())
    ax.set_aspect(1 / np.cos(np.radians(ally.mean())))
    ax.set_facecolor(BG)
    ax.set_axis_off()

    print(f"[{net['city']}] edges drawn: {len(segs):,} | with detour: {int(has.sum()):,} "
          f"| untested: {int((~has).sum()):,}")


def add_colorbar(fig, cax_rect, norm):
    '''Horizontal colour bar shared by every panel of a figure.'''
    cax = fig.add_axes(cax_rect)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=GEPHI),
                      cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors="#c8d0e8", labelsize=10, length=0)
    cb.set_label("Mean detour ratio (new travel time / old travel time)",
                 color="#c8d0e8", fontsize=11, labelpad=8)
    return cb

# %% [markdown]
# ## 9. Community maps#
# The communities found for each city. Every colour is one community; black lines are the sparse links connecting them.
# %%
for city in CITIES:
    plot_modules(results[city])

# %% [markdown]
# ## 10. Detour map for each city#
# Each city is drawn with **its own** colour scale here (so you can see its internal structure clearly). The combined figure in the next section uses a single shared scale so the two cities can be compared directly.
# %%
if SHOW_INDIVIDUAL_MAPS:
    for city in CITIES:
        res = results[city]
        norm = make_norm(res["detour"])

        fig, ax = plt.subplots(figsize=(16, 16), facecolor=BG)
        draw_detour_map(ax, res["net"], res["detour"], norm)
        fig.text(0.04, 0.955, f"{city.title()} road network: detour when a critical road is removed",
                 color="white", fontsize=22, fontweight="bold", ha="left")
        fig.text(0.04, 0.928, "Colour and thickness = mean shortest-path detour for trips that used the edge",
                 color="#8b95b3", fontsize=12, ha="left")
        add_colorbar(fig, [0.04, 0.06, 0.22, 0.012], norm)

        fig.savefig(f"{city}_detour_dark.png", dpi=DPI, facecolor=BG,
                    bbox_inches="tight", pad_inches=0.3)
        plt.show()

# %% [markdown]
# ## 11. Chennai and Bengaluru side by side#
# Both cities are drawn in one figure with a **common colour scale and a single legend (colour bar)**. The scale is computed from the pooled detour values of both cities (clipped to the `P_LOW`–`P_HIGH` percentiles), so the same colour means the same detour ratio in both panels. The figure is saved as PNG and SVG in `modularity_detour_output/`.
# %%
# One normalisation for both cities -> identical colours mean identical detours
from pathlib import Path

OUTPUT_DIR = Path("modularity_detour_output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

shared_norm = make_norm(*[results[c]["detour"] for c in CITIES])

fig, axes = plt.subplots(1, len(CITIES), figsize=(14 * len(CITIES), 15), facecolor=BG)
fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.10, wspace=0.03)

for ax, city in zip(np.atleast_1d(axes), CITIES):
    draw_detour_map(ax, results[city]["net"], results[city]["detour"],
                    shared_norm, lw_scale=1.15)
    ax.set_anchor("C")                  # centre each map in its panel (cities differ in shape)
    ax.apply_aspect()
    pos = ax.get_position()             # panel title at a common height, centred over the map
    fig.text((pos.x0 + pos.x1) / 2, 0.885, city.title(), color="white",
             fontsize=26, fontweight="bold", ha="center")

fig.suptitle("Detour when a critical road is removed", color="white",
             fontsize=30, fontweight="bold", x=0.5, y=0.975)
fig.text(0.5, 0.925,
         "Colour and thickness = mean shortest-path detour for trips that used the road  "
         "(one colour scale for both cities)",
         color="#8b95b3", fontsize=14, ha="center")

add_colorbar(fig, [0.36, 0.045, 0.28, 0.012], shared_norm)   # common legend

fig.savefig(OUTPUT_DIR / "chennai_bengaluru_detour_dark.png", dpi=DPI, facecolor=BG,
            bbox_inches="tight", pad_inches=0.3)
fig.savefig(OUTPUT_DIR / "chennai_bengaluru_detour_dark.svg", facecolor=BG,
            bbox_inches="tight", pad_inches=0.3)
plt.show()

# %% [markdown]
# ## 12. Save plots to directory
# %%
# ── Save every plotted output (and its data) into one folder ──
import shutil
from pathlib import Path

OUTPUT_DIR = Path("modularity_detour_output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for city in CITIES:
    res, net = results[city], results[city]["net"]
    module, es, ed = res["module"], net["es"], net["ed"]

    # 1) Community map (redrawn: section 9 only displayed it)
    intra = module[es] == module[ed]
    colours = plt.cm.tab20(module[es] % 20)
    colours[~intra] = (0, 0, 0, 1)                      # black = links between communities
    order = np.argsort(intra)[::-1]
    fig, ax = plt.subplots(figsize=(12, 12))
    net["gdf"].iloc[order].plot(ax=ax, color=colours[order].tolist(),
                                linewidth=np.where(intra, 0.4, 0.9)[order].tolist())
    ax.set_axis_off()
    ax.set_title(f"{city.title()} road communities; black = links between communities")
    fig.savefig(OUTPUT_DIR / f"{city}_communities.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # 2) Data behind the maps: per-edge detour values and the per-community summary
    pd.DataFrame({"u": [e[0] for e in net["edges"]], "v": [e[1] for e in net["edges"]],
                  "travel_time_s": net["w"], "module": res["edge_module"],
                  "mean_detour": res["detour"]}).to_csv(OUTPUT_DIR / f"{city}_edge_detours.csv", index=False)
    if not res["summary"].empty:                        # empty when results came from the cache
        res["summary"].to_csv(OUTPUT_DIR / f"{city}_module_summary.csv", index=False)

# 3) Detour maps (already saved by sections 10 and 11) -> copy into the same folder
for f in Path(".").glob("*_detour_dark.png"):
    shutil.copy(f, OUTPUT_DIR / f.name)

print("Saved to", OUTPUT_DIR.resolve())
for f in sorted(OUTPUT_DIR.iterdir()):
    print("  ", f.name)

