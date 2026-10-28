#Author: Omkar

# =============================================================================
# ANALYSIS 2: SMALL-WORLD AND PATH STRUCTURE OF THE ROAD NETWORKS
# Cities: Chennai and Bengaluru (OpenStreetMap drive network via osmnx)
#
# PURPOSE
#   Measure how "far apart" things are in each road network (shortest-path
#   length in hops and in metres, diameter) and how clustered it is, and test
#   whether the networks are small-world by comparing them with an
#   Erdos-Renyi (ER) graph of the same N and number of links.
#
# REASON
#   My robustness experiments (random failure vs betweenness attack) showed how
#   the networks break. This analysis explains WHY: a network whose paths are
#   long and that has no shortcuts or hubs depends on a few corridors, so
#   removing the high-betweenness nodes should hurt it a lot.
#   It also links to the course material: small-world networks have
#   L ~ ln N, while spatial/lattice-like networks have L ~ sqrt(N).
#
# EXPECTATION (hypotheses to check, not results)
#   1. Road networks are NOT small-world: their average path length is many
#      times longer than the ER network's (tens to ~100 hops against ~10 for ER).
#   2. Clustering is low but higher than ER's (short loops of a few streets), so
#      C/C_er is large while L/L_er is also large, giving a small-world index
#      sigma = (C/C_er) / (L/L_er) that is not far above 1 or even below it.
#   3. The diameter is huge compared with ER because there are no shortcuts.
#   4. On the maps, hop-distance colours grow in smooth rings from the source,
#      which is the signature of a spatial network. In ER they would all be
#      within a few hops.
#   5. The city with longer paths / less clustering should be the more
#      fragile one under targeted attack.
#
# Run:   pip install osmnx networkx numpy matplotlib
#        python road_path_structure.py
# =============================================================================

import os
import numpy as np
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

CITIES = {
    "Chennai":   "Chennai, Tamil Nadu, India",
    "Bengaluru": "Bengaluru, Karnataka, India",
}
SEED = 42
N_SOURCES = 200          # BFS/Dijkstra from this many random sources (sampling)
N_ITERATIONS = 2         # number of visualised "iterations" (random sources)
OUT_DIR = "path_structure_plots"

# ---- display settings ----
# For a zoomable window, uncomment the next two lines BEFORE importing pyplot
# (needs:  pip install pyqt5):
# import matplotlib
# matplotlib.use("QtAgg")
EDGE_ONLY = True       # True: colour the road segments, draw no node dots
LINE_W = 0.7           # road line width
DPI = 300              # saved PNG resolution (zoom into the saved file)
ZOOM = None            # e.g. (80.20, 80.30, 12.95, 13.10) = (lon_min, lon_max, lat_min, lat_max)
C_ROAD, C_ER = "#FF7A1A", "#2F6BFF"

rng = np.random.default_rng(SEED)
os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# 1. LOAD NETWORK (undirected, simple, largest component, with lengths + x,y)
# ----------------------------------------------------------------------------
def load_city(place):
    Gm = ox.graph_from_place(place, network_type="drive")
    U = nx.Graph()
    for u, v, d in Gm.edges(data=True):
        if u == v:
            continue
        L = d.get("length", 1.0)
        if U.has_edge(u, v):
            U[u][v]["length"] = min(U[u][v]["length"], L)
        else:
            U.add_edge(u, v, length=L)
    for n, d in Gm.nodes(data=True):
        if n in U:
            U.nodes[n]["x"], U.nodes[n]["y"] = d["x"], d["y"]
    U = U.subgraph(max(nx.connected_components(U), key=len)).copy()
    return nx.convert_node_labels_to_integers(U)


# ----------------------------------------------------------------------------
# 2. PATH STATISTICS (sampled sources, so it stays fast on 20k+ nodes)
# ----------------------------------------------------------------------------
def path_stats(G, weight=None):
    """Mean shortest path, and a lower bound on the diameter, from sampled sources."""
    nodes = np.array(G.nodes)
    srcs = rng.choice(nodes, min(N_SOURCES, len(nodes)), replace=False)
    all_d, ecc = [], 0
    for s in srcs:
        if weight is None:
            d = nx.single_source_shortest_path_length(G, s)
        else:
            d = nx.single_source_dijkstra_path_length(G, s, weight=weight)
        vals = np.array([v for v in d.values() if v > 0])
        all_d.append(vals)
        ecc = max(ecc, vals.max())
    all_d = np.concatenate(all_d)
    return all_d.mean(), ecc, all_d


def lcc(G):
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


# ----------------------------------------------------------------------------
# 3. MAP DRAWING
# ----------------------------------------------------------------------------
def edge_segments(G):
    return [[(G.nodes[u]["x"], G.nodes[u]["y"]), (G.nodes[v]["x"], G.nodes[v]["y"])]
            for u, v in G.edges]


def draw_iteration(ax, G, segs, source, it):
    """One iteration: pick a source, BFS hop distances, farthest node, shortest path."""
    hops = nx.single_source_shortest_path_length(G, source)
    far = max(hops, key=hops.get)
    path = nx.shortest_path(G, source, far)
    nodes = np.array(list(G.nodes))
    xs = np.array([G.nodes[n]["x"] for n in nodes])
    ys = np.array([G.nodes[n]["y"] for n in nodes])
    h = np.array([hops[n] for n in nodes])

    if EDGE_ONLY:   # colour every road segment by the mean hop distance of its two ends
        hv = np.array([(hops[u] + hops[v]) / 2 for u, v in G.edges])
        sc = LineCollection(segs, array=hv, cmap="viridis", linewidths=LINE_W, zorder=2)
        ax.add_collection(sc)
    else:
        ax.add_collection(LineCollection(segs, colors="#cccccc", linewidths=0.3, zorder=1))
        sc = ax.scatter(xs, ys, c=h, s=1.5, cmap="viridis", zorder=2)
    px = [G.nodes[n]["x"] for n in path]
    py = [G.nodes[n]["y"] for n in path]
    ax.plot(px, py, color="red", lw=2.2, zorder=3, label=f"longest hop-path ({hops[far]} hops)")
    ax.scatter(px[0], py[0], marker="*", s=180, color="red", edgecolor="k", zorder=4)
    ax.scatter(px[-1], py[-1], marker="X", s=90, color="red", edgecolor="k", zorder=4)
    ax.set_aspect(1 / np.cos(np.radians(ys.mean())))
    ax.autoscale()
    if ZOOM:
        ax.set_xlim(ZOOM[0], ZOOM[1]); ax.set_ylim(ZOOM[2], ZOOM[3])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Iteration {it}: source (star) -> farthest node (X)\n"
                 f"mean hops from source = {h.mean():.0f}, max = {h.max()}", fontsize=10)
    ax.legend(loc="lower left", fontsize=8, frameon=True)
    return sc


# ----------------------------------------------------------------------------
# 4. MAIN
# ----------------------------------------------------------------------------
def analyse(city, G):
    n, m = G.number_of_nodes(), G.number_of_edges()
    print(f"\n=== {city}: N={n:,}  M={m:,}  <k>={2*m/n:.2f} ===")

    ER = lcc(nx.gnm_random_graph(n, m, seed=SEED))

    L_hops, diam, hops_all = path_stats(G)
    L_m, _, dist_all = path_stats(G, weight="length")
    C = nx.average_clustering(G)

    L_er, diam_er, hops_er = path_stats(ER)
    C_er = nx.average_clustering(ER)
    sigma = (C / C_er) / (L_hops / L_er)

    res = dict(N=n, k=2 * m / n, L_hops=L_hops, L_km=L_m / 1000, diam=diam, C=C,
               L_er=L_er, diam_er=diam_er, C_er=C_er, sigma=sigma,
               sqrtN=np.sqrt(n), lnN=np.log(n), hops=hops_all, hops_er=hops_er)

    # ---- maps: N_ITERATIONS random sources, side by side ----
    segs = edge_segments(G)
    fig, axes = plt.subplots(1, N_ITERATIONS, figsize=(6.5 * N_ITERATIONS, 6.5))
    axes = np.atleast_1d(axes)
    for it, ax in enumerate(axes, 1):
        source = int(rng.choice(np.array(G.nodes)))
        sc = draw_iteration(ax, G, segs, source, it)
    fig.colorbar(sc, ax=axes, shrink=0.7, label="hop distance from source")
    fig.suptitle(f"{city}: hop distance from a random source (BFS), "
                 f"{N_ITERATIONS} iterations", fontsize=13, fontweight="bold")
    fig.savefig(os.path.join(OUT_DIR, f"{city}_hop_maps.png"), dpi=DPI, bbox_inches="tight")
    return res


def main():
    results = {}
    for city, place in CITIES.items():
        results[city] = analyse(city, load_city(place))

    # ---- path-length distributions: road vs ER ----
    fig, axes = plt.subplots(1, len(CITIES), figsize=(6.5 * len(CITIES), 4.5))
    for ax, (city, r) in zip(np.atleast_1d(axes), results.items()):
        bins = np.arange(0, max(r["hops"].max(), r["hops_er"].max()) + 2)
        ax.hist(r["hops"], bins=bins, density=True, color=C_ROAD, alpha=0.8,
                label=f"{city} roads (mean {r['L_hops']:.1f})")
        ax.hist(r["hops_er"], bins=bins, density=True, color=C_ER, alpha=0.8,
                label=f"ER same N, M (mean {r['L_er']:.1f})")
        ax.set_xlabel("shortest path length (hops)"); ax.set_ylabel("fraction of node pairs")
        ax.set_title(f"{city}: path length distribution", fontweight="bold")
        ax.legend(frameon=False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "path_length_distributions.png"), dpi=200)

    # ---- summary table ----
    print("\n" + "=" * 78)
    print(f"{'':<34}" + "".join(f"{c:>22}" for c in CITIES))
    print("-" * 78)
    rows = [("N (nodes)", "N", "{:,.0f}"), ("mean degree <k>", "k", "{:.2f}"),
            ("avg path length (hops)", "L_hops", "{:.1f}"),
            ("avg path length (km)", "L_km", "{:.2f}"),
            ("diameter (hops, lower bound)", "diam", "{:.0f}"),
            ("clustering C", "C", "{:.4f}"),
            ("ER: avg path (hops)", "L_er", "{:.1f}"),
            ("ER: diameter (hops)", "diam_er", "{:.0f}"),
            ("ER: clustering", "C_er", "{:.5f}"),
            ("small-world index sigma", "sigma", "{:.2f}"),
            ("reference: sqrt(N)", "sqrtN", "{:.0f}"),
            ("reference: ln N", "lnN", "{:.1f}")]
    for label, key, fmt in rows:
        print(f"{label:<34}" + "".join(f"{fmt.format(results[c][key]):>22}" for c in CITIES))
    print("=" * 78)
    print("Reading it: if avg path (hops) is close to sqrt(N)-scale and far above the ER value,")
    print("the network is lattice-like (spatial), not small-world, even if C is high.")
    print(f"\nPlots saved in ./{OUT_DIR}/")
    plt.show()


if __name__ == "__main__":
    main()