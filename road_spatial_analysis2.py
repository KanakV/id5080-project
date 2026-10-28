#Author: Omkar

# =============================================================================
# ANALYSIS 4: SPATIAL AND GEOGRAPHIC ANALYSIS OF THE ROAD NETWORKS
# Cities: Chennai and Bengaluru (OpenStreetMap drive network via osmnx)
#
# PURPOSE
#   Describe the geometry of each road network, not just its graph topology:
#     (a) street orientation entropy   - is the street layout grid-like or organic?
#     (b) circuity                     - how much longer is the real route than
#                                        the straight line between two places?
#     (c) intersection / street density per km2
#     (d) road-class mix               - share of main roads vs local streets
#
# REASON
#   The robustness experiments treat the road network as an abstract graph.
#   Roads are embedded in space, so geometry decides how redundant a network
#   really is. A grid has many equal alternative routes, so one blocked road
#   barely matters. An organic network with high circuity and a few main
#   corridors has fewer alternatives. These measures give a physical
#   explanation for why one city came out more or less robust than the other.
#
# EXPECTATION (hypotheses to check, not results)
#   1. Both cities have high orientation entropy (organic, grown over
#      centuries). The city with larger planned/gridded areas should show
#      lower entropy and a few strong peaks (N-S / E-W) in the polar plot.
#   2. Circuity is above 1 everywhere, roughly 1.2 to 1.5. The more circuitous
#      city has longer detours when a main road is blocked, so it should be
#      the more fragile one.
#   3. Intersection density is higher in the more compact, grid-like city.
#      More intersections per km2 means more alternative routes and more
#      robustness.
#   4. Most road length is local/residential streets; the small share of
#      primary/trunk roads carries most traffic. This matches the
#      "few high-betweenness corridors" seen in the betweenness attack.
#
# Run:   pip install osmnx networkx numpy matplotlib geopandas
#        python road_spatial_analysis.py
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
N_SOURCES = 60            # circuity sampling: random sources ...
TARGETS_PER_SOURCE = 5    # ... times random targets each
MIN_STRAIGHT_M = 500      # ignore very short pairs (circuity is noisy there)
N_ITERATIONS = 2          # number of visualised origin-destination iterations
VIS_DISTANCE_M = 5000     # aim for OD pairs about this far apart in the maps
OUT_DIR = "spatial_plots"

# ---- display settings ----
# For a zoomable window, uncomment the next two lines BEFORE importing pyplot
# (needs:  pip install pyqt5):
# import matplotlib
# matplotlib.use("QtAgg")
LINE_W = 0.7           # road line width in the maps
DPI = 300              # saved PNG resolution (zoom into the saved file)
ZOOM = None            # e.g. (80.20, 80.30, 12.95, 13.10) = (lon_min, lon_max, lat_min, lat_max)

rng = np.random.default_rng(SEED)
os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# 1. LOAD NETWORK
# ----------------------------------------------------------------------------
def load_city(place):
    Gm = ox.graph_from_place(place, network_type="drive")
    U = nx.Graph()
    for u, v, d in Gm.edges(data=True):
        if u == v:
            continue
        hw = d.get("highway", "other")
        hw = hw[0] if isinstance(hw, list) else hw
        L = d.get("length", 1.0)
        if U.has_edge(u, v):
            U[u][v]["length"] = min(U[u][v]["length"], L)
        else:
            U.add_edge(u, v, length=L, highway=hw)
    for n, d in Gm.nodes(data=True):
        if n in U:
            U.nodes[n]["x"], U.nodes[n]["y"] = d["x"], d["y"]
    U = U.subgraph(max(nx.connected_components(U), key=len)).copy()
    return nx.convert_node_labels_to_integers(U)


def area_km2(place):
    gdf = ox.geocode_to_gdf(place)
    return gdf.to_crs(gdf.estimate_utm_crs()).geometry.area.sum() / 1e6


# ----------------------------------------------------------------------------
# 2. GEOMETRY HELPERS
# ----------------------------------------------------------------------------
def straight_m(G, a, b):
    """Straight-line distance in metres between two nodes (equirectangular)."""
    xa, ya, xb, yb = G.nodes[a]["x"], G.nodes[a]["y"], G.nodes[b]["x"], G.nodes[b]["y"]
    dx = (xb - xa) * 111320 * np.cos(np.radians((ya + yb) / 2))
    dy = (yb - ya) * 110540
    return float(np.hypot(dx, dy))


def edge_bearings(G):
    """Bearing (0-180 deg, undirected) and length of every edge."""
    b, w = [], []
    for u, v, d in G.edges(data=True):
        dx = (G.nodes[v]["x"] - G.nodes[u]["x"]) * np.cos(np.radians(G.nodes[u]["y"]))
        dy = G.nodes[v]["y"] - G.nodes[u]["y"]
        b.append(np.degrees(np.arctan2(dx, dy)) % 180)
        w.append(d["length"])
    return np.array(b), np.array(w)


def orientation_entropy(bearings, weights, nbins=36):
    """Shannon entropy of the street-bearing histogram (same idea as osmnx:
    both directions counted, 10 degree bins centred on north)."""
    both = np.concatenate([bearings, (bearings + 180) % 360])
    wts = np.concatenate([weights, weights])
    shifted = (both + 180 / nbins) % 360
    counts, edges = np.histogram(shifted, bins=nbins, range=(0, 360), weights=wts)
    p = counts / counts.sum()
    H = -np.sum(p[p > 0] * np.log(p[p > 0]))
    return H, counts, edges


# ----------------------------------------------------------------------------
# 3. CIRCUITY (sampled)
# ----------------------------------------------------------------------------
def circuity_sample(G):
    nodes = np.array(G.nodes)
    vals = []
    for s in rng.choice(nodes, min(N_SOURCES, len(nodes)), replace=False):
        dist = nx.single_source_dijkstra_path_length(G, s, weight="length")
        for t in rng.choice(nodes, TARGETS_PER_SOURCE, replace=False):
            st = straight_m(G, s, t)
            if t != s and st >= MIN_STRAIGHT_M:
                vals.append(dist[t] / st)
    return np.array(vals)


def pick_od_pair(G):
    """Random source; target whose straight distance is closest to VIS_DISTANCE_M."""
    nodes = np.array(G.nodes)
    s = int(rng.choice(nodes))
    cands = rng.choice(nodes, min(400, len(nodes)), replace=False)
    t = int(min(cands, key=lambda c: abs(straight_m(G, s, int(c)) - VIS_DISTANCE_M)))
    return s, t


# ----------------------------------------------------------------------------
# 4. MAPS
# ----------------------------------------------------------------------------
def segments(G):
    return [[(G.nodes[u]["x"], G.nodes[u]["y"]), (G.nodes[v]["x"], G.nodes[v]["y"])]
            for u, v in G.edges]


def set_geo_aspect(ax, G):
    ys = np.array([G.nodes[n]["y"] for n in G.nodes])
    ax.set_aspect(1 / np.cos(np.radians(ys.mean())))
    ax.autoscale()
    if ZOOM:
        ax.set_xlim(ZOOM[0], ZOOM[1]); ax.set_ylim(ZOOM[2], ZOOM[3])
    ax.set_xticks([]); ax.set_yticks([])


def map_orientation(ax, G):
    b, _ = edge_bearings(G)
    lc = LineCollection(segments(G), array=b, cmap="hsv", linewidths=LINE_W)
    lc.set_clim(0, 180)
    ax.add_collection(lc)
    set_geo_aspect(ax, G)
    ax.set_title("Street orientation (colour = bearing, 0-180 deg)", fontsize=10)


def polar_orientation(ax, G):
    b, w = edge_bearings(G)
    H, counts, edges = orientation_entropy(b, w)
    theta = np.radians((edges[:-1] + edges[1:]) / 2 - 5)
    ax.bar(theta, counts / counts.sum(), width=np.radians(10), color="#7B3FE4", alpha=0.85)
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.set_title(f"Orientation entropy H = {H:.2f}\n(grid ~ 1.4, uniform = {np.log(36):.2f})",
                 fontsize=10)
    return H


def map_od_iteration(ax, G, segs, it):
    s, t = pick_od_pair(G)
    path = nx.shortest_path(G, s, t, weight="length")
    net = nx.shortest_path_length(G, s, t, weight="length")
    st = straight_m(G, s, t)
    ax.add_collection(LineCollection(segs, colors="#bbbbbb", linewidths=LINE_W * 0.6))
    ax.plot([G.nodes[n]["x"] for n in path], [G.nodes[n]["y"] for n in path],
            color="red", lw=2, label=f"network route {net/1000:.1f} km")
    ax.plot([G.nodes[s]["x"], G.nodes[t]["x"]], [G.nodes[s]["y"], G.nodes[t]["y"]],
            "k--", lw=1.5, label=f"straight line {st/1000:.1f} km")
    ax.scatter(G.nodes[s]["x"], G.nodes[s]["y"], marker="*", s=200, c="gold", edgecolor="k", zorder=5)
    ax.scatter(G.nodes[t]["x"], G.nodes[t]["y"], marker="X", s=100, c="gold", edgecolor="k", zorder=5)
    set_geo_aspect(ax, G)
    ax.set_title(f"Iteration {it}: circuity = {net/st:.2f}", fontsize=10)
    ax.legend(loc="lower left", fontsize=8)


# ----------------------------------------------------------------------------
# 5. MAIN
# ----------------------------------------------------------------------------
def analyse(city, place):
    G = load_city(place)
    n, m = G.number_of_nodes(), G.number_of_edges()
    A = area_km2(place)
    total_km = sum(d["length"] for _, _, d in G.edges(data=True)) / 1000
    inter = sum(1 for _, k in G.degree if k >= 3)

    b, w = edge_bearings(G)
    H, _, _ = orientation_entropy(b, w)
    circ = circuity_sample(G)

    cls = {}
    for _, _, d in G.edges(data=True):
        cls[d["highway"]] = cls.get(d["highway"], 0) + d["length"]
    cls = {k: v / 1000 for k, v in sorted(cls.items(), key=lambda kv: -kv[1])}

    # ---- figure: orientation map + polar plot + two OD iterations ----
    fig = plt.figure(figsize=(16, 15))
    ax1 = fig.add_subplot(2, 2, 1); map_orientation(ax1, G)
    ax2 = fig.add_subplot(2, 2, 2, projection="polar"); polar_orientation(ax2, G)
    segs = segments(G)
    for it in range(1, N_ITERATIONS + 1):
        map_od_iteration(fig.add_subplot(2, 2, 2 + it), G, segs, it)
    fig.suptitle(f"{city}: spatial analysis", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"{city}_spatial.png"), dpi=DPI)

    return dict(N=n, area=A, inter_density=inter / A, street_density=total_km / A,
                H=H, circ_mean=circ.mean(), circ_median=np.median(circ), circ=circ,
                classes=cls, total_km=total_km)


def main():
    R = {c: analyse(c, p) for c, p in CITIES.items()}

    # ---- circuity distribution ----
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for c, col in zip(R, ("#FF7A1A", "#0B8F62")):
        ax.hist(R[c]["circ"], bins=30, alpha=0.7, color=col, density=True,
                label=f"{c} (mean {R[c]['circ_mean']:.2f})")
    ax.set_xlabel("circuity = route distance / straight-line distance")
    ax.set_ylabel("density"); ax.legend(frameon=False)
    ax.set_title("Circuity of sampled origin-destination pairs", fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "circuity_comparison.png"), dpi=200)

    # ---- summary ----
    print("\n" + "=" * 62)
    print(f"{'':<34}" + "".join(f"{c:>14}" for c in CITIES))
    print("-" * 62)
    for label, key, fmt in [("area (km2)", "area", "{:,.0f}"),
                            ("road length (km)", "total_km", "{:,.0f}"),
                            ("intersections per km2", "inter_density", "{:.1f}"),
                            ("street length per km2 (km)", "street_density", "{:.1f}"),
                            ("orientation entropy H", "H", "{:.2f}"),
                            ("mean circuity", "circ_mean", "{:.3f}"),
                            ("median circuity", "circ_median", "{:.3f}")]:
        print(f"{label:<34}" + "".join(f"{fmt.format(R[c][key]):>14}" for c in CITIES))
    print("=" * 62)

    print("\nRoad-class share of total length (%):")
    classes = sorted({k for c in R for k in R[c]["classes"]},
                     key=lambda k: -sum(R[c]["classes"].get(k, 0) for c in R))[:8]
    for k in classes:
        print(f"  {k:<16}" + "".join(
            f"{100 * R[c]['classes'].get(k, 0) / R[c]['total_km']:>10.1f}" for c in CITIES))

    print(f"\nPlots saved in ./{OUT_DIR}/")
    plt.show()


if __name__ == "__main__":
    main()