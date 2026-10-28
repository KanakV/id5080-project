# Author: Omkar
# =============================================================================
# PEAK-HOUR TRAFFIC MODEL, BUSY-ROAD REMOVAL AND BOTTLENECK (NARROW-ROAD)
# ANALYSIS FOR THE CHENNAI AND BENGALURU ROAD NETWORKS
#
# PURPOSE
#   1. Build a simple traffic model on the OSM road network: off-peak vs
#      evening peak (e.g. Chennai 6-8 pm), with roads slowing down as they fill up.
#   2. "Strip off" (close) the busiest roads and measure how much trips slow
#      down and how many trips are cut off completely.
#   3. Find bottleneck roads: narrow roads (few lanes) that carry heavy peak
#      flow, and test what happens if they are widened.
#
# REASON
#   My earlier robustness experiments removed nodes at random or by betweenness
#   and measured the giant component. That treats every road as equal and
#   ignores traffic. In a real city, losing a busy road hurts even if the
#   network stays connected, because trips get much longer. This script adds
#   that missing piece: travel time under congestion.
#
# EXPECTATION (hypotheses to check, not results)
#   1. Peak hour makes trips noticeably slower than off-peak, and a small
#      share of roads is over capacity (flow/capacity > 1).
#   2. Closing the busiest roads hurts far more than closing the same number
#      of random roads: trip time rises steeply while the network is still
#      almost fully connected. Random closures should be mild.
#   3. Closing the most overloaded roads (highest flow/capacity) should be
#      even worse, because traffic spills onto neighbouring roads and
#      overloads them too.
#   4. Many overloaded roads are narrow (1-2 lanes). They are bottlenecks
#      because central, shortest-path routes pass through low-capacity roads.
#   5. Widening the few worst bottlenecks should give a clear drop in peak
#      trip time. The city that gains more has the more bottleneck-limited network.
#
# MODEL AND ASSUMPTIONS (state these in the report)
#   - OSM has NO traffic counts. Demand is simulated: random origin-destination
#     (OD) pairs routed along shortest-time paths. Flow on a road = number of
#     trips using it, which is similar to edge betweenness.
#   - Speed and lanes come from OSM tags (maxspeed, lanes) when present, else
#     defaults by road class (DEFAULT_SPEED, DEFAULT_LANES below).
#   - Capacity = lanes x CAP_PER_LANE vehicles/hour. The network is treated as
#     undirected (both directions share one flow and one capacity).
#   - Congested time uses the BPR formula t = t0 * (1 + 0.15 (flow/cap)^4).
#   - Traffic is assigned with the method of successive averages, so trips
#     re-route around congestion.
#   - Demand scale is relative: peak demand is set so that, if everyone used
#     free-flow shortest paths, the busiest 5% of roads would be at
#     TARGET_VC times capacity. Off-peak = OFFPEAK_FACTOR x peak demand.
#     So results compare scenarios and cities; absolute vehicle counts are
#     not real. For real numbers, replace this with measured data (Google /
#     TomTom travel times, or local traffic counts).
#   - Road closures are EDGE removals (closing road segments).
#
# Run:   pip install osmnx networkx numpy scipy matplotlib
#        python road_traffic_peak.py
# =============================================================================

import os
import re
import time
import numpy as np
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

# ----------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------
CITIES = {
    "Chennai":   "Chennai, Tamil Nadu, India",
    "Bengaluru": "Bengaluru, Karnataka, India",
}
SEED = 42
N_ORIGINS = 60              # random origins ...
DESTS_PER_ORIGIN = 15       # ... each with this many random destinations
MSA_ITERS = 3               # traffic assignment iterations (more = smoother, slower)
TARGET_VC = 3.0             # peak intensity (see header); raise for a harsher rush hour
OFFPEAK_FACTOR = 0.25       # off-peak demand as a fraction of peak
CAP_PER_LANE = 900          # vehicles / hour / lane
VC_CAP = 4.0                # cap on flow/capacity inside BPR (avoids blow-ups)
NARROW_LANES = 1            # roads with <= this many lanes count as "narrow"
BUSY_TOP_PCT = 95           # "busy" = flow above this percentile
REMOVE_FRACS = [0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10]   # fraction of roads closed
VIS_FRACS = (0.01, 0.05)    # the two "iterations" drawn on maps
RANDOM_RUNS = 2             # random-closure curves are averaged over this many runs
WIDEN_PCT = 0.01            # widen (x2 lanes) this share of roads (worst narrow bottlenecks)

# display
# For a zoomable window uncomment BEFORE plotting (pip install pyqt5):
# import matplotlib
# matplotlib.use("QtAgg")
LINE_W = 0.7
DPI = 300
ZOOM = None                 # e.g. (80.20, 80.30, 12.95, 13.10) = (lon_min, lon_max, lat_min, lat_max)
OUT_DIR = "traffic_plots"

DEFAULT_SPEED = {"motorway": 60, "trunk": 50, "primary": 40, "secondary": 35,
                 "tertiary": 30, "residential": 25, "unclassified": 25,
                 "living_street": 15, "service": 15}
DEFAULT_LANES = {"motorway": 3, "trunk": 3, "primary": 3, "secondary": 2,
                 "tertiary": 2, "residential": 1, "unclassified": 1,
                 "living_street": 1, "service": 1}

C_CITY = {"Chennai": "#FF7A1A", "Bengaluru": "#0B8F62"}
rng = np.random.default_rng(SEED)
os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# 1. LOAD AND BUILD THE NETWORK
# ----------------------------------------------------------------------------
def first_num(val, default=None):
    if isinstance(val, list):
        val = val[0] if val else None
    if val is None:
        return default
    m = re.search(r"\d+(\.\d+)?", str(val))
    if not m:
        return default
    x = float(m.group())
    return x * 1.609 if "mph" in str(val) else x


def load_city(place):
    Gm = ox.graph_from_place(place, network_type="drive")
    U = nx.Graph()
    for u, v, d in Gm.edges(data=True):
        if u == v or U.has_edge(u, v):
            continue
        hw = d.get("highway", "unclassified")
        hw = hw[0] if isinstance(hw, list) else hw
        hw = hw.replace("_link", "")
        U.add_edge(u, v, length=max(float(d.get("length", 1.0)), 1.0), highway=hw,
                   lanes=first_num(d.get("lanes")), maxspeed=first_num(d.get("maxspeed")))
    for n, d in Gm.nodes(data=True):
        if n in U:
            U.nodes[n]["x"], U.nodes[n]["y"] = d["x"], d["y"]
    U = U.subgraph(max(nx.connected_components(U), key=len)).copy()
    return nx.convert_node_labels_to_integers(U)


class Net:
    """Plain arrays describing the road network and the OD demand."""


def build_net(G):
    net = Net()
    net.n = G.number_of_nodes()
    edges = list(G.edges(data=True))
    net.m = len(edges)
    net.U = np.array([u for u, _, _ in edges])
    net.V = np.array([v for _, v, _ in edges])
    length = np.array([d["length"] for _, _, d in edges])
    hw = [d["highway"] for _, _, d in edges]
    speed = np.array([d["maxspeed"] or DEFAULT_SPEED.get(h, 25) for h, (_, _, d) in zip(hw, edges)])
    lanes = np.array([d["lanes"] or DEFAULT_LANES.get(h, 1) for h, (_, _, d) in zip(hw, edges)])
    net.t0 = length / (np.clip(speed, 5, None) / 3.6) / 60.0        # free-flow time, minutes
    net.cap = np.clip(lanes, 1, None) * float(CAP_PER_LANE)
    net.lanes = lanes
    net.hw = np.array(hw)
    net.length = length
    net.X = np.array([G.nodes[i]["x"] for i in range(net.n)])
    net.Y = np.array([G.nodes[i]["y"] for i in range(net.n)])
    net.eid = {}
    for k, (u, v) in enumerate(zip(net.U, net.V)):
        net.eid[int(u) * net.n + int(v)] = k
        net.eid[int(v) * net.n + int(u)] = k
    net.segs = np.stack((np.column_stack((net.X[net.U], net.Y[net.U])),
                         np.column_stack((net.X[net.V], net.Y[net.V]))), axis=1)
    # OD demand: same pairs for every scenario
    net.origins = rng.choice(net.n, min(N_ORIGINS, net.n), replace=False)
    net.dests = []
    for o in net.origins:
        pool = np.setdiff1d(rng.choice(net.n, DESTS_PER_ORIGIN + 5, replace=False), [o])
        net.dests.append(pool[:DESTS_PER_ORIGIN])
    return net


# ----------------------------------------------------------------------------
# 2. TRAFFIC MODEL
# ----------------------------------------------------------------------------
def bpr(t0, flow, cap):
    return t0 * (1.0 + 0.15 * np.minimum(flow / cap, VC_CAP) ** 4)


def route(net, t, active, w):
    """All-or-nothing: send w trips along each OD shortest-time path.
    Returns (edge flow, travel time of every OD pair; inf = disconnected)."""
    idx = np.flatnonzero(active)
    A = csr_matrix((t[idx], (net.U[idx], net.V[idx])), shape=(net.n, net.n))
    dist, pred = dijkstra(A, directed=False, indices=net.origins, return_predecessors=True)
    flow = np.zeros(net.m)
    od_t, k = [], 0
    for i, o in enumerate(net.origins):
        for d in net.dests[i]:
            dd = dist[i, d]
            od_t.append(dd)
            if not np.isfinite(dd):
                continue
            v = d
            while v != o:
                p = pred[i, v]
                flow[net.eid[int(p) * net.n + int(v)]] += w
                v = p
    return flow, np.array(od_t)


def assign(net, active, w, cap):
    """Method of successive averages. Returns (flow, congested times, OD times)."""
    flow = np.zeros(net.m)
    for k in range(1, MSA_ITERS + 1):
        aon, _ = route(net, bpr(net.t0, flow, cap), active, w)
        flow += (aon - flow) / k
    t = bpr(net.t0, flow, cap)
    _, od_t = route(net, t, active, w)
    return flow, t, od_t


def calibrate(net):
    """Peak demand per trip so that the busiest 5% of roads reach TARGET_VC x capacity
    under free-flow shortest paths."""
    f1, _ = route(net, net.t0, np.ones(net.m, bool), 1.0)
    vc = f1 / net.cap
    return TARGET_VC / np.percentile(vc[f1 > 0], 95)


# ----------------------------------------------------------------------------
# 3. CLOSURE EXPERIMENT
# ----------------------------------------------------------------------------
def closure_curve(net, order, w, base_mean_t):
    """For each closure fraction: relative mean trip time, % OD pairs cut off."""
    ratio, cut = [], []
    for f in REMOVE_FRACS:
        active = np.ones(net.m, bool)
        active[order[:int(round(f * net.m))]] = False
        _, _, od_t = assign(net, active, w, net.cap)
        ok = np.isfinite(od_t)
        ratio.append(od_t[ok].mean() / base_mean_t if ok.any() else np.nan)
        cut.append(100 * (1 - ok.mean()))
    return np.array(ratio), np.array(cut)


# ----------------------------------------------------------------------------
# 4. PLOTTING HELPERS
# ----------------------------------------------------------------------------
def prep_ax(ax, net, title):
    ax.set_aspect(1 / np.cos(np.radians(net.Y.mean())))
    ax.set_xlim(net.X.min(), net.X.max()); ax.set_ylim(net.Y.min(), net.Y.max())
    if ZOOM:
        ax.set_xlim(ZOOM[0], ZOOM[1]); ax.set_ylim(ZOOM[2], ZOOM[3])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10)


def draw_values(ax, net, values, vmin, vmax, cmap, lw):
    order = np.argsort(values)                       # draw the highest values on top
    lc = LineCollection(net.segs[order], array=values[order], cmap=cmap,
                        norm=Normalize(vmin, vmax), linewidths=lw[order])
    ax.add_collection(lc)
    return lc


def line_widths(flow):
    return LINE_W * (0.4 + 2.2 * flow / max(flow.max(), 1e-9))


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25)


# ----------------------------------------------------------------------------
# 5. ONE CITY
# ----------------------------------------------------------------------------
def analyse(city, G):
    t_start = time.time()
    net = build_net(G)
    print(f"\n=== {city}: {net.n:,} nodes, {net.m:,} road segments, "
          f"{len(net.origins) * DESTS_PER_ORIGIN} OD pairs ===")
    allon = np.ones(net.m, bool)

    # ---- peak vs off-peak ----
    scale = calibrate(net)
    w_peak, w_off = scale, scale * OFFPEAK_FACTOR
    f_off, t_off, od_off = assign(net, allon, w_off, net.cap)
    f_pk, t_pk, od_pk = assign(net, allon, w_peak, net.cap)
    vc_off, vc_pk = f_off / net.cap, f_pk / net.cap
    slowdown = od_pk.mean() / od_off.mean()
    over = vc_pk > 1
    narrow = net.lanes <= NARROW_LANES
    print(f"  mean trip time: off-peak {od_off.mean():.1f} min, peak {od_pk.mean():.1f} min "
          f"(x{slowdown:.2f})")
    print(f"  roads over capacity at peak: {100 * over.mean():.1f}% "
          f"| of those, narrow (<= {NARROW_LANES} lane): {100 * narrow[over].mean() if over.any() else 0:.0f}%")

    # ---- closure experiments (peak demand) ----
    base = od_pk.mean()
    order_busy = np.argsort(-f_pk)
    order_vc = np.argsort(-vc_pk)
    curves = {"Busiest roads (by flow)": closure_curve(net, order_busy, w_peak, base),
              "Most overloaded (flow/capacity)": closure_curve(net, order_vc, w_peak, base)}
    rr, rc = [], []
    for _ in range(RANDOM_RUNS):
        r, c = closure_curve(net, rng.permutation(net.m), w_peak, base)
        rr.append(r); rc.append(c)
    curves["Random roads"] = (np.mean(rr, axis=0), np.mean(rc, axis=0))

    # ---- bottlenecks: narrow & overloaded; widening what-if ----
    k = max(1, int(WIDEN_PCT * net.m))
    cand = np.flatnonzero(narrow)
    bott = cand[np.argsort(-vc_pk[cand])[:k]] if len(cand) else np.array([], int)
    cap2 = net.cap.copy(); cap2[bott] *= 2
    _, _, od_w = assign(net, allon, w_peak, cap2)
    gain = 100 * (1 - od_w.mean() / od_pk.mean())
    busy = f_pk >= np.percentile(f_pk[f_pk > 0], BUSY_TOP_PCT)
    print(f"  narrow AND busy roads: {int((narrow & busy).sum())} "
          f"({100 * (narrow & busy).sum() / max(busy.sum(), 1):.0f}% of the busy roads)")
    print(f"  widening the {len(bott)} worst narrow bottlenecks x2 cuts peak trip time by {gain:.1f}%")

    # ---- figure 1: load maps ----
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.5))
    lw = line_widths(f_pk)
    lc = draw_values(axes[0], net, vc_off, 0, 1.5, "RdYlGn_r", line_widths(f_off))
    prep_ax(axes[0], net, f"Off-peak: flow / capacity (mean trip {od_off.mean():.0f} min)")
    draw_values(axes[1], net, vc_pk, 0, 1.5, "RdYlGn_r", lw)
    prep_ax(axes[1], net, f"Evening peak: flow / capacity (mean trip {od_pk.mean():.0f} min)")
    fig.colorbar(lc, ax=axes[:2], shrink=0.6, label="flow / capacity (red = over capacity)")
    axes[2].add_collection(LineCollection(net.segs, colors="#cccccc", linewidths=LINE_W * 0.5))
    nb = np.flatnonzero(narrow & busy)
    axes[2].add_collection(LineCollection(net.segs[nb], colors="#d62728", linewidths=LINE_W * 2.2))
    prep_ax(axes[2], net, f"Bottlenecks: narrow roads (<= {NARROW_LANES} lane) among the "
                          f"top {100 - BUSY_TOP_PCT}% busiest ({len(nb)} roads, red)")
    fig.suptitle(f"{city}: modelled traffic load", fontsize=14, fontweight="bold")
    fig.savefig(os.path.join(OUT_DIR, f"{city}_load_maps.png"), dpi=DPI, bbox_inches="tight")

    # ---- figure 2: the closure process, two iterations on the map ----
    fig, axes = plt.subplots(1, len(VIS_FRACS), figsize=(8.5 * len(VIS_FRACS), 8))
    for it, (ax, f) in enumerate(zip(np.atleast_1d(axes), VIS_FRACS), 1):
        kk = int(round(f * net.m))
        ax.add_collection(LineCollection(net.segs, colors="#cccccc", linewidths=LINE_W * 0.5))
        ax.add_collection(LineCollection(net.segs[order_busy[:kk]], colors="#d62728",
                                         linewidths=LINE_W * 2.0))
        j = REMOVE_FRACS.index(f) if f in REMOVE_FRACS else None
        r = curves["Busiest roads (by flow)"]
        txt = f"trip time x{r[0][j]:.2f}, {r[1][j]:.1f}% trips cut off" if j is not None else ""
        prep_ax(ax, net, f"Iteration {it}: close the busiest {100 * f:.1f}% of roads "
                         f"({kk:,} segments, red)\n{txt}")
    fig.suptitle(f"{city}: busy-road closure", fontsize=14, fontweight="bold")
    fig.savefig(os.path.join(OUT_DIR, f"{city}_closure_maps.png"), dpi=DPI, bbox_inches="tight")

    # ---- figure 3: scatter + closure curves ----
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2))
    ax = axes[0]
    jitter = 1 + rng.uniform(-0.12, 0.12, net.m)
    ax.scatter(net.cap * jitter, f_pk, s=3, c=np.where(over, "#d62728", "#7f7f7f"), alpha=0.5)
    mx = max(net.cap.max(), f_pk.max())
    ax.plot([0, mx], [0, mx], "k--", lw=1, label="flow = capacity")
    ax.set_xlabel("road capacity (vehicles / hour, from lanes)")
    ax.set_ylabel("peak flow (modelled vehicles / hour)")
    ax.set_title("Capacity vs peak flow (red = over capacity)", fontweight="bold")
    ax.legend(frameon=False); clean(ax)
    cols = {"Busiest roads (by flow)": "#d62728", "Most overloaded (flow/capacity)": "#7B3FE4",
            "Random roads": "#2F6BFF"}
    for name, (r, c) in curves.items():
        axes[1].plot(100 * np.array(REMOVE_FRACS), r, "o-", color=cols[name], lw=2.2, label=name)
        axes[2].plot(100 * np.array(REMOVE_FRACS), c, "o-", color=cols[name], lw=2.2, label=name)
    axes[1].set_ylabel("mean peak trip time (relative to no closure)")
    axes[2].set_ylabel("% of OD pairs cut off")
    for a, ttl in zip(axes[1:], ("Trip-time increase", "Trips cut off")):
        a.set_xlabel("% of roads closed"); a.set_title(ttl, fontweight="bold")
        a.legend(frameon=False); clean(a)
    fig.suptitle(f"{city}: bottlenecks and closure effect", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"{city}_analysis.png"), dpi=200)

    print(f"  done in {time.time() - t_start:.0f}s")
    return dict(slowdown=slowdown, over=100 * over.mean(),
                narrow_over=100 * narrow[over].mean() if over.any() else 0.0,
                gain=gain, curves=curves, t_off=od_off.mean(), t_pk=od_pk.mean())


# ----------------------------------------------------------------------------
# 6. MAIN
# ----------------------------------------------------------------------------
def main():
    R = {}
    for city, place in CITIES.items():
        R[city] = analyse(city, load_city(place))

    # ---- comparison figure ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for city in R:
        axes[0].plot(100 * np.array(REMOVE_FRACS), R[city]["curves"]["Busiest roads (by flow)"][0],
                     "o-", lw=2.5, color=C_CITY[city], label=f"{city}: busiest roads")
        axes[0].plot(100 * np.array(REMOVE_FRACS), R[city]["curves"]["Random roads"][0],
                     "o--", lw=1.8, color=C_CITY[city], alpha=0.6, label=f"{city}: random roads")
    axes[0].set_xlabel("% of roads closed"); axes[0].set_ylabel("mean peak trip time (relative)")
    axes[0].set_title("Closure effect: Chennai vs Bengaluru", fontweight="bold")
    axes[0].legend(frameon=False); clean(axes[0])
    names = list(R)
    axes[1].bar(names, [R[c]["slowdown"] for c in names], color=[C_CITY[c] for c in names])
    axes[1].set_ylabel("peak trip time / off-peak trip time")
    axes[1].set_title("Peak-hour slowdown", fontweight="bold"); clean(axes[1])
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "chennai_vs_bengaluru_traffic.png"), dpi=200)

    # ---- summary ----
    j5 = REMOVE_FRACS.index(0.05)
    print("\n" + "=" * 70)
    print(f"{'':<46}" + "".join(f"{c:>12}" for c in R))
    print("-" * 70)
    rows = [("peak slowdown (peak / off-peak trip time)", lambda r: r["slowdown"], "{:.2f}"),
            ("roads over capacity at peak (%)", lambda r: r["over"], "{:.1f}"),
            ("  of which narrow roads (%)", lambda r: r["narrow_over"], "{:.0f}"),
            ("trip time x, 5% busiest roads closed",
             lambda r: r["curves"]["Busiest roads (by flow)"][0][j5], "{:.2f}"),
            ("trip time x, 5% overloaded roads closed",
             lambda r: r["curves"]["Most overloaded (flow/capacity)"][0][j5], "{:.2f}"),
            ("trip time x, 5% random roads closed",
             lambda r: r["curves"]["Random roads"][0][j5], "{:.2f}"),
            ("trips cut off, 5% busiest closed (%)",
             lambda r: r["curves"]["Busiest roads (by flow)"][1][j5], "{:.1f}"),
            ("peak trip time gain from widening bottlenecks (%)", lambda r: r["gain"], "{:.1f}")]
    for label, fn, fmt in rows:
        print(f"{label:<46}" + "".join(f"{fmt.format(fn(R[c])):>12}" for c in R))
    print("=" * 70)
    a, b = list(R)
    ta = R[a]["curves"]["Busiest roads (by flow)"][0][j5]
    tb = R[b]["curves"]["Busiest roads (by flow)"][0][j5]
    print(f"\nMore robust to closure of the busiest roads (smaller trip-time increase): "
          f"{a if ta < tb else b}")
    print(f"Plots saved in ./{OUT_DIR}/")
    plt.show()


if __name__ == "__main__":
    main()
