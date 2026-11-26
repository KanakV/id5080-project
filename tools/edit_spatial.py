"""
Edit distance to a robust street network: how many NEW LOCAL ROADS must be added to a city
centre before it reaches a fixed, absolute level of robustness?  Full write-up: EDIT_DISTANCE.md

1. Window   the ~9,000 junctions nearest the CBD (same size for every city).
2. Metric   Schneider robustness R = mean size of the largest connected piece while junctions
            are removed one by one (max 0.5), under
              random failure  (3 fixed random orders)                Albert et al. 2000
              degree attack   (most-connected junction first,        Holme et al. 2002
                               recomputed after every removal)
3. Target   the same absolute finish line for every city: a % of the R of a perfect square
            grid with the same number of junctions (the classic planned street layout).
4. Edits    only buildable roads: a new road may join two NEIGHBOURING junctions (Delaunay
            neighbours <= 300 m apart, i.e. the next junction in some direction) where both
            have fewer than 4 roads. Each road counts 1, so edit distance = roads added.
5. Choice   greedy: at each step add the roads that keep the most junctions connected during
            the collapse (see gain()); recompute R exactly after every batch.
            Baseline: the same number of roads chosen at random from the same candidates.

Usage:  python tools/edit_spatial.py [graphml ...]   (default: Bengaluru, Chennai)
Writes results/edit_spatial_report.txt, results/edit_spatial.csv, figures/edit_spatial_dark.png
"""
import csv
import os
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import LineCollection
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Delaunay

from graph_utils import SEED, load_simple, log, robustness, to_km

sys.stdout.reconfigure(encoding="utf-8")
WINDOW_N = 9000               # junctions in the central window
LOCAL_KM = 0.3                # longest new road
MAX_DEG = 4                   # new roads only at junctions with fewer than 4 roads
N_SNAP = 10                   # collapse snapshots per removal order (gain estimate)
CHECKPOINTS = 18
LEVELS = (0.6, 0.7, 0.8)      # finish lines, as a fraction of the grid's R
MODES = {"rand": "random failure", "deg": "degree attack"}
REF = "data/barcelona_city_barcelona_drive.graphml"   # reference line only, never a target
SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
CITY_COL = {"Bengaluru": "#3987e5", "Chennai": "#c98500"}
NEW_COL = "#ff8a3d"

paths = sys.argv[1:] or ["data/bengaluru_city_bengaluru_drive.graphml", "data/chennai_city_chennai_drive.graphml"]
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)


def candidates(G, nodes, xy):
    """All buildable new roads: Delaunay neighbours <= LOCAL_KM apart, both ends degree < MAX_DEG,
    not already joined. Returns index arrays (a, b)."""
    tri = Delaunay(xy).simplices
    p = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]]), axis=1), axis=0)
    p = p[np.hypot(*(xy[p[:, 0]] - xy[p[:, 1]]).T) <= LOCAL_KM]
    deg = np.array([G.degree(u) for u in nodes])
    p = p[(deg[p[:, 0]] < MAX_DEG) & (deg[p[:, 1]] < MAX_DEG)]
    p = np.array([q for q in p if not G.has_edge(nodes[q[0]], nodes[q[1]])], dtype=np.int64).reshape(-1, 2)
    return p[:, 0], p[:, 1]


def gain(a, b, info, n):
    """Estimated rise in R_rand + R_deg if road a-b is added, and the piece it reconnects.

    R is the area under the 'largest piece vs junctions removed' curve. Take N_SNAP snapshots
    of each collapse; at each, label the pieces that remain. If a and b both still stand in
    different pieces, the road joins them and the largest piece grows by
    max(size_a + size_b - largest, 0) (the smaller piece, when one side is the largest).
    Adding that up over the snapshots (times their spacing) estimates the change in the area."""
    A = info["A"]
    g = np.zeros(len(a))
    best = np.zeros(len(a))
    piece = np.arange(len(a), dtype=np.int64) + (1 << 40)
    snap = 0
    for m, orders in info["orders"].items():
        for o, S in zip(orders, info["curves"][m]):
            below = np.flatnonzero(S[:n] < 0.05)                 # snapshots cover the collapse only
            k_end = int(below[0]) if len(below) else n
            rank = np.empty(n, dtype=np.int64)
            rank[o] = np.arange(n)
            for k in np.unique(np.linspace(k_end / N_SNAP, k_end, N_SNAP).astype(int)):
                alive = np.flatnonzero(rank >= k)
                _, lab_alive = connected_components(A[alive][:, alive], directed=False)
                lab = np.full(n, -1)
                lab[alive] = lab_alive
                size = np.bincount(lab_alive)
                la, lb = lab[a], lab[b]
                ok = np.flatnonzero((la >= 0) & (lb >= 0) & (la != lb))
                grow = np.maximum(size[la[ok]] + size[lb[ok]] - size.max(), 0)
                term = (k_end / N_SNAP / n) * grow / n / len(orders)
                g[ok] += term
                upd = term > best[ok]
                small = np.where(size[la[ok]] <= size[lb[ok]], la[ok], lb[ok])
                piece[ok[upd]] = snap * (1 << 24) + small[upd]
                best[ok[upd]] = term[upd]
                snap += 1
    return g, piece


def greedy_batch(a, b, g, piece, k):
    """Up to k roads with the highest gain; one per junction and one per reconnected piece per
    batch, so a batch does not spend several roads rescuing the same cul-de-sac. Once no road
    has a positive gain, the rest are taken in random order (so greedy can use every road)."""
    out, used, pieces = [], set(), set()
    for i in np.lexsort((np.random.default_rng(SEED).random(len(g)), -g)):
        if len(out) == k:
            break
        if a[i] in used or b[i] in used:
            continue
        if g[i] > 0:
            if piece[i] in pieces:
                continue
            pieces.add(piece[i])
        out.append((int(a[i]), int(b[i])))
        used.update((a[i], b[i]))
    return out


def run(W, xy, rule):
    """Add roads by `rule` ('greedy' or 'random'); exact R after every checkpoint."""
    G = W.copy()
    nodes = list(G.nodes)
    n = len(nodes)
    a0, b0 = candidates(G, nodes, xy)
    marks = np.unique(np.round(np.concatenate([np.geomspace(10, len(a0), CHECKPOINTS),       # fine at the start
                                               np.linspace(0, len(a0), 21)[1:]])).astype(int))   # every 5%
    R, info = robustness(G)
    hist = [dict(rule=rule, added=0, **{f"R_{m}": R[m] for m in MODES})]
    built = []
    shuffled = list(np.random.default_rng(SEED).permutation(len(a0)))
    for mk in marks:
        need = mk - len(built)
        if rule == "greedy":
            a, b = candidates(G, nodes, xy)
            g, piece = gain(a, b, info, n)
            batch = greedy_batch(a, b, g, piece, need)
        else:
            batch = []
            while shuffled and len(batch) < need:
                i = shuffled.pop()
                u, v = nodes[a0[i]], nodes[b0[i]]
                if G.degree(u) < MAX_DEG and G.degree(v) < MAX_DEG and not G.has_edge(u, v):
                    batch.append((int(a0[i]), int(b0[i])))
                    G.add_edge(u, v)                  # so later degree checks see it
        if not batch:
            break
        for u, v in batch:
            G.add_edge(nodes[u], nodes[v])
        built += batch
        R, info = robustness(G)
        hist.append(dict(rule=rule, added=len(built), **{f"R_{m}": R[m] for m in MODES}))
    return hist, built


def roads_to_reach(hist, ceil, level):
    """New roads at the first checkpoint where BOTH modes are >= level x ceiling (None if never)."""
    for h in hist:
        if all(h[f"R_{m}"] >= level * ceil[m] for m in MODES):
            return h["added"]
    return None


# ---------------------------------------------------------------- reference ceiling + Barcelona
side = int(round(np.sqrt(WINDOW_N)))
ceil, _ = robustness(nx.grid_2d_graph(side, side))
log(f"ceiling ({side}x{side} grid): " + "  ".join(f"{m}={ceil[m]:.3f}" for m in MODES))
B = load_simple(REF)
R_bcn, _ = robustness(B)
del B

res, rows = {}, []
for p in paths:
    tag = os.path.splitext(os.path.basename(p))[0]
    city = tag.split("_")[0].capitalize()
    S = load_simple(p)
    lon = np.array([d["x"] for _, d in S.nodes(data=True)])
    lat = np.array([d["y"] for _, d in S.nodes(data=True)])
    f = f"results/{tag}_assortativity.npz"            # CBD = centroid of the top-betweenness junctions
    cbd = np.load(f, allow_pickle=True)["cbd"] if os.path.exists(f) else np.array([lon.mean(), lat.mean()])
    P = to_km(lon, lat, *cbd)
    dist = np.hypot(*P.T)
    r_win = float(np.sort(dist)[min(WINDOW_N, len(dist)) - 1])
    W = S.subgraph([u for u, d in zip(S.nodes, dist) if d <= r_win])
    W = W.subgraph(max(nx.connected_components(W), key=len)).copy()
    pos = dict(zip(S.nodes, map(tuple, P)))
    del S
    xy = np.array([pos[u] for u in W.nodes])
    m0 = W.number_of_edges()
    log(f"{city}: window {r_win:.2f} km, n={W.number_of_nodes():,}, m={m0:,}")
    out = {}
    for rule in ("greedy", "random"):
        hist, built = run(W, xy, rule)
        out[rule] = (hist, built)
        log(f"  {rule:6s}: " + " | ".join(f"+{h['added']} " + " ".join(f"{m}={h['R_' + m] / ceil[m]:.0%}" for m in MODES)
                                          for h in hist[::4] + hist[-1:]))
        rows += [dict(city=city, **h) for h in hist]
    res[city] = dict(m0=m0, r_win=r_win, W=W, pos=pos, xy=xy, runs=out)

with open("results/edit_spatial.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

# ---------------------------------------------------------------- report
pc = lambda R, m: f"{R:.3f} ({100 * R / ceil[m]:3.0f}%)"
rep = ("================= Edit distance to a robust street network =================\n"
       f"Window: {WINDOW_N:,} junctions nearest the CBD.  New roads: neighbouring junctions (Delaunay, <= "
       f"{1000 * LOCAL_KM:.0f} m),\nboth ends with < {MAX_DEG} roads.  R = Schneider robustness (max 0.5).  "
       f"Ceiling = R of a {side}x{side} square grid.\n\n--- Starting robustness: R (% of ceiling) ---\n"
       f"{'':24s}{'random failure':>18s}{'degree attack':>18s}\n"
       f"{'ceiling (grid)':24s}" + "".join(f"{ceil[m]:>18.3f}" for m in MODES) + "\n")
for city, r in res.items():
    h0 = r["runs"]["greedy"][0][0]
    rep += f"{city + ' window':24s}" + "".join(f"{pc(h0['R_' + m], m):>18s}" for m in MODES) + "\n"
rep += (f"{'Barcelona (whole city)':24s}" + "".join(f"{pc(R_bcn[m], m):>18s}" for m in MODES) + "   (reference only)\n"
        "\n--- Edit distance: new roads per 1000 existing roads until BOTH modes reach the level ---\n"
        f"{'level':>8s}" + "".join(f"{c + ' ' + rule:>22s}" for c in res for rule in ("greedy", "random")) + "\n")
for L in LEVELS:
    cells = []
    for city, r in res.items():
        for rule in ("greedy", "random"):
            k = roads_to_reach(r["runs"][rule][0], ceil, L)
            cells.append("not reached" if k is None else f"{1000 * k / r['m0']:.0f}")
    rep += f"{100 * L:>7.0f}%" + "".join(f"{c:>22s}" for c in cells) + "\n"
rep += "\n--- Highest level reached when every buildable road is used ---\n"
for city, r in res.items():
    for rule in ("greedy", "random"):
        h = r["runs"][rule][0][-1]
        rep += (f"  {city:10s} {rule:6s}: {h['added']:,} roads ({1000 * h['added'] / r['m0']:.0f} per 1000)  "
                + "  ".join(f"{MODES[m]} {100 * h['R_' + m] / ceil[m]:.0f}%" for m in MODES) + "\n")
print(rep)
with open("results/edit_spatial_report.txt", "w", encoding="utf-8") as fh:
    fh.write(rep)

# ---------------------------------------------------------------- figure
cities = list(res)
reached = [L for L in LEVELS if all(roads_to_reach(res[c]["runs"]["greedy"][0], ceil, L) is not None for c in cities)]
HEAD = max(reached) if reached else LEVELS[0]
fig = plt.figure(figsize=(18, 9.5), dpi=160)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(2, len(cities) + 1, width_ratios=[1.25] + [1] * len(cities), hspace=0.38, wspace=0.12)


def style(ax, title):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=10)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.set_title(title, color=INKD, fontsize=13, loc="left")


ax = fig.add_subplot(gs[0, 0])
style(ax, "A  Edit distance: new local roads per 1000 existing roads\n    until random failure AND degree attack "
          "reach the level")
x = np.arange(len(LEVELS))
wbar = 0.8 / len(cities)
for i, city in enumerate(cities):
    r = res[city]
    vals = [roads_to_reach(r["runs"]["greedy"][0], ceil, L) for L in LEVELS]
    for j, v in enumerate(vals):
        xv = x[j] + (i - (len(cities) - 1) / 2) * wbar
        if v is None:
            ax.text(xv, 5, "not\nreached", ha="center", va="bottom", color=SUBD, fontsize=9)
            continue
        v = 1000 * v / r["m0"]
        ax.bar(xv, v, width=wbar * 0.9, color=CITY_COL[city], label=city if j == 0 else None)
        ax.text(xv, v, f"{v:.0f}", ha="center", va="bottom", color=INKD, fontsize=10)
ax.set_xticks(x)
ax.set_xticklabels([f"{100 * L:.0f}% of grid" for L in LEVELS], color=INKD, fontsize=11)
ax.grid(True, axis="y", color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.legend(frameon=False, labelcolor=INKD, fontsize=11, loc="upper left")

ax = fig.add_subplot(gs[1, 0])
style(ax, "B  Robustness as roads are added\n    (weaker of the two modes)")
for city in cities:
    r = res[city]
    for rule, ls in (("greedy", "-"), ("random", "--")):
        hist = r["runs"][rule][0]
        ax.plot([1000 * h["added"] / r["m0"] for h in hist],
                [100 * min(h[f"R_{m}"] / ceil[m] for m in MODES) for h in hist], ls, color=CITY_COL[city], lw=2,
                marker="o" if rule == "greedy" else None, ms=3,
                label=f"{city}, {'best roads first' if rule == 'greedy' else 'random roads'}")
for L in LEVELS:
    ax.axhline(100 * L, color=GRID, lw=1)
bw = 100 * min(R_bcn[m] / ceil[m] for m in MODES)
ax.axhline(bw, color="#d55181", ls=":", lw=1.2)
ax.text(1, bw + 1, "Barcelona today (reference)", color="#d55181", fontsize=9)
ax.set_xscale("symlog", linthresh=10)
ax.set_xlim(left=0)
ax.set_xlabel("new roads per 1000 existing roads", color=SUBD)
ax.set_ylabel("% of grid robustness", color=SUBD)
ax.grid(True, axis="x", color=GRID, lw=0.6)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9, loc="upper left")

for j, city in enumerate(cities):
    ax = fig.add_subplot(gs[:, j + 1])
    r = res[city]
    W, pos = r["W"], r["pos"]
    k = roads_to_reach(r["runs"]["greedy"][0], ceil, HEAD)
    built = r["runs"]["greedy"][1][:k]
    nodes = list(W.nodes)
    style(ax, f"{city}: the {len(built):,} new roads to reach {100 * HEAD:.0f}%")
    ax.add_collection(LineCollection([[pos[u], pos[v]] for u, v in W.edges], colors="#4a4a46", linewidths=0.35))
    ax.add_collection(LineCollection([[pos[nodes[u]], pos[nodes[v]]] for u, v in built], colors=NEW_COL,
                                     linewidths=1.3, zorder=3))
    ax.set_aspect("equal")
    ax.autoscale()
    ax.set_xticks([])
    ax.set_yticks([])

fig.suptitle("How many new local roads make the city centre robust?  Same absolute target for both cities: "
             "a % of a perfect grid's robustness\n(new road = between neighbouring junctions, ≤ 300 m;  robustness = "
             "random failure and targeted attack on the most-connected junctions)", color=INKD, fontsize=13)
fig.savefig("figures/edit_spatial_dark.png", facecolor=SURFACE, bbox_inches="tight")
log("wrote results/edit_spatial_report.txt, results/edit_spatial.csv, figures/edit_spatial_dark.png")
