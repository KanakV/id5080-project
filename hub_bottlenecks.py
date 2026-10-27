"""
Bottlenecks between named hubs of a city, three ways per hub pair:

  1. Effective resistance R_ij (resistance distance): the network as a resistor grid.
     Low R = many parallel routes; high R = traffic squeezes through few roads.
     Unit conductance per segment (R_unit) and conductance = 1/length (R_len, in km).
     Computed exactly: grounded Laplacian, one sparse LU, one solve per hub.
  2. Spectral bottleneck of the corridor between the two hubs (nodes inside the ellipse
     |p-a| + |p-b| <= 1.3 d_ab + 2 km): lambda_2 of the normalized Laplacian and the
     Fiedler sweep-cut conductance phi; 'splits' = the cut separates the two hubs, i.e.
     the corridor's worst bottleneck lies between them.
  3. Max-flow min-cut between the hubs' 1 km neighbourhoods (unit capacity per segment):
     the fewest road segments whose loss disconnects one hub area from the other.

Hub coordinates are approximate centres, snapped to the nearest junction; hubs more than
1.5 km from the network (outside the municipal boundary) are dropped.

Usage:  python hub_bottlenecks.py [path/to/graph.graphml]
Writes results/<tag>_hub_bottlenecks.csv, results/<tag>_hub_bottlenecks_report.txt,
figures/<tag>_hub_bottlenecks_dark.png
"""
import csv
import itertools
import os
import sys
import time

import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as sp
from matplotlib.collections import LineCollection
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

from graph_utils import classify_cuts, load_simple, log, set_min_cut, sweep_cut, to_km

sys.stdout.reconfigure(encoding="utf-8")
GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
tag = os.path.splitext(os.path.basename(GRAPH))[0]
CITY = tag.split("_")[0].capitalize()
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

HUBS = {  # (lat, lon), approximate centres
    "bengaluru": {"MG Road": (12.9756, 77.6066), "Whitefield": (12.9698, 77.7500),
                  "Electronic City": (12.8452, 77.6602), "Koramangala": (12.9352, 77.6245),
                  "Hebbal": (13.0358, 77.5970), "Peenya": (13.0285, 77.5197),
                  "Yeshwanthpur": (13.0237, 77.5503), "Jayanagar": (12.9250, 77.5938)},
    "chennai": {"T Nagar": (13.0418, 80.2341), "Guindy": (13.0067, 80.2206),
                "Sholinganallur": (12.9010, 80.2279), "Anna Nagar": (13.0850, 80.2101),
                "Ambattur": (13.1143, 80.1548), "Tambaram": (12.9249, 80.1000),
                "Chennai Central": (13.0827, 80.2757), "Velachery": (12.9815, 80.2180)},
    "barcelona": {"Pl. Catalunya": (41.3870, 2.1700), "Sagrada Familia": (41.4036, 2.1744),
                  "Barceloneta": (41.3809, 2.1897), "Sants": (41.3790, 2.1400),
                  "Diagonal Mar": (41.4030, 2.2040), "Camp Nou": (41.3809, 2.1228),
                  "Gracia": (41.4026, 2.1564), "Sant Andreu": (41.4354, 2.1897)},
}
SNAP_MAX_KM, NBHD_KM = 1.5, 1.0
SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
ORANGE = "#ff8a3d"

S = load_simple(GRAPH)
nodes = list(S.nodes)
idx = {u: i for i, u in enumerate(nodes)}
lon = np.array([S.nodes[u]["x"] for u in nodes])
lat = np.array([S.nodes[u]["y"] for u in nodes])
lon0, lat0 = lon.mean(), lat.mean()
P = to_km(lon, lat, lon0, lat0)
tree = cKDTree(P)

hubs = []
for name, (la, lo) in HUBS[tag.split("_")[0]].items():
    q = to_km([lo], [la], lon0, lat0)[0]
    d, i = tree.query(q)
    if d > SNAP_MAX_KM:
        log(f"dropping {name}: nearest junction {d:.1f} km away (outside the city graph)")
        continue
    hubs.append((name, int(i)))
    log(f"hub {name:16s} -> node {nodes[i]} ({d * 1000:.0f} m away)")
H = len(hubs)
hix = np.array([i for _, i in hubs])
# roads leaving each hub's 1 km area: a pair's min cut can never exceed the smaller of the two
area = [{nodes[i] for i in tree.query_ball_point(P[h], NBHD_KM)} for _, h in hubs]
area_out = [sum(1 for _ in nx.edge_boundary(S, a_)) for a_ in area]
for (n_, _), k_ in zip(hubs, area_out):
    log(f"  {n_:16s}: {k_} road segments leave its {NBHD_KM:.0f} km area")

# ---------------------------------------------------------------- 1. effective resistance
def resistance_matrix(weight):
    A = nx.to_scipy_sparse_array(S, nodelist=nodes, weight=weight, format="csr", dtype=float)
    if weight == "length":
        A.data = 1.0 / np.maximum(A.data / 1000.0, 1e-3)    # conductance = 1 / length (km)
    L = laplacian(A).tocsc()
    g = next(i for i in range(len(nodes)) if i not in set(hix))   # ground a non-hub node
    keep = np.ones(len(nodes), bool)
    keep[g] = False
    Lg = L[keep][:, keep].tocsc()
    pos_in_g = np.cumsum(keep) - 1
    t0 = time.time()
    lu = splu(Lg)
    B = np.zeros((Lg.shape[0], H))
    B[pos_in_g[hix], np.arange(H)] = 1.0
    X = lu.solve(B)[pos_in_g[hix]]                    # H x H block of Lg^-1
    log(f"effective resistance ({weight or 'unit'}) via LU [{time.time() - t0:.0f}s]")
    d = np.diag(X)
    return d[:, None] + d[None, :] - 2 * X


R_unit = resistance_matrix(None)
R_len = resistance_matrix("length")

# ---------------------------------------------------------------- 2 + 3. corridors and min cuts
rows = []
for a, b in itertools.combinations(range(H), 2):
    (na, ia), (nb, ib) = hubs[a], hubs[b]
    pa, pb = P[ia], P[ib]
    dab = float(np.hypot(*(pa - pb)))
    ell = np.hypot(*(P - pa).T) + np.hypot(*(P - pb).T) <= 1.3 * dab + 2.0
    C = S.subgraph([nodes[i] for i in np.flatnonzero(ell)])
    C = C.subgraph(max(nx.connected_components(C), key=len)).copy()
    cn = list(C.nodes)
    Ac = nx.to_scipy_sparse_array(C, nodelist=cn, weight=None, format="csr", dtype=float)
    sw = sweep_cut(Ac)
    ci = {u: i for i, u in enumerate(cn)}
    splits = (nodes[ia] in ci and nodes[ib] in ci and sw["side"][ci[nodes[ia]]] != sw["side"][ci[nodes[ib]]])

    src, snk = area[a], area[b]
    common = src & snk
    t0 = time.time()
    cut, cut_edges = set_min_cut(S, src - common, snk - common)
    limit = (f"{na} access roads" if cut == area_out[a] else f"{nb} access roads" if cut == area_out[b] else "")
    rows.append(dict(hub_a=na, hub_b=nb, dist_km=dab, R_unit=R_unit[a, b], R_len_km=R_len[a, b],
                     corridor_n=len(cn), lam2=sw["lam2"], n_lam2=len(cn) * sw["lam2"], phi=sw["phi_bal"],
                     sweep_cut_edges=sw["n_cut"], splits=bool(splits), min_cut=cut, cut_limited_by=limit,
                     src_n=len(src), snk_n=len(snk)))
    log(f"{na} - {nb}: d={dab:.1f} km  R={R_unit[a, b]:.2f}  phi={sw['phi_bal']:.2e}  "
        f"min cut={cut}  [{time.time() - t0:.0f}s]")

# Where does each min cut sit? If a hub's cut is the same number to (nearly) every other hub,
# the choke point is in that hub's own catchment, not in the corridor between the pair.
where = classify_cuts([(r_["hub_a"], r_["hub_b"], r_["min_cut"]) for r_ in rows])
for r_ in rows:
    r_["cut_limited_by"] = r_["cut_limited_by"] or where[(r_["hub_a"], r_["hub_b"])]

with open(f"results/{tag}_hub_bottlenecks.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

# ---------------------------------------------------------------- report
rep = (f"================= Hub-to-hub bottlenecks: {tag} =================\n"
       f"Graph: giant component, simple undirected, n = {len(nodes):,}, m = {S.number_of_edges():,}\n"
       f"Hubs ({H}): " + ", ".join(n for n, _ in hubs) + "\n\n"
       f"Road segments leaving each hub's {NBHD_KM:.0f} km area (an upper limit on every min cut involving it):\n  "
       + ", ".join(f"{n_} {k_}" for (n_, _), k_ in zip(hubs, area_out)) + "\n\n"
       f"{'pair':36s}{'km':>6s}{'R_unit':>8s}{'R_len':>8s}{'n·λ₂':>8s}{'φ':>10s}{'splits':>7s}{'min cut':>8s}"
       f"  limited by\n")
for r_ in sorted(rows, key=lambda r_: r_["min_cut"]):
    rep += (f"{r_['hub_a'] + ' – ' + r_['hub_b']:36s}{r_['dist_km']:6.1f}{r_['R_unit']:8.2f}{r_['R_len_km']:8.2f}"
            f"{r_['n_lam2']:8.2f}{r_['phi']:10.2e}{'yes' if r_['splits'] else 'no':>7s}{r_['min_cut']:8d}"
            f"  {r_['cut_limited_by']}\n")
Ru = np.array([r_["R_unit"] for r_ in rows])
dk = np.array([r_["dist_km"] for r_ in rows])
mc = np.array([r_["min_cut"] for r_ in rows])
slope, icpt = np.polyfit(np.log(dk), Ru, 1)
rep += (f"\nR_unit vs distance: R ≈ {icpt:.2f} + {slope:.2f}·ln(d km)  (a 2-D grid gives R ∝ ln d)\n"
        f"Median R_unit {np.median(Ru):.2f}, median min cut {np.median(mc):.0f} segments, "
        f"median corridor φ {np.median([r_['phi'] for r_ in rows]):.2e}\n"
        f"Pairs whose corridor bottleneck lies between the hubs (splits = yes): "
        f"{sum(r_['splits'] for r_ in rows)}/{len(rows)}\n"
        f"Where the min cut sits: in one hub's own catchment (same cut to nearly every other hub) "
        f"{sum(r_['cut_limited_by'].startswith('near ') for r_ in rows)}/{len(rows)}, "
        f"in the corridor between the hubs {sum(r_['cut_limited_by'] == 'between hubs' for r_ in rows)}/{len(rows)}, "
        f"at a hub's 1 km access roads {sum(r_['cut_limited_by'].endswith('access roads') for r_ in rows)}/{len(rows)}\n"
        f"Choked hubs: " + (", ".join(sorted({r_['cut_limited_by'][5:] for r_ in rows
                                                if r_['cut_limited_by'].startswith('near ')})) or "none") + "\n")
resid = Ru - (icpt + slope * np.log(dk))
worst = np.argsort(-resid)[:3]
rep += "Most constricted pairs for their distance (largest R above the trend):\n"
for i in worst:
    rep += f"  {rows[i]['hub_a']} – {rows[i]['hub_b']}: R = {Ru[i]:.2f} (+{resid[i]:.2f} over trend), min cut {mc[i]}\n"
print(rep)
with open(f"results/{tag}_hub_bottlenecks_report.txt", "w", encoding="utf-8") as f:
    f.write(rep)

# ---------------------------------------------------------------- figure
log("plotting ...")
names = [n for n, _ in hubs]
fig = plt.figure(figsize=(21, 13), dpi=160)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.6], hspace=0.35, wspace=0.3)
SEQ = mpl.colors.LinearSegmentedColormap.from_list("seq", ["#184f95", "#3987e5", "#9ec5f4"])
SEQ_R = SEQ.reversed()


def heat(ax, M, title, fmt, cmap):
    ax.set_facecolor(SURFACE)
    Mm = np.ma.masked_where(np.eye(H, dtype=bool), M)
    im = ax.imshow(Mm, cmap=cmap)
    for i in range(H):
        for j in range(H):
            if i != j:
                ax.text(j, i, format(M[i, j], fmt), ha="center", va="center", fontsize=8, color=INKD,
                        path_effects=[mpl.patheffects.withStroke(linewidth=2, foreground="#1a1a19")])
    ax.set_xticks(range(H))
    ax.set_yticks(range(H))
    ax.set_xticklabels(names, rotation=45, ha="right", color=INKD, fontsize=9)
    ax.set_yticklabels(names, color=INKD, fontsize=9)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.set_title(title, color=INKD, fontsize=12)
    cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cb.ax.tick_params(colors=SUBD, labelsize=8)


def mat(key):
    M = np.zeros((H, H))
    for r_ in rows:
        a, b = names.index(r_["hub_a"]), names.index(r_["hub_b"])
        M[a, b] = M[b, a] = r_[key]
    return M


heat(fig.add_subplot(gs[0, 0]), R_unit, "Effective resistance R (unit conductance)\nhigher = fewer parallel routes",
     ".2f", SEQ_R)
heat(fig.add_subplot(gs[0, 1]), mat("min_cut"), "Min cut between 1 km hub areas\n(road segments; lower = tighter bottleneck)",
     ".0f", SEQ)
heat(fig.add_subplot(gs[1, 0]), 1e3 * mat("phi"), "Corridor Fiedler conductance φ × 10³\n(lower = sharper bottleneck)",
     ".2f", SEQ)
heat(fig.add_subplot(gs[1, 1]), mat("n_lam2"), "Corridor n·λ₂ (size-normalised)\n(lower = weaker connectivity)", ".1f", SEQ)

ax = fig.add_subplot(gs[:, 2])
ax.set_facecolor(SURFACE)
segs = np.array([[P[idx[u]], P[idx[v]]] for u, v in S.edges])
ax.add_collection(LineCollection(segs, colors="#2f2f2c", linewidths=0.2, zorder=1))
norm = mpl.colors.Normalize(Ru.min(), Ru.max())
for r_ in rows:
    pa, pb = P[hubs[names.index(r_["hub_a"])][1]], P[hubs[names.index(r_["hub_b"])][1]]
    ax.plot(*zip(pa, pb), color=SEQ_R(norm(r_["R_unit"])), lw=1 + 3 * norm(r_["R_unit"]), zorder=2, alpha=0.9)
halo = [mpl.patheffects.withStroke(linewidth=3, foreground=SURFACE)]
for n_, i in hubs:
    ax.scatter(*P[i], s=90, color=ORANGE, zorder=4, edgecolors=SURFACE, linewidths=1.5)
    ax.annotate(n_, P[i], xytext=(7, 7), textcoords="offset points", color=INKD, fontsize=10, fontweight="bold",
                zorder=5, path_effects=halo)
ax.set_aspect("equal")
ax.autoscale()
ax.set_xticks([])
ax.set_yticks([])
for s_ in ax.spines.values():
    s_.set_color(SPINE)
ax.margins(0.08)                                      # room for hub labels at the edges
sm = mpl.cm.ScalarMappable(norm=norm, cmap=SEQ_R)
cb = fig.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.6, pad=0.03)
cb.set_label("effective resistance R between the pair (thicker = more constricted)", color=SUBD)
cb.ax.tick_params(colors=SUBD)
ax.set_title("Hub pairs coloured by resistance", color=INKD, fontsize=12)
fig.suptitle(f"Hub-to-hub bottlenecks — {CITY}", color=INKD, fontsize=15)
fig.savefig(f"figures/{tag}_hub_bottlenecks_dark.png", bbox_inches="tight", facecolor=SURFACE)
log(f"done: results/{tag}_hub_bottlenecks.csv, _report.txt, figures/{tag}_hub_bottlenecks_dark.png")
