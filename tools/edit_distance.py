"""
Edit distance from a city street network to four target networks, measured on
the robustness / flow-capacity metrics that make a network "good":

  barcelona : the real Barcelona drive network (Cerda grid; reference optimum)
  er        : Erdos-Renyi G(n, m) on the city's own node positions
  sf        : scale-free (Chung-Lu, gamma = 2.5, same mean degree) on the same positions
  ringspoke : concentric rings joined by radial spokes (robust to node failure)

Exact graph edit distance is NP-hard at this size, so two bounds are reported:

  LOWER bound (exact): edges that MUST be added / removed just to turn the city's
    degree sequence into the target's (sorted matching minimises the L1 gap):
      adds >= sum max(t_i - d_i, 0) / 2,   removes >= sum max(d_i - t_i, 0) / 2
  UPPER bound (constructive): greedy edits on a central window of the city (same n
    as Barcelona) until the metric vector is within 10% of the target's gap:
      barcelona : add short links (<= 300 m) with the worst detour, dead ends first;
                  remove links only at nodes of degree >= 5
      ringspoke : add ring links (every 1 km around the CBD) and radial spoke links
      er        : random rewiring (remove a random edge, add a random pair)
      sf        : preferential rewiring (remove a random edge, attach prop. to degree)
    Adds and removes are counted separately; metrics are logged at ~15 checkpoints.

Metrics (size-intensive): Schneider robustness R under random failure and adaptive
degree attack, bridge share, n*lambda_2, N-S/E-W min cut / sqrt(n), length-weighted
efficiency <d_euclid / d_network>, median detour, assortativity r, mean degree.

Usage:  python edit_distance.py [graphml] [--targets barcelona,er,sf,ringspoke]
                                [--ref data/barcelona_city_barcelona_drive.graphml] [--window-n N]
Writes results/<tag>_edit_distance_report.txt, results/<tag>_edit_distance.npz,
figures/<tag>_edit_distance_dark.png, figures/<tag>_edit_map_dark.png
"""
import argparse
import csv
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.collections import LineCollection
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

from graph_utils import (SEED, adaptive_degree_order, band_min_cut, giant_curve, lambda2, load_simple, log,
                         to_km)

ap = argparse.ArgumentParser()
ap.add_argument("graph", nargs="?", default="data/bengaluru_city_bengaluru_drive.graphml")
ap.add_argument("--targets", default="barcelona,er,sf,ringspoke")
ap.add_argument("--ref", default="data/barcelona_city_barcelona_drive.graphml")
ap.add_argument("--window-n", type=int, default=0, help="window size (default: Barcelona's n)")
ap.add_argument("--checkpoints", type=int, default=14)
args = ap.parse_args()
sys.stdout.reconfigure(encoding="utf-8")

tag = os.path.splitext(os.path.basename(args.graph))[0]
CITY = tag.split("_")[0].capitalize()
TARGETS = args.targets.split(",")
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
COLORS = {"city": "#3987e5", "barcelona": "#d55181", "er": "#d95926", "sf": "#9085e9", "ringspoke": "#008300"}
MARKERS = {"barcelona": "D", "er": "s", "sf": "v", "ringspoke": "o"}
LABELS = {"barcelona": "Barcelona", "er": "Erdős–Rényi", "sf": "scale-free", "ringspoke": "ring–spoke"}

# robustness + flow capacity: these define "reached the target"
RF = ["rob_rand", "rob_attack", "bridges", "n_lam2", "mincut", "efficiency"]
BETTER = {"rob_rand": 1, "rob_attack": 1, "bridges": -1, "n_lam2": 1, "mincut": 1, "efficiency": 1}
# "optimal" targets: matching OR beating them counts; reference models: match them (two-sided)
ONE_SIDED = {"barcelona", "ringspoke"}
ALL = RF + ["detour", "r", "k_mean"]
MNAME = {"rob_rand": "robustness R, random failure", "rob_attack": "robustness R, degree attack",
         "bridges": "bridge share (lower = more loops)", "n_lam2": "n·λ₂ (spectral connectivity)",
         "mincut": "min cut / √n (flow capacity)", "efficiency": "efficiency ⟨d_eucl / d_net⟩",
         "detour": "median detour d_net / d_eucl", "r": "assortativity r", "k_mean": "mean degree ⟨k⟩"}


# ---------------------------------------------------------------- metrics
def metrics(G, pos, n_src=24, planar_like=True):
    """G: simple undirected nx.Graph with edge 'length' (m); pos: node -> (x, y) km.
    planar_like=False for randomised graphs (LU fill-in makes shift-invert too costly)."""
    nodes = list(G.nodes)
    n_all = len(nodes)
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, format="csr")
    rng = np.random.default_rng(SEED)
    rob_rand = np.mean([giant_curve(A.indptr, A.indices, rng.permutation(n_all))[:n_all].mean() for _ in range(3)])
    rob_attack = giant_curve(A.indptr, A.indices, adaptive_degree_order(A.indptr, A.indices))[:n_all].mean()
    k_mean = 2 * G.number_of_edges() / n_all

    H = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    hn = list(H.nodes)
    n, m = len(hn), H.number_of_edges()
    Ah = nx.to_scipy_sparse_array(H, nodelist=hn, weight=None, format="csr").astype(float)
    bridges = sum(1 for _ in nx.bridges(H)) / m
    n_lam2 = n * lambda2(Ah, planar_like=planar_like)
    xy = np.array([pos[u] for u in hn])
    mincut = (band_min_cut(H, xy[:, 1])[0] + band_min_cut(H, xy[:, 0])[0]) / 2 / np.sqrt(n)

    Aw = nx.to_scipy_sparse_array(H, nodelist=hn, weight="length", format="csr")
    src = rng.choice(n, size=min(n_src, n), replace=False)
    D = dijkstra(Aw, indices=src)
    E = np.hypot(xy[src, None, 0] - xy[None, :, 0], xy[src, None, 1] - xy[None, :, 1]) * 1000
    ok = E > 300
    efficiency = float(np.mean(E[ok] / D[ok]))
    detour = float(np.median(D[ok] / E[ok]))

    deg = np.asarray(Ah.sum(axis=1)).ravel()
    co = Ah.tocoo()
    ki, kj = deg[co.row], deg[co.col]                      # both directions of every edge
    num = (ki * kj).mean() - ki.mean() ** 2
    den = (ki ** 2).mean() - ki.mean() ** 2
    r = float(num / den) if den > 0 else 0.0
    return dict(rob_rand=rob_rand, rob_attack=rob_attack, bridges=bridges, n_lam2=n_lam2, mincut=mincut,
                efficiency=efficiency, detour=detour, r=r, k_mean=k_mean)


def gap(x, t, x0, one_sided=False, noise=None):
    """Normalised distance over the RF metrics: 1 at the start, 0 at the target.
    one_sided: a metric already better than the target contributes 0.
    noise: per-metric sampling noise of a random target; metrics where the city already
    starts within 2x that noise count as matched and are left out."""
    d = []
    for k in RF:
        if abs(x0[k] - t[k]) <= max(1e-9, 2 * (noise or {}).get(k, 0.0)):
            continue
        short = (t[k] - x[k]) * BETTER[k] if one_sided else abs(x[k] - t[k])
        d.append(max(short, 0.0) / abs(x0[k] - t[k]))
    return float(np.mean(d)) if d else 0.0


def add_edge(G, pos, u, v):
    G.add_edge(u, v, length=1000 * float(np.hypot(*np.subtract(pos[u], pos[v]))))


# ---------------------------------------------------------------- city window
S = load_simple(args.graph)
lon = np.array([d["x"] for _, d in S.nodes(data=True)])
lat = np.array([d["y"] for _, d in S.nodes(data=True)])
cnodes = list(S.nodes)
cbd = None
for f in (f"results/{tag}_assortativity.npz",):
    if os.path.exists(f):
        cbd = np.load(f, allow_pickle=True)["cbd"]
if cbd is None and os.path.exists(f"results/{tag}_top_betweenness.csv"):
    with open(f"results/{tag}_top_betweenness.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))[:8]
    cbd = np.array([np.mean([float(r_["lon"]) for r_ in rows]), np.mean([float(r_["lat"]) for r_ in rows])])
if cbd is None:
    cbd = np.array([lon.mean(), lat.mean()])
P = to_km(lon, lat, *cbd)
pos_city = dict(zip(cnodes, map(tuple, P)))

REF = None
if "barcelona" in TARGETS or not args.window_n:
    REF = load_simple(args.ref)
win_n = args.window_n or REF.number_of_nodes()
is_ref = os.path.abspath(args.graph) == os.path.abspath(args.ref)
if win_n >= len(cnodes) or is_ref:
    W0 = S.copy()
    R_win = float(np.hypot(*P.T).max())
else:
    dist = np.hypot(*P.T)
    R_win = float(np.sort(dist)[win_n - 1])
    keep = [u for u, dd in zip(cnodes, dist) if dd <= R_win]
    W0 = S.subgraph(keep)
    W0 = W0.subgraph(max(nx.connected_components(W0), key=len)).copy()
n0, m0 = W0.number_of_nodes(), W0.number_of_edges()
log(f"{CITY} window: radius {R_win:.2f} km around CBD, n={n0:,}, m={m0:,}")
wnodes = list(W0.nodes)

# ---------------------------------------------------------------- targets
targets, tdeg, floor_sample, noise = {}, {}, {}, {}
rng = np.random.default_rng(SEED)
for T in TARGETS:
    log(f"building target {T} ...")
    if T == "barcelona":
        B = REF
        blon = np.array([d["x"] for _, d in B.nodes(data=True)])
        blat = np.array([d["y"] for _, d in B.nodes(data=True)])
        bP = to_km(blon, blat, blon.mean(), blat.mean())
        Gt, post = B, dict(zip(B.nodes, map(tuple, bP)))
    elif T in ("er", "sf"):
        if T == "er":
            Gs = nx.gnm_random_graph(n0, m0, seed=SEED)
        else:
            w = (np.arange(n0) + 10.0) ** (-1 / (2.5 - 1))
            w *= 2 * m0 / w.sum()
            Gs = nx.Graph(nx.expected_degree_graph(w.tolist(), seed=SEED, selfloops=False))
        perm = rng.permutation(n0)
        if T == "sf":                                   # rewiring samples endpoints by these weights
            sf_weight = np.empty(n0)
            sf_weight[perm] = w
        Gt = nx.relabel_nodes(Gs, {i: wnodes[perm[i]] for i in range(n0)})
        Gt.add_nodes_from(wnodes)
        post = pos_city
        for u, v in Gt.edges:
            Gt[u][v]["length"] = 1000 * float(np.hypot(*np.subtract(post[u], post[v])))
    else:  # ringspoke: R rings x S spokes + hub, roughly square cells, same disc
        Rr = max(2, int(round(np.sqrt(n0 / np.pi))))
        Ss = max(6, int(round(np.pi * Rr)))
        Gt, post = nx.Graph(), {"c": (0.0, 0.0)}
        for i in range(Rr):
            for j in range(Ss):
                rad, th = R_win * (i + 1) / Rr, 2 * np.pi * j / Ss
                post[(i, j)] = (rad * np.cos(th), rad * np.sin(th))
        for i in range(Rr):
            for j in range(Ss):
                Gt.add_edge((i, j), (i, (j + 1) % Ss))                  # ring
                Gt.add_edge((i, j), (i + 1, j) if i + 1 < Rr else (i, j))  # spoke
                if i == 0:
                    Gt.add_edge("c", (0, j))
        Gt.remove_edges_from(nx.selfloop_edges(Gt))
        for u, v in Gt.edges:
            Gt[u][v]["length"] = 1000 * float(np.hypot(*np.subtract(post[u], post[v])))
        log(f"  ring-spoke: {Rr} rings x {Ss} spokes, n={Gt.number_of_nodes():,}, m={Gt.number_of_edges():,}")
    targets[T] = metrics(Gt, post, planar_like=T not in ("er", "sf"))
    tdeg[T] = np.array([d for _, d in Gt.degree()])
    if T in ("er", "sf"):   # more samples of the same model = how close "identical" can get
        samples = []
        for s in (1, 2, 3, 4):
            if T == "er":
                G2 = nx.gnm_random_graph(n0, m0, seed=SEED + s)
            else:
                G2 = nx.Graph(nx.expected_degree_graph(w.tolist(), seed=SEED + s, selfloops=False))
            perm2 = rng.permutation(n0)
            G2 = nx.relabel_nodes(G2, {i: wnodes[perm2[i]] for i in range(n0)})
            G2.add_nodes_from(wnodes)
            for u, v in G2.edges:
                G2[u][v]["length"] = 1000 * float(np.hypot(*np.subtract(post[u], post[v])))
            samples.append(metrics(G2, post, planar_like=False))
        noise[T] = {k: max(abs(s_[k] - targets[T][k]) for s_ in samples) for k in ALL}
        floor_sample[T] = samples
    log(f"  {T}: " + "  ".join(f"{k}={v:.4g}" for k, v in targets[T].items()))


# ---------------------------------------------------------------- lower bound (degree sequences)
def degree_lb(d, t_dist):
    """Exact lower bound on adds/removes to turn degree sequence d into one drawn from t_dist."""
    d = np.sort(d)[::-1].astype(float)
    t = np.sort(np.quantile(t_dist, (np.arange(len(d)) + 0.5) / len(d), method="inverted_cdf"))[::-1]
    return int(np.ceil(np.clip(t - d, 0, None).sum() / 2)), int(np.ceil(np.clip(d - t, 0, None).sum() / 2))


deg_full = np.array([dd for _, dd in S.degree()])
deg_win = np.array([dd for _, dd in W0.degree()])
lb = {T: dict(full=degree_lb(deg_full, tdeg[T]), win=degree_lb(deg_win, tdeg[T])) for T in TARGETS}


# ---------------------------------------------------------------- move sets (upper bound)
def spatial_adds(G, pos, b, maxd_km=0.30, max_deg=4):
    """Short links with the worst current detour; dead-end closures weighted x2."""
    nodes = list(G.nodes)
    ix = {u: i for i, u in enumerate(nodes)}
    xy = np.array([pos[u] for u in nodes])
    deg = np.array([G.degree(u) for u in nodes])
    pairs = cKDTree(xy).query_pairs(maxd_km, output_type="ndarray")
    pairs = pairs[(deg[pairs[:, 0]] < max_deg) & (deg[pairs[:, 1]] < max_deg)]
    pairs = np.array([p for p in pairs if not G.has_edge(nodes[p[0]], nodes[p[1]])])
    if len(pairs) == 0:
        return []
    Aw = nx.to_scipy_sparse_array(G, nodelist=nodes, weight="length", format="csr")
    srcs = np.unique(pairs[:, 0])
    net = np.empty(len(pairs))
    where = {s: np.flatnonzero(pairs[:, 0] == s) for s in srcs}
    for k in range(0, len(srcs), 400):
        chunk = srcs[k:k + 400]
        Dm = dijkstra(Aw, indices=chunk, limit=3000.0)
        for row, s in enumerate(chunk):
            net[where[s]] = Dm[row, pairs[where[s], 1]]
    euc = 1000 * np.hypot(*(xy[pairs[:, 0]] - xy[pairs[:, 1]]).T)
    score = np.minimum(net, 3000.0) / np.maximum(euc, 1.0)
    score *= np.where((deg[pairs[:, 0]] == 1) | (deg[pairs[:, 1]] == 1), 2.0, 1.0)
    out, used = [], set()
    for p in pairs[np.argsort(-score)]:
        u, v = nodes[p[0]], nodes[p[1]]
        if u in used or v in used:
            continue
        out.append((u, v))
        used.update((u, v))
        if len(out) == b:
            break
    return out


def ring_spoke_candidates(G, pos, ring_km=0.5, n_spokes=24, tol=0.12, step_km=0.25, max_gap=0.5):
    """Ring links every ring_km around the CBD and radial links along n_spokes rays.
    Each ring (spoke) keeps ONE node per step_km of arc (of length) — the one nearest
    the ideal curve — so consecutive picks trace a single clean ring, not a zig-zag
    across the tolerance band. Ordered inner rings / inner spoke segments first."""
    nodes = list(G.nodes)
    xy = np.array([pos[u] for u in nodes])
    rad, th = np.hypot(*xy.T), np.arctan2(xy[:, 1], xy[:, 0])

    def chain(on, coord, off, key, closed):
        if len(on) < 2:
            return
        bins = np.floor(coord[on] / step_km).astype(int)
        pick = {}
        for i, b in zip(on, bins):
            if b not in pick or off[i] < off[pick[b]]:
                pick[b] = i
        seq = [pick[b] for b in sorted(pick)]
        pairs = list(zip(seq, seq[1:] + seq[:1])) if closed and len(seq) > 2 else list(zip(seq[:-1], seq[1:]))
        for a, b in pairs:
            if np.hypot(*(xy[a] - xy[b])) <= max_gap:
                cands.append((key(a), nodes[a], nodes[b]))

    cands = []
    for k in range(1, int(rad.max() // ring_km) + 1):
        R = k * ring_km
        off = np.abs(rad - R)
        on = np.flatnonzero(off < tol)
        chain(on, (th + np.pi) * R, off, lambda a, k=k: k, closed=True)       # arc length along ring
    for j in range(n_spokes):
        ang = 2 * np.pi * j / n_spokes
        perp = np.abs(-np.sin(ang) * xy[:, 0] + np.cos(ang) * xy[:, 1])
        along = np.cos(ang) * xy[:, 0] + np.sin(ang) * xy[:, 1]
        on = np.flatnonzero((perp < tol) & (along > 0))
        chain(on, along, perp, lambda a, along=along: along[a] / ring_km, closed=False)
    cands.sort(key=lambda c: c[0])
    return [(u, v) for _, u, v in cands if not G.has_edge(u, v)]


ring_spoke_edits = [0]          # edits spent on rings/spokes before densification began


def reach_D(T, x0):
    """Gap that counts as 'reached': 0.10, or 1.5x the gap between two samples of the
    same random model (a fully rewired graph cannot beat the model's own noise)."""
    if T in floor_sample:
        return max(0.10, 1.5 * max(gap(s_, targets[T], x0, noise=noise[T]) for s_ in floor_sample[T]))
    return 0.10


def run_target(T):
    G = W0.copy()
    x0 = metrics(G, pos_city)
    t = targets[T]
    budget = 2 * m0 if T == "er" else m0             # 2m edits = every edge rewired once
    marks = np.unique(np.round(np.geomspace(20, budget, args.checkpoints)).astype(int))
    hist = [dict(edits=0, adds=0, removes=0, D=1.0, **x0)]
    adds = removes = 0
    added, removed = [], []
    rs_queue = ring_spoke_candidates(G, pos_city) if T == "ringspoke" else None
    rr = np.random.default_rng(SEED + 7)
    originals = [el for el in W0.edges]                 # rewire original city edges first
    rr.shuffle(originals)
    exhausted = False
    for mk in marks:
        need = mk - (adds + removes)
        if need <= 0:
            continue
        if T == "er":
            nodes = list(G.nodes)
            for _ in range(min(need // 2, len(originals))):
                u, v = originals.pop()
                if not G.has_edge(u, v):
                    continue
                G.remove_edge(u, v)
                removed.append((u, v))
                removes += 1
                while True:
                    a, b = nodes[rr.integers(len(nodes))], nodes[rr.integers(len(nodes))]
                    if a != b and not G.has_edge(a, b):
                        break
                add_edge(G, pos_city, a, b)
                added.append((a, b))
                adds += 1
        elif T == "ringspoke":
            # phase 1: the ring roads and spokes themselves; phase 2: local loop closures
            take = [(u, v) for u, v in rs_queue[:need] if not G.has_edge(u, v)]
            rs_queue = rs_queue[need:]
            ring_spoke_edits[0] += len(take)
            if len(take) < need:
                take += spatial_adds(G, pos_city, need - len(take))
            for u, v in take:
                if not G.has_edge(u, v):
                    add_edge(G, pos_city, u, v)
            added += take
            adds += len(take)
            exhausted = len(take) < need // 4
        else:  # barcelona
            n_rm = 0
            if (tdeg[T] >= 5).mean() < np.mean([d >= 5 for _, d in G.degree()]):
                hi = [(u, v) for u, v in G.edges if max(G.degree(u), G.degree(v)) >= 5]
                br = {frozenset(e) for e in nx.bridges(G)}
                hi = [e for e in hi if frozenset(e) not in br][:max(1, need // 10)]
                G.remove_edges_from(hi)
                removed += hi
                removes += len(hi)
                n_rm = len(hi)
            take = spatial_adds(G, pos_city, need - n_rm)
            for u, v in take:
                add_edge(G, pos_city, u, v)
            added += take
            adds += len(take)
            exhausted = len(take) < (need - n_rm) // 4
        x = metrics(G, pos_city, planar_like=T != "er" or adds < 0.05 * m0)
        D = gap(x, t, x0, one_sided=T in ONE_SIDED, noise=noise.get(T))
        hist.append(dict(edits=adds + removes, adds=adds, removes=removes, D=D, **x))
        log(f"  {T:9s} edits={adds + removes:6,d} (+{adds:,} / -{removes:,})  D={D:.3f}  "
            + "  ".join(f"{k}={x[k]:.3g}" for k in RF))
        if D <= reach_D(T, x0) or exhausted:
            break
    return hist, added, removed, exhausted


def run_sf():
    """Preferential rewiring: both ends of each new edge drawn with probability proportional
    to the target's Chung-Lu weight, so a full rewire lands on the target model."""
    G = W0.copy()
    x0 = metrics(G, pos_city)
    t = targets["sf"]
    marks = np.unique(np.round(np.geomspace(20, 2 * m0, args.checkpoints)).astype(int))
    hist = [dict(edits=0, adds=0, removes=0, D=1.0, **x0)]
    adds = removes = 0
    added, removed = [], []
    rr = np.random.default_rng(SEED + 11)
    nodes = list(G.nodes)                                # same order as wnodes
    p = sf_weight / sf_weight.sum()
    originals = list(W0.edges)                           # rewire original city edges first
    rr.shuffle(originals)
    for mk in marks:
        need = (mk - (adds + removes)) // 2
        if need <= 0:
            continue
        for _ in range(min(need, len(originals))):
            u, v = originals.pop()
            if not G.has_edge(u, v):
                continue
            G.remove_edge(u, v)
            removed.append((u, v))
            removes += 1
            while True:
                a, b = (nodes[i] for i in rr.choice(len(nodes), size=2, p=p))
                if a != b and not G.has_edge(a, b):
                    break
            add_edge(G, pos_city, a, b)
            added.append((a, b))
            adds += 1
        x = metrics(G, pos_city, planar_like=adds < 0.05 * m0)
        D = gap(x, t, x0, noise=noise["sf"])
        hist.append(dict(edits=adds + removes, adds=adds, removes=removes, D=D, **x))
        log(f"  sf        edits={adds + removes:6,d} (+{adds:,} / -{removes:,})  D={D:.3f}  "
            + "  ".join(f"{k}={x[k]:.3g}" for k in RF))
        if D <= reach_D("sf", x0):
            break
    return hist, added, removed, False


runs = {}
for T in TARGETS:
    log(f"greedy edits toward {T} ...")
    runs[T] = run_sf() if T == "sf" else run_target(T)

# ---------------------------------------------------------------- report
x0 = runs[TARGETS[0]][0][0]
rep = (f"================= Edit distance to target networks: {tag} =================\n"
       f"Window: {R_win:.2f} km around the CBD (centroid of top-betweenness junctions), "
       f"n = {n0:,}, m = {m0:,} (full city n = {S.number_of_nodes():,}, m = {S.number_of_edges():,})\n\n"
       f"--- Metric vectors (window start vs targets at about the same n) ---\n"
       f"{'metric':38s}{CITY + ' window':>18s}" + "".join(f"{LABELS[T]:>14s}" for T in TARGETS) + "\n")
for k in ALL:
    rep += f"{MNAME[k]:38s}{x0[k]:>18.4g}" + "".join(f"{targets[T][k]:>14.4g}" for T in TARGETS) + "\n"

rep += ("\n--- LOWER bound: degree-sequence edit distance (exact; sorted matching) ---\n"
        "Edges that must change for the city to have the target's degree mix at all; a true lower\n"
        "bound on graph edit distance to ANY graph with the target's degree distribution.\n"
        f"{'target':14s}{'full city: adds':>18s}{'removes':>10s}{'per 1000 edges':>16s}"
        f"{'window: adds':>15s}{'removes':>10s}\n")
for T in TARGETS:
    (fa, fr), (wa, wr) = lb[T]["full"], lb[T]["win"]
    rep += (f"{LABELS[T]:14s}{fa:>18,d}{fr:>10,d}{1000 * (fa + fr) / S.number_of_edges():>16.1f}"
            f"{wa:>15,d}{wr:>10,d}\n")

rep += ("\n--- UPPER bound: edits that reach the target's robustness & flow capacity (constructive) ---\n"
        "A different criterion from the degree bound: matching the 6 metrics can need fewer edits than\n"
        "matching the degree mix (e.g. adding loops raises robustness without making the degrees equal).\n"
        "reached = normalised gap D over the 6 robustness/flow metrics <= 0.10 (D = 1 at the start).\n"
        "Barcelona / ring-spoke: one-sided (beating the target on a metric counts as reaching it).\n"
        "ER / scale-free: two-sided (the metric must match the model, better or worse).\n"
        f"{'target':14s}{'reached':>9s}{'edits':>8s}{'adds':>8s}{'removes':>9s}{'per 1000 edges':>16s}"
        f"{'final D':>9s}{'deg-LB (win)':>14s}\n")
summary = {}
for T in TARGETS:
    hist, added, removed, exhausted = runs[T]
    h = hist[-1]
    thr = reach_D(T, x0)
    reached = h["D"] <= thr
    summary[T] = (reached, h)
    lbw = sum(lb[T]["win"])
    rep += (f"{LABELS[T]:14s}{'yes' if reached else 'no':>9s}{h['edits']:>8,d}{h['adds']:>8,d}{h['removes']:>9,d}"
            f"{1000 * h['edits'] / m0:>16.1f}{h['D']:>9.3f}{lbw:>14,d}"
            + (f"   (threshold = model noise floor {thr:.2f})" if thr > 0.10 else "")
            + ("   (move set exhausted)" if exhausted and not reached else "") + "\n")
if "ringspoke" in TARGETS:
    rep += (f"\nRing-spoke: {ring_spoke_edits[0]:,} of the {runs['ringspoke'][0][-1]['adds']:,} added links are the rings "
            f"(every 0.5 km) and 24 spokes themselves;\nthe rest are local loop closures needed to approach "
            f"the target's zero bridges and <k> = {targets['ringspoke']['k_mean']:.1f}.\n")
rep += ("\nNotes: an upper bound is only certified when 'reached' = yes; otherwise the row shows how far the\n"
        "move set got. The trivial ceiling is m + m_target (delete everything, build the target).\n"
        "Lower bounds see only degrees: two networks with the same degree mix but different wiring\n"
        "get a lower bound of 0, so the true distance can sit far above it.\n")
print(rep)
with open(f"results/{tag}_edit_distance_report.txt", "w", encoding="utf-8") as f:
    f.write(rep)

np.savez(f"results/{tag}_edit_distance.npz",
         targets=np.array(TARGETS), metrics=np.array(ALL),
         target_vectors=np.array([[targets[T][k] for k in ALL] for T in TARGETS]),
         **{f"hist_{T}": np.array([[h[k] for k in ["edits", "adds", "removes", "D"] + ALL] for h in runs[T][0]])
            for T in TARGETS},
         lb_full=np.array([lb[T]["full"] for T in TARGETS]), lb_win=np.array([lb[T]["win"] for T in TARGETS]),
         n0=n0, m0=m0, R_win=R_win, cbd=cbd)

# ---------------------------------------------------------------- figure: metrics vs edits
log("plotting ...")
panels = RF + ["detour", "r", "D"]
fig, axs = plt.subplots(3, 3, figsize=(17, 13), dpi=170)
fig.patch.set_facecolor(SURFACE)
for ax, k in zip(axs.ravel(), panels):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=9)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for T in TARGETS:
        hist = runs[T][0]
        ex = [h["edits"] for h in hist]
        ax.plot(ex, [h[k] for h in hist], "-", color=COLORS[T], lw=2, marker=MARKERS[T], ms=6, label=LABELS[T])
        if k != "D":
            ax.axhline(targets[T][k], color=COLORS[T], ls="--", lw=1.2, alpha=0.9)
        if summary[T][0]:
            ax.scatter([ex[-1]], [hist[-1][k]], s=140, facecolors="none", edgecolors=INKD, lw=1.5, zorder=5)
    if k == "D":
        ax.axhline(0.10, color=SUBD, ls=":", lw=1.2)
        ax.text(ax.get_xlim()[0], 0.12, " reached below 0.10 (ER / SF: below their noise floor)",
                color=SUBD, fontsize=8)
    if k == "n_lam2":
        ax.set_yscale("log")
    ax.set_xscale("symlog", linthresh=20)
    ax.set_xlim(left=0)
    ax.set_title(MNAME.get(k, "normalised gap D to target (6 metrics)"), color=INKD, fontsize=11)
    ax.set_xlabel("edits (edges added + removed)", color=SUBD, fontsize=9)
axs[0, 0].legend(frameon=False, labelcolor=INKD, fontsize=9, loc="best",
                 title="solid: city being edited\ndashed: target value", title_fontsize=8)
axs[0, 0].get_legend().get_title().set_color(SUBD)
fig.suptitle(f"{CITY} → target networks: metrics vs number of edge edits  "
             f"(central window, n = {n0:,}, m = {m0:,}; circle = target reached)", color=INKD, fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(f"figures/{tag}_edit_distance_dark.png", facecolor=SURFACE)
plt.close(fig)

# ---------------------------------------------------------------- figure: where the edits are
map_targets = [T for T in ("barcelona", "ringspoke") if T in TARGETS]
if map_targets:
    fig, axs = plt.subplots(1, len(map_targets), figsize=(9 * len(map_targets), 9), dpi=170, squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    base = np.array([[pos_city[u], pos_city[v]] for u, v in W0.edges])
    for ax, T in zip(axs[0], map_targets):
        hist, added, removed, _ = runs[T]
        ax.set_facecolor(SURFACE)
        ax.add_collection(LineCollection(base, colors="#3a3a37", linewidths=0.4, zorder=1))
        if added:
            ax.add_collection(LineCollection([[pos_city[u], pos_city[v]] for u, v in added],
                                             colors="#ff8a3d", linewidths=1.4, zorder=3))
        if removed:
            ax.add_collection(LineCollection([[pos_city[u], pos_city[v]] for u, v in removed],
                                             colors="#3987e5", linewidths=1.8, zorder=4))
        ax.set_aspect("equal")
        ax.autoscale()
        ax.set_xticks([])
        ax.set_yticks([])
        for s_ in ax.spines.values():
            s_.set_color(SPINE)
        ax.legend(handles=[plt.Line2D([], [], color="#ff8a3d", lw=2, label=f"added ({len(added):,})"),
                           plt.Line2D([], [], color="#3987e5", lw=2, label=f"removed ({len(removed):,})")],
                  loc="lower left", frameon=False, labelcolor=INKD, fontsize=10)
        ax.set_title(f"Edits toward {LABELS[T]} (D = {hist[-1]['D']:.2f})", color=INKD, fontsize=12)
    fig.suptitle(f"{CITY} central window — where the greedy edits go", color=INKD, fontsize=14)
    fig.savefig(f"figures/{tag}_edit_map_dark.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
log(f"done: results/{tag}_edit_distance_report.txt + 2 figures")
