# =============================================================================
# HOW THE CHENNAI AND BENGALURU ROAD NETWORKS GREW OVER TIME (2018 -> 2022)
# Uses historical OpenStreetMap snapshots (Overpass "attic" data via osmnx).
#
# PURPOSE
#   Compare snapshots of each city and answer six questions:
#     1. Growth: how many nodes / roads / km were added?
#     2. Preferential attachment: do new roads attach to well-connected nodes?
#     3. Structure over time: degree distribution, clustering, path length
#        (against an ER graph of the same size).
#     4. Robustness over time: does the network get more or less robust as it
#        grows? (same random vs betweenness attack as road_robustness.py)
#     5. Where growth happens: map of the new roads, growth by distance from
#        the city centre, orientation of new vs old streets.
#     6. Centrality change: do the same corridors stay the most important?
#
# REASON
#   My robustness study looked at one moment in time. Real road networks grow,
#   and the Barabasi-Albert model says growth by preferential attachment
#   creates hubs and scale-free structure. Snapshots let me TEST that growth
#   rule on real data, and see if growth makes each city more or less robust.
#
# EXPECTATION (hypotheses to check, not results)
#   1. Growth is mostly new streets and infill, with mean degree changing little
#      (roads cannot have more than ~4-6 links per junction).
#   2. Preferential attachment is weak or absent: alpha in Pi(k) ~ k^alpha is
#      well below 1 (close to 0), because new roads join NEARBY junctions in
#      space, not the biggest hubs. Road networks stay non-scale-free.
#   3. Degree distribution stays narrow (peak at k = 1, 3 and 4) in all years;
#      path length stays much longer than ER's.
#   4. Robustness changes only a little between years. If it improves, new
#      roads are adding alternative routes.
#   5. New roads concentrate in the outer rings (sprawl) more than in the centre.
#   6. Betweenness rankings are stable (high Spearman correlation): the main
#      corridors keep their role even as new roads appear.
#
# BIG CAVEAT (state it in the report)
#   OSM edits mix real construction with better MAPPING. India was less
#   completely mapped in 2018 than in 2022, so part of the "growth" is roads
#   that already existed but were added to the map later. The script prints
#   a warning when growth looks very large, and reports how much new road
#   lies near the old network (infill, often mapping completion) vs far from
#   it (expansion, more likely real).
#   Overpass history only goes back to Sept 2012.
#
# Run:   pip install osmnx networkx numpy scipy matplotlib
#        python road_growth_over_time.py
#   First run downloads 3 snapshots per city in small tiles (slow: it can take
#   an hour or more per city). Each tile is cached in ./snapshot_cache/tiles, so
#   if it stops, just run it again and it resumes where it left off.
# =============================================================================

import os
import pickle
import time
import numpy as np
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

# ----------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------
CITIES = {
    "Chennai":   "Chennai, Tamil Nadu, India",
    "Bengaluru": "Bengaluru, Karnataka, India",
}
SNAPSHOTS = {2018: "2018-01-01T00:00:00Z",
             2020: "2020-01-01T00:00:00Z",
             2022: "2022-01-01T00:00:00Z"}
SEED = 42
BC_K = 300              # betweenness is estimated from this many sampled sources
RANDOM_RUNS = 5         # random-failure curve averaged over this many runs
PATH_SOURCES = 100      # sampled BFS sources for path length / diameter
KMAX = 6                # largest degree shown in kernels / P(k)
FIT_KMIN = 3            # fit alpha only on junction degrees k >= this
MIN_COUNT = 30          # need at least this many nodes of a degree to use it
NEW_EDGE_M = 20         # an edge is "truly new" if > this many metres from any old road
INFILL_M = 300          # new road within this distance of the old network = infill
RING_EDGES_KM = [0, 3, 6, 10, 15, 25, 1e9]   # distance rings from the city centre
TOP_PCT = 1.0           # "top corridors" = top this % of nodes by betweenness
GROWTH_WARN = 0.40      # warn if edges grow by more than this fraction
TILE_DEG = 0.03         # download tile size in degrees (~6.5 km); smaller = safer, slower
MAX_RETRIES = 5         # retries per tile when the Overpass server drops the connection

# display
# For a zoomable window uncomment BEFORE plotting (pip install pyqt5):
# import matplotlib
# matplotlib.use("QtAgg")
LINE_W = 0.5
DPI = 300
ZOOM = None             # e.g. (80.20, 80.30, 12.95, 13.10) = (lon_min, lon_max, lat_min, lat_max)
OUT_DIR = "growth_plots"
CACHE_DIR = "snapshot_cache"

GRID = np.linspace(0, 1, 401)
C_CITY = {"Chennai": "#FF7A1A", "Bengaluru": "#0B8F62"}
rng = np.random.default_rng(SEED)
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)
TILE_DIR = os.path.join(CACHE_DIR, "tiles")
os.makedirs(TILE_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# 1. DOWNLOAD HISTORICAL SNAPSHOTS (cached on disk)
# ----------------------------------------------------------------------------
def to_simple(Gm):
    """osmnx MultiDiGraph -> simple undirected graph with length, x, y."""
    U = nx.Graph()
    for u, v, d in Gm.edges(data=True):
        if u == v:
            continue
        L = float(d.get("length", 1.0))
        if U.has_edge(u, v):
            U[u][v]["length"] = min(U[u][v]["length"], L)
        else:
            U.add_edge(u, v, length=L)
    for n, d in Gm.nodes(data=True):
        if n in U:
            U.nodes[n]["x"], U.nodes[n]["y"] = d["x"], d["y"]
    return U


def tile_polygons(place):
    """Split the city boundary into small tiles (a big single query makes the
    Overpass server drop the connection)."""
    from shapely.geometry import box, Polygon, MultiPolygon
    from shapely.ops import unary_union
    geom = ox.geocode_to_gdf(place).geometry
    poly = geom.union_all() if hasattr(geom, "union_all") else geom.unary_union
    x0, y0, x1, y1 = poly.bounds
    tiles = []
    for i, x in enumerate(np.arange(x0, x1, TILE_DEG)):
        for j, y in enumerate(np.arange(y0, y1, TILE_DEG)):
            cell = box(x, y, min(x + TILE_DEG, x1), min(y + TILE_DEG, y1)).intersection(poly)
            parts = [g for g in getattr(cell, "geoms", [cell]) if isinstance(g, (Polygon, MultiPolygon))]
            if not parts:
                continue
            cell = unary_union(parts)
            if cell.area > 1e-7:
                tiles.append(((i, j), cell))
    return tiles


def download_tile(poly):
    """One small Overpass request with retries (the server sometimes cuts the
    connection or is busy)."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            Gm = ox.graph_from_polygon(poly, network_type="drive", simplify=False,
                                       retain_all=True, truncate_by_edge=True)
            return to_simple(Gm)
        except Exception as e:
            if "InsufficientResponse" in type(e).__name__:
                return nx.Graph()

            response = getattr(e, "response", None)
            status = getattr(response, "status_code", None)

            print(f"      attempt {attempt}/{MAX_RETRIES} failed")
            print(f"      Error type: {type(e).__name__}")
            print(f"      HTTP status: {status}")
            print(f"      Error details: {str(e)[:500]}")

            if attempt < MAX_RETRIES:
                wait = 30 * attempt
                print(f"      Retrying in {wait}s")
                time.sleep(wait)
    raise RuntimeError("Overpass kept failing for this tile. Try again later, make TILE_DEG "
                       "smaller, or check https://overpass-api.de/api/status")


def download_raw(city, place, year):
    """Un-simplified drive network as it was on SNAPSHOTS[year], downloaded tile by tile.
    Every finished tile is saved to disk, so a crash can be resumed by re-running."""
    ox.settings.use_cache = False          # make sure each date is really downloaded
    ox.settings.requests_timeout = 600

# Use an alternative Overpass server with historical (Attic) data
    ox.settings.overpass_url = "https://overpass.private.coffee/api"
    ox.settings.overpass_settings = (
        f'[out:json][timeout:{{timeout}}][date:"{SNAPSHOTS[year]}"]{{maxsize}}')
    tiles = tile_polygons(place)
    print(f"  {city} {year}: {len(tiles)} tiles")
    U = nx.Graph()
    for n_done, ((i, j), poly) in enumerate(tiles, 1):
        path = os.path.join(TILE_DIR, f"{city}_{year}_{i}_{j}.pkl")
        if os.path.exists(path):
            with open(path, "rb") as f:
                T = pickle.load(f)
        else:
            T = download_tile(poly)
            with open(path, "wb") as f:
                pickle.dump(T, f)
        U.add_nodes_from(T.nodes(data=True))
        U.add_edges_from(T.edges(data=True))
        if n_done % 5 == 0 or n_done == len(tiles):
            print(f"    tile {n_done}/{len(tiles)}")
    return U


def get_snapshot(city, place, year):
    path = os.path.join(CACHE_DIR, f"{city}_{year}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    G = download_raw(city, place, year)
    with open(path, "wb") as f:
        pickle.dump(G, f)
    return G


# ----------------------------------------------------------------------------
# 2. CONTRACT DEGREE-2 NODES (so nodes = junctions and dead ends, edges = streets)
# ----------------------------------------------------------------------------
def contract(G_raw):
    H = G_raw.copy()
    for n in [n for n in G_raw if G_raw.degree(n) == 2]:
        if n not in H or H.degree(n) != 2:
            continue
        a, b = list(H.neighbors(n))
        if a == b or H.has_edge(a, b):
            continue
        H.add_edge(a, b, length=H[a][n]["length"] + H[n][b]["length"])
        H.remove_node(n)
    H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
    return H


# ----------------------------------------------------------------------------
# 3. STRUCTURE METRICS
# ----------------------------------------------------------------------------
def sampled_paths(G, n_src):
    nodes = list(G.nodes)
    src = rng.choice(len(nodes), min(n_src, len(nodes)), replace=False)
    vals, ecc = [], 0
    for i in src:
        d = np.array([v for v in nx.single_source_shortest_path_length(G, nodes[i]).values() if v > 0])
        vals.append(d); ecc = max(ecc, d.max())
    return np.concatenate(vals).mean(), ecc


def structure_metrics(G):
    n, m = G.number_of_nodes(), G.number_of_edges()
    deg = np.array([d for _, d in G.degree()])
    Pk = np.array([(deg == k).mean() for k in range(1, KMAX + 1)])
    L, diam = sampled_paths(G, PATH_SOURCES)
    ER = nx.gnm_random_graph(n, m, seed=SEED)
    ER = ER.subgraph(max(nx.connected_components(ER), key=len)).copy()
    L_er, _ = sampled_paths(ER, PATH_SOURCES)
    return dict(N=n, M=m, km=sum(d["length"] for *_, d in G.edges(data=True)) / 1000,
                k=2 * m / n, Pk=Pk, C=nx.average_clustering(G), L=L, diam=diam,
                L_er=L_er, C_er=nx.average_clustering(ER))


# ----------------------------------------------------------------------------
# 4. ROBUSTNESS (random vs betweenness attack; same metric as road_robustness.py)
# ----------------------------------------------------------------------------
def giant_curve(adj, order):
    n = len(adj)
    parent = np.arange(n); size = np.ones(n, dtype=np.int64)
    present = np.zeros(n, dtype=bool)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    giant = np.zeros(n + 1, dtype=np.int64); best = 0
    for r in range(n, 0, -1):
        giant[r] = best
        v = order[r - 1]; present[v] = True
        for w in adj[v]:
            if present[w]:
                a, b = find(v), find(w)
                if a != b:
                    if size[a] < size[b]:
                        a, b = b, a
                    parent[b] = a; size[a] += size[b]
        best = max(best, size[find(v)])
    giant[0] = best
    return giant


def to_grid(giant):
    n = len(giant) - 1
    f = np.arange(n + 1) / n
    rel = giant / giant[0]
    below = np.where(rel < 0.01)[0]
    return np.interp(GRID, f, rel), giant.mean() / n, (f[below[0]] if len(below) else 1.0)


def robustness(G, bc):
    nodes = list(G.nodes); idx = {v: i for i, v in enumerate(nodes)}; n = len(nodes)
    adj = [[idx[w] for w in G.neighbors(v)] for v in nodes]
    out = {}
    runs = [to_grid(giant_curve(adj, rng.permutation(n))) for _ in range(RANDOM_RUNS)]
    out["random"] = (np.mean([r[0] for r in runs], axis=0), np.mean([r[1] for r in runs]),
                     np.mean([r[2] for r in runs]))
    order = np.array([idx[v] for v in sorted(nodes, key=lambda v: -bc[v])])
    out["attack"] = to_grid(giant_curve(adj, order))
    return out


# ----------------------------------------------------------------------------
# 5. PREFERENTIAL ATTACHMENT KERNEL
# ----------------------------------------------------------------------------
def ekey(u, v):
    return (u, v) if u < v else (v, u)


def fit_alpha(ks, pi, nk):
    sel = (ks >= FIT_KMIN) & (nk >= MIN_COUNT) & (pi > 0)
    if sel.sum() < 2:
        return np.nan
    return np.polyfit(np.log(ks[sel]), np.log(pi[sel]), 1)[0]


def attach_kernel(G0, G1):
    """Pi(k): average number of new links gained by an old node that had k links."""
    d0, d1 = dict(G0.degree()), dict(G1.degree())
    E0 = {ekey(u, v) for u, v in G0.edges}
    new_link = {}
    for u, v in G1.edges:
        if ekey(u, v) in E0:
            continue
        if u in d0 and v not in d0:
            new_link[u] = new_link.get(u, 0) + 1
        elif v in d0 and u not in d0:
            new_link[v] = new_link.get(v, 0) + 1
    old = [n for n in d0 if n in d1]
    k0 = np.array([d0[n] for n in old])
    gain = np.array([max(d1[n] - d0[n], 0) for n in old])
    newl = np.array([new_link.get(n, 0) for n in old])
    ks = np.arange(1, KMAX + 1)
    nk = np.array([(k0 == k).sum() for k in ks])
    safe = np.maximum(nk, 1)
    pi_gain = np.array([gain[k0 == k].sum() for k in ks]) / safe
    pi_new = np.array([newl[k0 == k].sum() for k in ks]) / safe
    return dict(ks=ks, nk=nk, pi_gain=pi_gain, pi_new=pi_new,
                a_gain=fit_alpha(ks, pi_gain, nk), a_new=fit_alpha(ks, pi_new, nk))


# ----------------------------------------------------------------------------
# 6. NEW ROADS: WHERE, HOW FAR FROM THE OLD NETWORK, WHICH ORIENTATION
# ----------------------------------------------------------------------------
def edge_arrays(G, edges):
    """Coordinates (x1, y1, x2, y2) and length of each edge."""
    if len(edges) == 0:
        z = np.zeros(0)
        return z, z, z, z, z
    e = np.array(edges)
    x1 = np.array([G.nodes[a]["x"] for a in e[:, 0]]); y1 = np.array([G.nodes[a]["y"] for a in e[:, 0]])
    x2 = np.array([G.nodes[b]["x"] for b in e[:, 1]]); y2 = np.array([G.nodes[b]["y"] for b in e[:, 1]])
    L = np.array([G[a][b]["length"] for a, b in edges])
    return x1, y1, x2, y2, L


def orientation_entropy(x1, y1, x2, y2, L, nbins=36):
    b = np.degrees(np.arctan2((x2 - x1) * np.cos(np.radians(y1)), y2 - y1)) % 180
    both = np.concatenate([b, (b + 180) % 360]); w = np.concatenate([L, L])
    counts, _ = np.histogram((both + 180 / nbins) % 360, bins=nbins, range=(0, 360), weights=w)
    p = counts / counts.sum()
    return -np.sum(p[p > 0] * np.log(p[p > 0]))


def new_edge_analysis(G0, G1):
    lat0 = np.mean([d["y"] for _, d in G0.nodes(data=True)])
    cx = np.median([d["x"] for _, d in G0.nodes(data=True)])
    cy = np.median([d["y"] for _, d in G0.nodes(data=True)])
    mx = lambda x: x * 111320 * np.cos(np.radians(lat0))
    my = lambda y: y * 110540

    e0 = list(G0.edges)
    x1, y1, x2, y2, L0 = edge_arrays(G0, e0)
    cnt = np.maximum(1, np.ceil(L0 / 10).astype(int))                  # points every ~10 m
    rep = np.repeat(np.arange(len(L0)), cnt)
    t = (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt) + 0.5) / np.repeat(cnt, cnt)
    px = x1[rep] + (x2 - x1)[rep] * t; py = y1[rep] + (y2 - y1)[rep] * t
    tree = cKDTree(np.column_stack((mx(px), my(py))))

    E0 = {ekey(u, v) for u, v in e0}
    cand = [e for e in G1.edges if ekey(*e) not in E0]
    a1, b1, a2, b2, Ln = edge_arrays(G1, cand)
    mid_x, mid_y = (a1 + a2) / 2, (b1 + b2) / 2
    dist, _ = tree.query(np.column_stack((mx(mid_x), my(mid_y))))
    truly = dist > NEW_EDGE_M
    infill = truly & (dist <= INFILL_M)

    def ring_km(xx, yy):
        d = np.hypot((xx - cx) * 111320 * np.cos(np.radians(lat0)), (yy - cy) * 110540) / 1000
        return np.digitize(d, RING_EDGES_KM) - 1
    r_old = ring_km((x1 + x2) / 2, (y1 + y2) / 2)
    r_new = ring_km(mid_x[truly], mid_y[truly])
    nr = len(RING_EDGES_KM) - 1
    old_len = np.array([L0[r_old == i].sum() for i in range(nr)]) / 1000
    new_len = np.array([Ln[truly][r_new == i].sum() for i in range(nr)]) / 1000

    keep = np.flatnonzero(truly)
    return dict(new_km=Ln[truly].sum() / 1000, old_km=L0.sum() / 1000,
                infill_share=Ln[infill].sum() / max(Ln[truly].sum(), 1e-9),
                ring_growth=100 * new_len / np.maximum(old_len, 1e-9), old_len=old_len,
                H_old=orientation_entropy(x1, y1, x2, y2, L0),
                H_new=orientation_entropy(a1[keep], b1[keep], a2[keep], b2[keep], Ln[keep])
                if len(keep) > 10 else np.nan,
                mean_len_old=L0.mean(), mean_len_new=Ln[truly].mean() if truly.any() else np.nan,
                segs_new=np.stack((np.column_stack((a1[keep], b1[keep])),
                                   np.column_stack((a2[keep], b2[keep]))), axis=1),
                segs_all=np.stack((np.column_stack((x1, y1)), np.column_stack((x2, y2))), axis=1))


# ----------------------------------------------------------------------------
# 7. PLOT HELPERS
# ----------------------------------------------------------------------------
def prep_ax(ax, X, Y, title):
    ax.set_aspect(1 / np.cos(np.radians(np.mean(Y))))
    ax.set_xlim(np.min(X), np.max(X)); ax.set_ylim(np.min(Y), np.max(Y))
    if ZOOM:
        ax.set_xlim(ZOOM[0], ZOOM[1]); ax.set_ylim(ZOOM[2], ZOOM[3])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10)


def clean(ax):
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def graph_segs(G):
    return [[(G.nodes[u]["x"], G.nodes[u]["y"]), (G.nodes[v]["x"], G.nodes[v]["y"])] for u, v in G.edges]


# ----------------------------------------------------------------------------
# 8. ONE CITY
# ----------------------------------------------------------------------------
def analyse(city, place):
    t_start = time.time()
    print(f"\n=== {city} ===")
    years = sorted(SNAPSHOTS)
    raw = {y: get_snapshot(city, place, y) for y in years}
    sizes = {y: (raw[y].number_of_nodes(), raw[y].number_of_edges()) for y in years}
    for y in years:
        print(f"  {y}: {sizes[y][0]:,} nodes, {sizes[y][1]:,} edges (raw, un-simplified)")
    if len({sizes[y] for y in years}) == 1:
        print("  WARNING: all snapshots are identical. Historical queries probably did not work "
              "(check osmnx version / delete ./snapshot_cache and retry).")
    g = sizes[years[-1]][1] / max(sizes[years[0]][1], 1) - 1
    if g > GROWTH_WARN:
        print(f"  NOTE: edges grew {100 * g:.0f}% from {years[0]} to {years[-1]}. "
              f"A large part is probably better mapping, not new construction.")

    con = {y: contract(raw[y]) for y in years}
    st = {y: structure_metrics(con[y]) for y in years}

    bc, rob = {}, {}
    for y in years:
        print(f"  betweenness + robustness {y} ...")
        bc[y] = nx.betweenness_centrality(con[y], k=min(BC_K, con[y].number_of_nodes()), seed=SEED)
        rob[y] = robustness(con[y], bc[y])

    pairs = list(zip(years[:-1], years[1:])) + ([(years[0], years[-1])] if len(years) > 2 else [])
    pa = {p: attach_kernel(raw[p[0]], raw[p[1]]) for p in pairs}
    ne = new_edge_analysis(raw[years[0]], raw[years[-1]])

    y0, y1 = years[0], years[-1]
    common = sorted(set(con[y0]) & set(con[y1]))
    rho = spearmanr([bc[y0][v] for v in common], [bc[y1][v] for v in common]).correlation
    topn = lambda y: set(sorted(con[y], key=lambda v: -bc[y][v])[:max(1, int(TOP_PCT / 100 * con[y].number_of_nodes()))])
    t0s, t1s = topn(y0), topn(y1)
    jac = len(t0s & t1s) / len(t0s | t1s)

    cols = plt.cm.plasma(np.linspace(0.15, 0.8, len(years)))

    # ---------- Fig 1: growth + structure ----------
    fig, ax = plt.subplots(2, 3, figsize=(16, 8))
    ax[0, 0].plot(years, [st[y]["N"] for y in years], "o-", color=C_CITY[city], lw=2.5)
    ax[0, 0].set_title("Nodes (junctions + dead ends)", fontweight="bold")
    ax[0, 1].plot(years, [st[y]["km"] for y in years], "o-", color=C_CITY[city], lw=2.5)
    ax[0, 1].set_title("Total road length (km)", fontweight="bold")
    w = 0.8 / len(years)
    for i, y in enumerate(years):
        ax[0, 2].bar(np.arange(1, KMAX + 1) + (i - (len(years) - 1) / 2) * w, st[y]["Pk"], w,
                     color=cols[i], label=str(y))
    ax[0, 2].set_title("Degree distribution P(k)", fontweight="bold")
    ax[0, 2].set_xlabel("degree k"); ax[0, 2].legend(frameon=False)
    ax[1, 0].plot(years, [st[y]["k"] for y in years], "o-", color=C_CITY[city], lw=2.5)
    ax[1, 0].set_title("Mean degree <k>", fontweight="bold")
    ax[1, 1].plot(years, [st[y]["L"] for y in years], "o-", color=C_CITY[city], lw=2.5, label="road")
    ax[1, 1].plot(years, [st[y]["L_er"] for y in years], "o--", color="#2F6BFF", lw=2, label="ER (same N, M)")
    ax[1, 1].set_title("Average path length (hops)", fontweight="bold"); ax[1, 1].legend(frameon=False)
    ax[1, 2].plot(years, [st[y]["C"] for y in years], "o-", color=C_CITY[city], lw=2.5, label="road")
    ax[1, 2].plot(years, [st[y]["C_er"] for y in years], "o--", color="#2F6BFF", lw=2, label="ER (same N, M)")
    ax[1, 2].set_title("Clustering coefficient", fontweight="bold"); ax[1, 2].legend(frameon=False)
    for a in ax.flat:
        clean(a)
        if a is not ax[0, 2]:
            a.set_xticks(years)
    fig.suptitle(f"{city}: growth and structure over time", fontsize=14, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f"{city}_1_structure_over_time.png"), dpi=200)

    # ---------- Fig 2: preferential attachment ----------
    p = (y0, y1) if (y0, y1) in pa else pairs[-1]
    K = pa[p]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for a, key, akey, ttl in ((ax[0], "pi_gain", "a_gain", "Links gained by old nodes"),
                              (ax[1], "pi_new", "a_new", "Links from NEW nodes to old nodes")):
        ok = (K["pi_gain"] > 0) if key == "pi_gain" else (K["pi_new"] > 0)
        a.loglog(K["ks"][ok], K[key][ok], "o-", color=C_CITY[city], lw=2.5, ms=8,
                 label=f"data (alpha = {K[akey]:.2f})")
        if ok.any():
            k0 = K["ks"][ok][-1]; v0 = K[key][ok][-1]
            kk = np.array([1, KMAX])
            a.loglog(kk, v0 * (kk / k0) ** 1, "k--", lw=1, label="linear PA (alpha = 1)")
            a.loglog(kk, np.full(2, v0), ":", color="gray", lw=1.5, label="no PA (alpha = 0)")
        a.set_xlabel("old degree k"); a.set_ylabel(r"$\Pi(k)$: new links per node")
        a.set_title(f"{ttl}, {p[0]} -> {p[1]}", fontweight="bold"); a.legend(frameon=False); clean(a)
    fig.suptitle(f"{city}: attachment kernel (fit on k >= {FIT_KMIN})", fontsize=13, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f"{city}_2_preferential_attachment.png"), dpi=200)

    # ---------- Fig 3: robustness over time ----------
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    for a, mode, ttl in ((ax[0], "random", "Random failure"), (ax[1], "attack", "Betweenness attack")):
        for i, y in enumerate(years):
            a.plot(GRID, rob[y][mode][0], color=cols[i], lw=2.5, label=str(y))
        a.set_xlabel("fraction of nodes removed  f"); a.set_ylabel(r"$P_\infty(f)/P_\infty(0)$")
        a.set_title(ttl, fontweight="bold"); a.legend(frameon=False); clean(a)
    for mode, ls in (("random", "o-"), ("attack", "s--")):
        ax[2].plot(years, [rob[y][mode][1] for y in years], ls, lw=2.3, color=C_CITY[city], label=mode)
    ax[2].set_ylabel("robustness R (area under curve)"); ax[2].set_xticks(years)
    ax[2].set_title("R over time (higher = more robust)", fontweight="bold")
    ax[2].legend(frameon=False); clean(ax[2])
    fig.suptitle(f"{city}: robustness over time", fontsize=14, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f"{city}_3_robustness_over_time.png"), dpi=200)

    # ---------- Fig 4: where the new roads are ----------
    X = ne["segs_all"][:, :, 0].ravel(); Y = ne["segs_all"][:, :, 1].ravel()
    fig, ax = plt.subplots(1, 2, figsize=(17, 8), gridspec_kw={"width_ratios": [1.6, 1]})
    ax[0].add_collection(LineCollection(ne["segs_all"], colors="#cccccc", linewidths=LINE_W))
    ax[0].add_collection(LineCollection(ne["segs_new"], colors="#d62728", linewidths=LINE_W * 2.2))
    prep_ax(ax[0], X, Y, f"New roads {y0} -> {y1} (red): {ne['new_km']:,.0f} km, "
                         f"{100 * ne['infill_share']:.0f}% within {INFILL_M} m of the old network")
    labels = [f"{a:.0f}-{b:.0f}" if b < 1e8 else f">{a:.0f}" for a, b in zip(RING_EDGES_KM[:-1], RING_EDGES_KM[1:])]
    ax[1].bar(labels, ne["ring_growth"], color=C_CITY[city])
    ax[1].set_xlabel("distance from city centre (km)"); ax[1].set_ylabel(f"new road length / {y0} road length (%)")
    ax[1].set_title("Growth by distance from centre", fontweight="bold"); clean(ax[1])
    fig.suptitle(f"{city}: where the network grew", fontsize=14, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f"{city}_4_new_roads.png"), dpi=DPI, bbox_inches="tight")

    # ---------- Fig 5: main corridors in first and last year ----------
    fig, ax = plt.subplots(1, 2, figsize=(17, 8))
    for a, y, top in ((ax[0], y0, t0s), (ax[1], y1, t1s)):
        G = con[y]
        a.add_collection(LineCollection(graph_segs(G), colors="#cccccc", linewidths=LINE_W))
        hot = [[(G.nodes[u]["x"], G.nodes[u]["y"]), (G.nodes[v]["x"], G.nodes[v]["y"])]
               for u, v in G.edges if u in top or v in top]
        a.add_collection(LineCollection(hot, colors="#d62728", linewidths=LINE_W * 3))
        xs = [d["x"] for _, d in G.nodes(data=True)]; ys = [d["y"] for _, d in G.nodes(data=True)]
        prep_ax(a, xs, ys, f"{y}: roads at the top {TOP_PCT:g}% most central junctions (red)")
    fig.suptitle(f"{city}: main corridors. Spearman rank correlation of betweenness = {rho:.2f}, "
                 f"overlap of top sets (Jaccard) = {jac:.2f}", fontsize=13, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f"{city}_5_corridors.png"), dpi=DPI, bbox_inches="tight")

    print(f"  done in {time.time() - t_start:.0f}s")
    return dict(years=years, st=st, rob=rob, pa=pa, pairs=pairs, ne=ne, rho=rho, jac=jac,
                growth=100 * g)


# ----------------------------------------------------------------------------
# 9. MAIN
# ----------------------------------------------------------------------------
def main():
    R = {c: analyse(c, p) for c, p in CITIES.items()}
    years = R[next(iter(R))]["years"]

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for c in R:
        for mode, ls in (("random", "o-"), ("attack", "s--")):
            ax[0].plot(years, [R[c]["rob"][y][mode][1] for y in years], ls, lw=2.3,
                       color=C_CITY[c], label=f"{c}: {mode}")
    ax[0].set_xticks(years); ax[0].set_ylabel("robustness R"); ax[0].legend(frameon=False)
    ax[0].set_title("Robustness over time", fontweight="bold"); clean(ax[0])
    names = list(R)
    ax[1].bar(names, [R[c]["pa"][R[c]["pairs"][-1]]["a_gain"] for c in names], color=[C_CITY[c] for c in names])
    ax[1].axhline(1, color="k", ls="--", lw=1); ax[1].text(0.5, 1.02, "linear PA", ha="center", fontsize=9)
    ax[1].set_ylabel("attachment exponent alpha"); ax[1].set_title("Preferential attachment exponent", fontweight="bold")
    clean(ax[1])
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "cities_over_time.png"), dpi=200)

    print("\n" + "=" * 78)
    print(f"{'':<42}" + "".join(f"{c:>18}" for c in R))
    print("-" * 78)
    for y in years:
        for label, fn in (("N (nodes)", lambda r: f"{r['st'][y]['N']:,}"),
                          ("road km", lambda r: f"{r['st'][y]['km']:,.0f}"),
                          ("mean degree", lambda r: f"{r['st'][y]['k']:.2f}"),
                          ("avg path (hops) / ER", lambda r: f"{r['st'][y]['L']:.0f} / {r['st'][y]['L_er']:.0f}"),
                          ("R random", lambda r: f"{r['rob'][y]['random'][1]:.3f}"),
                          ("R betweenness attack", lambda r: f"{r['rob'][y]['attack'][1]:.3f}"),
                          ("f_c betweenness attack", lambda r: f"{r['rob'][y]['attack'][2]:.3f}")):
            print(f"{label + ' ' + str(y):<42}" + "".join(f"{fn(R[c]):>18}" for c in R))
        print("-" * 78)
    for p in R[names[0]]["pairs"]:
        print(f"{'PA exponent alpha ' + str(p[0]) + '->' + str(p[1]):<42}" +
              "".join(f"{R[c]['pa'][p]['a_gain']:>18.2f}" for c in R))
    print(f"{'new road km (truly new)':<42}" + "".join(f"{R[c]['ne']['new_km']:>18,.0f}" for c in R))
    print(f"{'  share near old network (infill)':<42}" + "".join(f"{100 * R[c]['ne']['infill_share']:>17.0f}%" for c in R))
    print(f"{'orientation entropy old / new':<42}" + "".join(
        f"{R[c]['ne']['H_old']:>11.2f} /{R[c]['ne']['H_new']:>5.2f}" for c in R))
    print(f"{'betweenness Spearman first vs last':<42}" + "".join(f"{R[c]['rho']:>18.2f}" for c in R))
    print(f"{'top-corridor overlap (Jaccard)':<42}" + "".join(f"{R[c]['jac']:>18.2f}" for c in R))
    print("=" * 78)
    print("alpha near 1 = preferential attachment; near 0 = new links attach regardless of degree.")
    print("Remember: part of the growth may be mapping effort, not construction.")
    print(f"\nPlots saved in ./{OUT_DIR}/")
    plt.show()


if __name__ == "__main__":
    main()