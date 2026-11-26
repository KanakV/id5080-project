"""
Shared helpers for edit_spatial.py, hub_bottlenecks.py and hub_compare.py.

Copied (not imported) from the base scripts so those stay untouched:
  sweep_cut    <- expander_bottleneck.py (Fiedler sweep)
  giant_curve  <- percolation.py (reverse Newman-Ziff union-find)
  adaptive_degree_order <- percolation.py (lazy-heap peeling)
"""
import heapq
import time

import networkx as nx
import numpy as np
import osmnx as ox
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from networkx.algorithms.flow import preflow_push
from scipy.sparse import csgraph

SEED = 42
KM_PER_DEG = 111.2


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_simple(graphml):
    """Giant component of the graphml as a simple undirected nx.Graph.
    Nodes keep x/y (lon/lat); each edge keeps the shortest parallel 'length' (m)."""
    log(f"loading {graphml}")
    G = ox.load_graphml(graphml)
    S = nx.Graph()
    S.add_nodes_from((u, {"x": d["x"], "y": d["y"]}) for u, d in G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        if u == v:
            continue
        L = float(d.get("length", 0.0))
        if not S.has_edge(u, v) or L < S[u][v]["length"]:
            S.add_edge(u, v, length=L)
    del G
    S = S.subgraph(max(nx.connected_components(S), key=len)).copy()
    log(f"giant component (simple, undirected): n={S.number_of_nodes():,}  m={S.number_of_edges():,}")
    return S


def to_km(lon, lat, lon0, lat0):
    """Equirectangular projection to km around (lon0, lat0) — fine at city scale."""
    return np.column_stack([(np.asarray(lon) - lon0) * KM_PER_DEG * np.cos(np.radians(lat0)),
                            (np.asarray(lat) - lat0) * KM_PER_DEG])


def sweep_cut(A, min_side=0.10):
    """Fiedler vector + sweep cut of adjacency A (csr, connected).
    Returns dict: lam2, phi (best), phi_bal (best with >= min_side of volume each side),
    side (bool mask of the balanced cut's side A), n_cut (edges crossing it)."""
    n = A.shape[0]
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = csgraph.laplacian(A.astype(np.float64), normed=True).tocsc()
    if n < 400:   # small corridors: dense solve is exact and cheaper than Lanczos
        vals, vecs = np.linalg.eigh(L.toarray())
    else:
        vals, vecs = spla.eigsh(L, k=2, sigma=-1e-2, which="LM")
    lam2 = float(vals[1])
    fiedler = vecs[:, 1] / np.sqrt(np.maximum(deg, 1))
    order = np.argsort(fiedler)
    rank = np.empty(n, dtype=np.int64)
    rank[order] = np.arange(n)
    indptr, indices = A.indptr, A.indices
    vol_tot = float(deg.sum())
    cut = vol = 0.0
    phi = np.full(n - 1, np.inf)
    for step, i in enumerate(order[:-1]):
        nbrs = indices[indptr[i]:indptr[i + 1]]
        inside = np.count_nonzero(rank[nbrs] <= step)      # neighbours already swept in
        cut += len(nbrs) - 2 * inside
        vol += deg[i]
        phi[step] = cut / min(vol, vol_tot - vol)
    volcum = np.cumsum(deg[order])
    balanced = (volcum[:-1] >= min_side * vol_tot) & (volcum[:-1] <= (1 - min_side) * vol_tot)
    best_bal = int(np.flatnonzero(balanced)[np.argmin(phi[balanced])]) if balanced.any() else int(np.argmin(phi))
    side = np.zeros(n, dtype=bool)
    side[order[:best_bal + 1]] = True
    coo = sp.triu(A, 1).tocoo()
    n_cut = int(np.count_nonzero(side[coo.row] != side[coo.col]))
    return dict(lam2=lam2, phi=float(phi.min()), phi_bal=float(phi[best_bal]), side=side, n_cut=n_cut)


def giant_curve(indptr, indices, removal_order):
    """removal_order: node indices, first-removed first.
    Returns S_rel[k] = giant component fraction after removing the first k nodes."""
    n = len(indptr) - 1
    parent = np.arange(n)
    size = np.ones(n, dtype=np.int64)
    present = np.zeros(n, dtype=bool)

    def find(a):
        root = a
        while parent[root] != root:
            root = parent[root]
        while parent[a] != root:
            parent[a], a = root, parent[a]
        return root

    best = 0
    S_rel = np.empty(n + 1)
    S_rel[n] = 0.0
    for pos in range(n - 1, -1, -1):          # re-add in reverse removal order
        v = removal_order[pos]
        present[v] = True
        for w in indices[indptr[v]:indptr[v + 1]]:
            if present[w]:
                ra, rb = find(v), find(w)
                if ra != rb:
                    if size[ra] < size[rb]:
                        ra, rb = rb, ra
                    parent[rb] = ra
                    size[ra] += size[rb]
        best = max(best, size[find(v)])
        S_rel[pos] = best / n
    return S_rel


def adaptive_degree_order(indptr, indices):
    """Removal order that always takes the currently highest-degree node."""
    n = len(indptr) - 1
    livedeg = np.diff(indptr).copy()
    alive = np.ones(n, dtype=bool)
    heap = [(-int(d), i) for i, d in enumerate(livedeg)]
    heapq.heapify(heap)
    order = np.empty(n, dtype=np.int64)
    k = 0
    while heap:
        d, v = heapq.heappop(heap)
        if not alive[v] or -d != livedeg[v]:
            continue
        alive[v] = False
        order[k] = v
        k += 1
        for w in indices[indptr[v]:indptr[v + 1]]:
            if alive[w]:
                livedeg[w] -= 1
                heapq.heappush(heap, (-int(livedeg[w]), w))
    return order


def set_min_cut(G, src_set, snk_set):
    """Unit-capacity max-flow min-cut between two node sets of G (simple undirected)."""
    F = nx.Graph()
    F.add_edges_from(G.edges, capacity=1)
    F.add_edges_from((("__SRC__", u) for u in src_set), capacity=np.inf)
    F.add_edges_from((("__SNK__", u) for u in snk_set), capacity=np.inf)
    val, (P1, _) = nx.minimum_cut(F, "__SRC__", "__SNK__", flow_func=preflow_push)
    P1 = set(P1) - {"__SRC__"}
    return int(val), list(nx.edge_boundary(G, P1))


def classify_cuts(pairs):
    """pairs: list of (hub_a, hub_b, min_cut). Says where each pair's min cut sits.
    A hub is 'choked' at c when its cut to every other hub (except at most one) equals c:
    then the choke point is in that hub's own catchment, whoever is on the other end.
    Returns {(hub_a, hub_b): 'near <hub>' | 'between hubs'}."""
    from collections import Counter
    hubs = {h for a, b, _ in pairs for h in (a, b)}
    choke = {}
    for h in hubs:
        cuts = [c for a, b, c in pairs if h in (a, b)]
        c, k = Counter(cuts).most_common(1)[0]
        if len(cuts) >= 3 and k >= len(cuts) - 1:
            choke[h] = c
    out = {}
    for a, b, c in pairs:
        ca, cb = choke.get(a, np.inf), choke.get(b, np.inf)
        if c == ca and ca <= cb:
            out[(a, b)] = f"near {a}"
        elif c == cb:
            out[(a, b)] = f"near {b}"
        else:
            out[(a, b)] = "between hubs"
    return out


def robustness(G, n_rand=3):
    """Schneider R (mean giant-component fraction over the whole removal curve; max 0.5) under
    random failure (mean of n_rand fixed-seed orders) and adaptive degree attack.
    Returns ({'rand': R, 'deg': R}, dict(A, orders, curves)) — the extras feed edit_spatial.py."""
    nodes = list(G.nodes)
    n = len(nodes)
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, format="csr")
    rng = np.random.default_rng(SEED)
    orders = {"rand": [rng.permutation(n) for _ in range(n_rand)],
              "deg": [adaptive_degree_order(A.indptr, A.indices)]}
    curves = {m: [giant_curve(A.indptr, A.indices, o) for o in os_] for m, os_ in orders.items()}
    R = {m: float(np.mean([c[:n].mean() for c in cs])) for m, cs in curves.items()}
    return R, dict(A=A, orders=orders, curves=curves)


# named hubs for hub_bottlenecks.py / hub_compare.py (approximate centres, picked by hand)
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


def snap_hubs(S, city_key, max_km=1.5):
    """Snap HUBS[city_key] to the nearest node of S (as hub_bottlenecks.py does); hubs more
    than max_km from any node are dropped. Returns [(name, node, snap_km)]."""
    from scipy.spatial import cKDTree
    nodes = list(S.nodes)
    lon = np.array([S.nodes[u]["x"] for u in nodes])
    lat = np.array([S.nodes[u]["y"] for u in nodes])
    lon0, lat0 = lon.mean(), lat.mean()
    tree = cKDTree(to_km(lon, lat, lon0, lat0))
    out = []
    for name, (la, lo) in HUBS[city_key].items():
        d, i = tree.query(to_km([lo], [la], lon0, lat0)[0])
        if d <= max_km:
            out.append((name, nodes[i], float(d)))
    return out


def hub_road_km(S, hubs):
    """Shortest road distance (km, Dijkstra on 'length') between every pair of snapped hubs."""
    from scipy.sparse.csgraph import dijkstra
    nodes = list(S.nodes)
    ix = {u: i for i, u in enumerate(nodes)}
    A = nx.to_scipy_sparse_array(S, nodelist=nodes, weight="length", format="csr")
    src = [ix[n] for _, n, _ in hubs]
    D = dijkstra(A, indices=src)[:, src] / 1000.0
    return {(hubs[a][0], hubs[b][0]): float(D[a, b]) for a in range(len(hubs)) for b in range(len(hubs)) if a != b}
