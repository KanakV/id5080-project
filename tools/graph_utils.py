"""
Shared helpers for edit_distance.py and hub_bottlenecks.py.

Copied (not imported) from the base scripts so those stay untouched:
  lambda2      <- model_comparison.py
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


def lambda2(A, planar_like=True):
    """Second-smallest eigenvalue of the normalized Laplacian of adjacency A (csr)."""
    L = csgraph.laplacian(A.astype(np.float64), normed=True)
    try:
        if planar_like:   # good separators -> cheap LU factorization, shift-invert
            vals = spla.eigsh(L.tocsc(), k=2, sigma=-1e-2, which="LM", return_eigenvectors=False)
            return float(np.sort(vals)[1])
        # expanders: Lanczos on 2I - L for the top of the spectrum (no factorization)
        M = 2.0 * sp.identity(L.shape[0], format="csr") - L.tocsr()
        vals = spla.eigsh(M, k=2, which="LA", return_eigenvectors=False, maxiter=8000, tol=1e-7)
        return float(2.0 - np.sort(vals)[0])
    except Exception as e:
        log(f"  lambda2 failed ({type(e).__name__}); retrying with shift-invert")
        try:
            vals = spla.eigsh(L.tocsc(), k=2, sigma=-1e-2, which="LM", return_eigenvectors=False)
            return float(np.sort(vals)[1])
        except Exception:
            return float("nan")


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


def band_min_cut(G, coord, q=0.10):
    """Unit-capacity max-flow min-cut between the top and bottom q-bands of `coord`
    (array aligned with list(G.nodes)). Returns (cut value, cut edges)."""
    nodes = list(G.nodes)
    hi, lo = np.quantile(coord, 1 - q), np.quantile(coord, q)
    return set_min_cut(G, [nodes[i] for i in np.flatnonzero(coord >= hi)],
                       [nodes[i] for i in np.flatnonzero(coord <= lo)])


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
