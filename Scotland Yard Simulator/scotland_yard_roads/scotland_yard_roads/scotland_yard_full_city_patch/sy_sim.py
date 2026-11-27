#!/usr/bin/env python3
"""
Road-only Scotland Yard pursuit-evasion simulator for network science experiments.

One Mr. X vs several detectives on ANY undirected graph. Every edge is a road and
every player moves along one edge per turn. No tickets, no transport types.

Examples
--------
  python sy_sim.py --graph taxi --games 2000
  python sy_sim.py --graph er --n 199 --games 2000 --instances 10
  python sy_sim.py --graph ws --games 2000 --detectives-strategy strategic
  python sy_sim.py --graph edgelist:my_chennai.csv --games 1000 --max-rounds 40
  python sy_sim.py --graph osm:"Adyar, Chennai, India" --games 500
  python sy_sim.py --graph taxi --replay replay.json      # save one game trajectory

Outputs (in --out, default ./results): games.csv, graph_metrics.csv, summary.json
"""
import argparse
import csv
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, shortest_path

HERE = os.path.dirname(os.path.abspath(__file__))
ORIGINAL_STARTS = [13, 26, 29, 34, 50, 53, 91, 94, 103, 112, 117, 132, 138, 141,
                   155, 174, 197, 198]            # labels from the original repo
ORIGINAL_REVEALS = (3, 8, 13, 18, 24)             # reveal rounds of the board game
ORIGINAL_ROUNDS = 24


# ----------------------------------------------------------------------------
# Graph construction
# ----------------------------------------------------------------------------
def largest_component(G):
    comp = max(nx.connected_components(G), key=len)
    return G.subgraph(comp).copy()


def relabel_zero_based(G):
    """Return (graph with nodes 0..n-1, list mapping index -> original label)."""
    labels = list(G.nodes())
    H = nx.relabel_nodes(G, {l: i for i, l in enumerate(labels)}, copy=True)
    return H, labels


def read_edge_file(path):
    G = nx.Graph()
    with open(path, newline="") as f:
        for row in csv.reader(f):
            row = [c.strip() for c in row if c.strip() != ""]
            if len(row) < 2:
                continue
            try:
                a, b = int(row[0]), int(row[1])
            except ValueError:          # header line
                continue
            if a != b:
                G.add_edge(a, b)
    return G


def ball_subgraph(G, size, rng, start=None):
    """Connected subgraph of ~`size` nodes grown by BFS from a random node.
    Handy for cutting a city into pieces comparable to the 199-node board."""
    if start is None:
        start = rng.choice(list(G.nodes()))
    seen, frontier = [start], [start]
    mark = {start}
    while frontier and len(seen) < size:
        nxt = []
        for u in frontier:
            for v in G.neighbors(u):
                if v not in mark:
                    mark.add(v)
                    seen.append(v)
                    nxt.append(v)
                    if len(seen) >= size:
                        break
            if len(seen) >= size:
                break
        frontier = nxt
    return G.subgraph(seen).copy()


def graph_from_osm(place):
    try:
        import osmnx as ox
    except ImportError:
        sys.exit("osmnx is not installed: pip install osmnx")
    D = ox.graph_from_place(place, network_type="drive")
    G = nx.Graph(ox.convert.to_undirected(D))
    G.remove_edges_from(nx.selfloop_edges(G))
    G = nx.Graph(G.edges())              # keep topology only
    return largest_component(G)


def build_graph(args, instance):
    """Return (networkx graph with arbitrary labels, name, fixed_starts_or_None)."""
    rng = random.Random(args.seed * 1000 + instance)
    spec = args.graph
    if spec == "taxi":
        G = read_edge_file(os.path.join(HERE, "data", "board_edges.csv"))
        return G, "taxi_board", ORIGINAL_STARTS
    if spec == "er":
        n = args.n
        p = args.mean_degree / (n - 1)
        for _ in range(200):
            G = nx.gnp_random_graph(n, p, seed=rng.randrange(2**31))
            if nx.is_connected(G):
                return G, "er", None
        G = largest_component(G)           # fall back: giant component
        return G, "er", None
    if spec == "ws":
        G = nx.connected_watts_strogatz_graph(args.n, args.k, args.p,
                                              tries=200, seed=rng.randrange(2**31))
        return G, "watts_strogatz", None
    if spec == "ba":
        G = nx.barabasi_albert_graph(args.n, args.m, seed=rng.randrange(2**31))
        return G, "barabasi_albert", None
    if spec.startswith("edgelist:"):
        G = largest_component(read_edge_file(spec.split(":", 1)[1]))
        if args.subgraph:
            G = ball_subgraph(G, args.subgraph, rng)
        return G, os.path.splitext(os.path.basename(spec.split(":", 1)[1]))[0], None
    if spec.startswith("osm:"):
        G = graph_from_osm(spec.split(":", 1)[1])
        if args.subgraph:
            G = ball_subgraph(G, args.subgraph, rng)
        return G, "osm", None
    sys.exit(f"unknown --graph {spec!r}")


def graph_metrics(G, rng, sample=300):
    """Compute exact metrics on small graphs and sampled estimates on large graphs.

    Full-city street graphs can contain 70k--160k nodes. Sampling fewer source
    nodes keeps map loading responsive while still reporting useful estimates.
    """
    n, m = G.number_of_nodes(), G.number_of_edges()
    deg = np.fromiter((d for _, d in G.degree()), dtype=np.int32, count=n)
    if n <= 1500:
        diam = nx.diameter(G)
        apl = nx.average_shortest_path_length(G)
        bc = nx.betweenness_centrality(G)
    else:
        path_sample = min(18 if n > 50000 else 40 if n > 10000 else sample, n)
        centrality_sample = min(10 if n > 50000 else 20 if n > 10000 else sample, n)
        srcs = rng.sample(list(G.nodes()), path_sample)
        ecc, tot, cnt = 0, 0, 0
        for s in srcs:
            d = nx.single_source_shortest_path_length(G, s)
            ecc = max(ecc, max(d.values()))
            tot += sum(d.values())
            cnt += len(d) - 1
        diam, apl = ecc, tot / max(1, cnt)  # sampled diameter is a lower bound
        bc = nx.betweenness_centrality(G, k=centrality_sample, seed=1)
    try:
        assort = nx.degree_assortativity_coefficient(G)
    except Exception:
        assort = float("nan")
    return {
        "n": n, "m": m, "mean_degree": round(float(deg.mean()), 4),
        "max_degree": int(deg.max()), "min_degree": int(deg.min()),
        "degree_std": round(float(deg.std()), 4),
        "density": round(nx.density(G), 8),
        "avg_clustering": round(nx.average_clustering(G), 4),
        "transitivity": round(nx.transitivity(G), 4),
        "diameter": int(diam), "avg_path_length": round(float(apl), 4),
        "degree_assortativity": round(float(assort), 4),
        "n_bridges": sum(1 for _ in nx.bridges(G)),
        "n_articulation_points": sum(1 for _ in nx.articulation_points(G)),
        "max_betweenness": round(max(bc.values()), 5),
        "pct_degree1": round(float((deg == 1).mean()), 4),
    }


# ----------------------------------------------------------------------------
# Game engine
# ----------------------------------------------------------------------------
class Board:
    """Zero-based graph plus bounded shortest-path lookups for large networks."""

    def __init__(self, G):
        self.G = G
        self.n = G.number_of_nodes()
        self.adj = [list(G.neighbors(i)) for i in range(self.n)]
        self.deg = np.fromiter((len(a) for a in self.adj), dtype=np.int32, count=self.n)
        self.A = csr_matrix(nx.to_scipy_sparse_array(G, nodelist=range(self.n), dtype=np.int8))
        self._full = None
        self._rows = {}
        # All-pairs shortest paths are convenient on the original board, but
        # would consume prohibitive memory on a whole-city road graph.
        if self.n <= 3000:
            self._full = shortest_path(self.A, unweighted=True).astype(np.int32)
        self._row_limit = 48 if self.n > 10000 else 2000

    def dist(self, v):
        if self._full is not None:
            return self._full[v]
        r = self._rows.get(v)
        if r is None:
            r = shortest_path(self.A, unweighted=True, indices=v).astype(np.int32)
            if len(self._rows) >= self._row_limit:
                # Keep the cache bounded: each city distance row can be >0.5 MB.
                self._rows.clear()
            self._rows[v] = r
        return r

    def distance_to_any(self, sources):
        """Distance from every node to its nearest source, in one graph traversal."""
        sources = np.asarray(list(sources), dtype=np.int32)
        if sources.size == 0:
            return np.full(self.n, self.n + 1, dtype=np.int32)
        if sources.size == 1:
            return self.dist(int(sources[0]))
        return dijkstra(self.A, directed=False, indices=sources,
                        unweighted=True, min_only=True).astype(np.int32)


def reveal_rounds_for(max_rounds, explicit=None):
    if explicit:
        return sorted(set(explicit))
    return sorted({max(1, round(max_rounds * r / ORIGINAL_ROUNDS)) for r in ORIGINAL_REVEALS})


# ---- Mr. X strategies -------------------------------------------------------
def x_random(board, x, dets, rng):
    opts = [v for v in board.adj[x] if v not in dets]
    return rng.choice(opts) if opts else None


def x_evade(board, x, dets, rng):
    """Move to the neighbor farthest from the nearest detective.

    Uses a single multi-source shortest-path pass, which scales to full-city
    networks much better than calculating a separate distance field per detective.
    """
    opts = [v for v in board.adj[x] if v not in dets]
    if not opts:
        return None
    nearest = board.distance_to_any(dets)
    best_dist = max(int(nearest[v]) for v in opts)
    candidates = [v for v in opts if int(nearest[v]) == best_dist]
    best_deg = max(int(board.deg[v]) for v in candidates)
    return rng.choice([v for v in candidates if int(board.deg[v]) == best_deg])


# ---- Detective strategies ---------------------------------------------------
def det_random(board, i, dets, belief, last_seen, rng):
    others = set(dets[:i] + dets[i + 1:])
    opts = [v for v in board.adj[dets[i]] if v not in others]
    return rng.choice(opts) if opts else dets[i]


def det_belief(board, i, dets, belief, last_seen, rng):
    """Scale-friendly pursuit policy.

    Once Mr. X has been spotted, move along a shortest path toward the last
    reported location. Before the first sighting, favor connected high-degree
    junctions to explore the network. This policy intentionally avoids repeated
    all-candidate distance calculations that make city-sized games impractical.
    """
    others = set(dets[:i] + dets[i + 1:])
    opts = [v for v in board.adj[dets[i]] if v not in others]
    if not opts:
        return dets[i]
    if last_seen is not None:
        d = board.dist(last_seen)
        best = min(int(d[v]) for v in opts)
        candidates = [v for v in opts if int(d[v]) == best]
        # When equally near, prefer junctions with more onward choices.
        best_deg = max(int(board.deg[v]) for v in candidates)
        return rng.choice([v for v in candidates if int(board.deg[v]) == best_deg])
    best_deg = max(int(board.deg[v]) for v in opts)
    return rng.choice([v for v in opts if int(board.deg[v]) == best_deg])


X_STRATS = {"strategic": x_evade, "random": x_random}
D_STRATS = {"strategic": det_belief, "random": det_random}


def play_game(board, cfg, seed, starts=None, record=False):
    rng = random.Random(seed)
    n = board.n
    k = cfg["detectives"]
    T = cfg["max_rounds"]
    reveals = set(cfg["reveals"])
    xs, ds = X_STRATS[cfg["x_strategy"]], D_STRATS[cfg["det_strategy"]]

    pool = starts if starts else list(range(n))
    spawn = rng.sample(pool, k + 1)
    x, dets = spawn[0], spawn[1:]
    x_start = x

    belief_mask = np.ones(n, dtype=bool)
    belief_mask[dets] = False
    last_seen = None
    track = {"x": [x], "dets": [list(dets)], "revealed": [], "belief": []} if record else None
    collect = bool(cfg.get("collect"))
    x_path, det_visits = [x], []
    caught_round, caught_node = None, None
    outcome_reason = "survived"

    for t in range(1, T + 1):
        # --- Mr. X moves
        nx_ = xs(board, x, dets, rng)
        if nx_ is None:                              # boxed in by detectives
            caught_round, caught_node, outcome_reason = t, x, "no_move"
            break
        x = nx_
        x_path.append(x)
        # --- detectives' information
        if t in reveals:
            belief_mask[:] = False
            belief_mask[x] = True
            last_seen = x
            if record:
                track["revealed"].append(t)
        else:
            belief_mask = (board.A @ belief_mask.astype(np.int8)) > 0
            belief_mask[dets] = False
        belief_nodes = np.flatnonzero(belief_mask)
        if belief_nodes.size == 0:                   # should not happen; reset
            belief_mask[:] = True
            belief_mask[dets] = False
            belief_nodes = np.flatnonzero(belief_mask)
        belief = {"all": belief_nodes, "remaining": belief_nodes.copy()}
        if record:
            # A full-city belief set can contain tens of thousands of nodes per
            # round. It is not legible in the animation and bloats the response.
            track["belief"].append(belief_nodes.tolist() if n <= 10000 else [])
        # --- detectives move one by one (cannot share a square)
        for i in range(k):
            dets[i] = ds(board, i, dets, belief, last_seen, rng)
            if dets[i] == x:
                caught_round, caught_node, outcome_reason = t, x, "captured"
                break
        if record:
            track["x"].append(x)
            track["dets"].append(list(dets))
        if collect:
            det_visits.extend(dets)
        if caught_round:
            break
        belief_mask[dets] = False

    res = {
        "x_wins": int(caught_round is None),
        "rounds": caught_round if caught_round else T,
        "reason": outcome_reason,
        "x_start": x_start,
        "caught_node": caught_node if caught_node is not None else -1,
        "n_reveals_before_end": sum(1 for r in reveals if r <= (caught_round or T)),
    }
    if collect:
        res["x_path"] = x_path
        res["det_visits"] = det_visits
    if record:
        res["track"] = track
    return res


# ----------------------------------------------------------------------------
# Batch running
# ----------------------------------------------------------------------------
_B = {}


def _init_worker(edges, n, cfg, starts):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    _B["board"], _B["cfg"], _B["starts"] = Board(G), cfg, starts


def _run_one(seed):
    return play_game(_B["board"], _B["cfg"], seed, _B["starts"])


def wilson(p, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph", default="taxi",
                    help="taxi | er | ws | ba | edgelist:PATH | osm:\"Place, Country\"")
    ap.add_argument("--n", type=int, default=199, help="nodes (er/ws/ba)")
    ap.add_argument("--mean-degree", type=float, default=3.31, help="ER mean degree")
    ap.add_argument("--k", type=int, default=4, help="WS neighbours (even)")
    ap.add_argument("--p", type=float, default=0.1, help="WS rewiring prob")
    ap.add_argument("--m", type=int, default=2, help="BA edges per new node")
    ap.add_argument("--subgraph", type=int, default=0,
                    help="cut a connected ~N-node piece out of an edgelist/osm graph")
    ap.add_argument("--instances", type=int, default=1,
                    help="number of graph instances (for random/subgraph models)")
    ap.add_argument("--games", type=int, default=1000, help="games per instance")
    ap.add_argument("--detectives", type=int, default=6)
    ap.add_argument("--max-rounds", type=int, default=ORIGINAL_ROUNDS)
    ap.add_argument("--reveals", type=int, nargs="*", help="override reveal rounds")
    ap.add_argument("--x-strategy", choices=X_STRATS, default="strategic")
    ap.add_argument("--detectives-strategy", choices=D_STRATS, default="strategic")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results")
    ap.add_argument("--tag", default="", help="suffix for output file names")
    ap.add_argument("--replay", help="also save one recorded game to this JSON file")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    cfg = {"detectives": args.detectives, "max_rounds": args.max_rounds,
           "reveals": reveal_rounds_for(args.max_rounds, args.reveals),
           "x_strategy": args.x_strategy, "det_strategy": args.detectives_strategy}
    tag = f"_{args.tag}" if args.tag else ""
    games_path = os.path.join(args.out, f"games{tag}.csv")
    metrics_path = os.path.join(args.out, f"graph_metrics{tag}.csv")
    t0 = time.time()
    all_rows, metric_rows, name = [], [], None

    for inst in range(args.instances):
        G, name, label_starts = build_graph(args, inst)
        if G.number_of_nodes() < args.detectives + 2:
            sys.exit("graph too small for that many detectives")
        H, labels = relabel_zero_based(G)
        starts = None
        if label_starts:
            idx = {l: i for i, l in enumerate(labels)}
            starts = [idx[s] for s in label_starts if s in idx]
            if len(starts) < args.detectives + 1:
                starts = None
        met = graph_metrics(H, random.Random(inst))
        met.update({"graph": name, "instance": inst})
        metric_rows.append(met)
        edges = list(H.edges())
        seeds = [args.seed * 10_000_019 + inst * 1_000_003 + g for g in range(args.games)]
        with Pool(args.workers, initializer=_init_worker,
                  initargs=(edges, H.number_of_nodes(), cfg, starts)) as pool:
            res = pool.map(_run_one, seeds, chunksize=max(1, args.games // (args.workers * 8)))
        for g, r in enumerate(res):
            r.update({"graph": name, "instance": inst, "game": g})
            all_rows.append(r)
        wins = sum(r["x_wins"] for r in res)
        print(f"[{name} #{inst}] n={met['n']} m={met['m']}  Mr. X wins {wins}/{len(res)} "
              f"({wins/len(res):.1%})  mean rounds {np.mean([r['rounds'] for r in res]):.2f}")

        if args.replay and inst == 0:
            rec = play_game(Board(H), cfg, seeds[0], starts, record=True)
            edges_out = [[labels[a], labels[b]] for a, b in edges]
            rec["track"] = {k: ([labels[v] for v in vals] if k == "x" else
                               [[labels[v] for v in step] for step in vals] if k == "dets" else vals)
                            for k, vals in rec["track"].items()}
            with open(args.replay, "w") as f:
                json.dump({"graph": name, "edges": edges_out, "config": cfg, "game": rec}, f)
            print("replay saved to", args.replay)

    fields = ["graph", "instance", "game", "x_wins", "rounds", "reason", "x_start",
              "caught_node", "n_reveals_before_end"]
    with open(games_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows([{k: r[k] for k in fields} for r in all_rows])
    with open(metrics_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(metric_rows[0].keys()))
        w.writeheader()
        w.writerows(metric_rows)

    N = len(all_rows)
    wins = sum(r["x_wins"] for r in all_rows)
    lo, hi = wilson(wins / N, N)
    summary = {"graph": name, "games_total": N, "x_win_rate": wins / N,
               "x_win_rate_95ci": [lo, hi], "detective_win_rate": 1 - wins / N,
               "mean_rounds": float(np.mean([r["rounds"] for r in all_rows])),
               "config": cfg, "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(args.out, f"summary{tag}.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
