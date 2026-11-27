#!/usr/bin/env python3
"""
Visual road-only Scotland Yard.  Run:  python app.py   (opens http://127.0.0.1:5000)

Main menu -> pick a map -> watch one game, or run hundreds/thousands of games and
explore the analysis.  Game logic lives in sy_sim.py (same engine as the CLI).
"""
import csv
import glob
import io
import os
import random
import socket
import threading
import webbrowser
from argparse import Namespace
from multiprocessing import Pool, freeze_support

import networkx as nx
import numpy as np
from flask import Flask, Response, jsonify, request

import sy_sim as sy

ROOT = os.path.dirname(os.path.abspath(__file__))
MAPS_DIR = os.path.join(ROOT, "maps")
app = Flask(__name__, static_folder=os.path.join(ROOT, "static"), static_url_path="/static")

STATE = {"g": None, "batch": None}
LOCK = threading.Lock()

BUILTIN = [
    {"id": "taxi", "name": "Original board", "random": False,
     "desc": "Scotland Yard's 199 stations, taxi (road) edges only."},
    {"id": "er", "name": "Erdős–Rényi", "random": True,
     "desc": "Random graph: every pair of nodes linked independently."},
    {"id": "ws", "name": "Watts–Strogatz", "random": True,
     "desc": "Small world: a ring lattice with a few rewired shortcuts."},
    {"id": "ba", "name": "Power-law (Barabási–Albert)", "random": True,
     "desc": "Preferential attachment: a few huge hubs, many leaves."},
]


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


app.json.default = _json_default


CITY_CACHE = {}


def pretty_name(stem):
    parts = stem.split("_")
    city = parts[0].title()
    scale = parts[1] if len(parts) > 1 else ""
    rest = " ".join(parts[2:]).title()
    if scale == "city":
        return f"{city} (whole city)"
    if scale == "nbhd":
        return f"{city}: {rest}"
    if scale == "metro":
        return f"{city} (30 km metro area)"
    return stem.replace("_", " ").title()


def list_maps():
    out = list(BUILTIN)
    for f in sorted(glob.glob(os.path.join(MAPS_DIR, "*_edges.csv"))):
        stem = os.path.basename(f)[:-len("_edges.csv")]
        has_xy = os.path.exists(os.path.join(MAPS_DIR, stem + "_nodes.csv"))
        out.append({"id": "file:" + stem, "name": pretty_name(stem), "random": False,
                    "big": True, "whole_network": True,
                    "desc": "Entire connected road network, with intersections placed at their map coordinates."
                    if has_xy else "Entire road network (no coordinates, so it is laid out automatically)."})
    return out


def read_coords(stem):
    path = os.path.join(MAPS_DIR, stem + "_nodes.csv")
    if not os.path.exists(path):
        return None
    d = {}
    with open(path, newline="") as f:
        for row in csv.reader(f):
            try:
                d[int(row[0])] = (float(row[1]), float(row[2]))      # lon, lat
            except (ValueError, IndexError):
                continue
    return d


def load_city(stem):
    if stem not in CITY_CACHE:
        G = sy.largest_component(sy.read_edge_file(os.path.join(MAPS_DIR, stem + "_edges.csv")))
        coords = read_coords(stem)
        CITY_CACHE[stem] = (G, coords if coords else None)
    return CITY_CACHE[stem]


def normalise(pos):
    P = np.array(pos, dtype=float)
    P -= P.min(axis=0)
    span = P.max(axis=0).max() or 1.0
    P /= span
    P += (1 - P.max(axis=0)) / 2          # centre, keeping the true aspect ratio
    return P


def board_geometry():
    """True station positions on the board picture (pixels of the cropped image)."""
    pos, size = {}, (1510, 1130)
    with open(os.path.join(sy.HERE, "data", "board_nodes.csv"), newline="") as f:
        for r in csv.DictReader(f):
            pos[int(r["label"])] = (float(r["x_px"]), float(r["y_px"]))
            size = (int(r["image_width"]), int(r["image_height"]))
    return pos, size


def make_graph(map_id, n, seed, subgraph, area=""):
    starts_labels, area_label = None, ""
    box, bg, geo_pos = [1.0, 1.0], None, None
    if map_id == "taxi":
        G = sy.read_edge_file(os.path.join(sy.HERE, "data", "board_edges.csv"))
        starts_labels, name = sy.ORIGINAL_STARTS, "Original Scotland Yard board"
        px, (iw, ih) = board_geometry()
        geo_pos = {l: (px[l][0] / iw, px[l][1] / iw) for l in G.nodes()}
        box = [1.0, ih / iw]
        bg = {"url": "/static/boards/scotland_yard.jpg", "aspect": iw / ih}
    elif map_id in ("er", "ws", "ba"):
        a = Namespace(graph=map_id, n=n, mean_degree=3.31, k=4, p=0.1, m=2, seed=seed, subgraph=0)
        G, _, _ = sy.build_graph(a, 0)
        name = next(m["name"] for m in BUILTIN if m["id"] == map_id)
    elif map_id.startswith("file:"):
        stem = map_id[5:]
        full, coords = load_city(stem)
        name = pretty_name(stem)
        # Real city maps are always played as their complete largest connected
        # road network. The old subgraph/area parameters are accepted by the API
        # for backward compatibility, but deliberately ignored.
        G = full
        if coords and all(i in coords for i in G.nodes()):
            lat0 = np.mean([coords[i][1] for i in G.nodes()])
            k = float(np.cos(np.radians(lat0)))
            geo_pos = {i: ((coords[i][0]) * k, -coords[i][1]) for i in G.nodes()}
    else:
        raise ValueError("unknown map")
    n_nodes = G.number_of_nodes()
    # The city converter already relabels nodes to 0..n-1. Reusing that graph
    # avoids a costly second NetworkX copy for 70k--160k node city maps.
    if n_nodes and all(isinstance(v, (int, np.integer)) for v in G.nodes()) and min(G.nodes()) == 0 and max(G.nodes()) == n_nodes - 1:
        H, labels = G, list(range(n_nodes))
    else:
        H, labels = sy.relabel_zero_based(G)
    if map_id == "taxi":
        P = np.array([geo_pos[l] for l in labels])
    elif geo_pos:
        P = normalise([geo_pos[l] for l in labels])
    elif H.number_of_nodes() <= 400:
        L = nx.kamada_kawai_layout(H)
        P = normalise([L[i] for i in range(H.number_of_nodes())])
    else:
        L = nx.spring_layout(H, seed=1, iterations=40)
        P = normalise([L[i] for i in range(H.number_of_nodes())])
    starts = None
    if starts_labels:
        idx = {l: i for i, l in enumerate(labels)}
        starts = [idx[s] for s in starts_labels if s in idx]
    metrics = sy.graph_metrics(H, random.Random(0))
    return {
        "H": H, "labels": labels, "starts": starts, "board": sy.Board(H), "name": name,
        "map_id": map_id, "x": P[:, 0], "y": P[:, 1], "metrics": metrics, "node_metrics": None,
        "box": box, "bg": bg, "area": area_label,
    }


def graph_payload(g):
    n = g["H"].number_of_nodes()
    diam = g["metrics"]["diameter"]
    return {
        "name": g["name"], "area": g["area"], "n": n, "labels": g["labels"],
        "x": [round(float(v), 5) for v in g["x"]], "y": [round(float(v), 5) for v in g["y"]],
        "edges": [list(e) for e in g["H"].edges()], "metrics": g["metrics"],
        "box": g["box"], "bg": g["bg"],
        "show_labels": bool(g["map_id"] == "taxi"),
        # City-wide games need more turns than the original 199-station board,
        # but an estimated graph diameter is only a sampled lower bound. Keep a
        # sensible default while allowing the UI to raise it up to 500 rounds.
        "suggested_rounds": (80 if n > 50000 else 60 if n > 10000 else
                             max(24, min(160, int(round(1.1 * diam))))),
        # Same city default for Chennai and Bengaluru so their results can be
        # compared under identical settings. Users can change this in the UI.
        "suggested_detectives": (30 if n > 10000 else 6),
        "whole_network": bool(g["map_id"].startswith("file:")),
    }


def cfg_from(d, g):
    n = g["H"].number_of_nodes()
    T = max(1, min(500, int(d.get("max_rounds", 24))))
    k = max(1, min(100, int(d.get("detectives", 6)), n - 2))
    xs = d.get("x_strategy", "strategic")
    ds = d.get("det_strategy", "strategic")
    if xs not in sy.X_STRATS or ds not in sy.D_STRATS:
        raise ValueError("bad strategy")
    return {"detectives": k, "max_rounds": T, "reveals": sy.reveal_rounds_for(T),
            "x_strategy": xs, "det_strategy": ds}


@app.route("/")
def index():
    return app.send_static_file("index.html")


@app.route("/api/maps")
def api_maps():
    return jsonify(list_maps())


@app.route("/api/graph", methods=["POST"])
def api_graph():
    d = request.get_json(force=True)
    try:
        g = make_graph(d["map"], int(d.get("n", 199)), int(d.get("seed", 1)),
                       int(d.get("subgraph", 0) or 0), d.get("area", ""))
    except Exception as e:                      # noqa
        return jsonify({"error": str(e)}), 400
    with LOCK:
        STATE["g"] = g
    return jsonify(graph_payload(g))


@app.route("/api/play", methods=["POST"])
def api_play():
    g = STATE["g"]
    if g is None:
        return jsonify({"error": "no graph loaded"}), 400
    d = request.get_json(force=True)
    cfg = cfg_from(d, g)
    seed = random.randrange(1, 10**9)
    rec = sy.play_game(g["board"], cfg, seed, g["starts"], record=True)
    return jsonify({"config": cfg, "result": {k: v for k, v in rec.items() if k != "track"},
                    "track": rec["track"]})


def node_metrics(g):
    if g["node_metrics"] is None:
        H, n = g["H"], g["H"].number_of_nodes()
        deg = np.fromiter((H.degree(i) for i in range(n)), dtype=float, count=n)
        if n <= 1500:
            bc = nx.betweenness_centrality(H)
            cc = nx.closeness_centrality(H)
        else:
            # Sampling avoids an excessively expensive centrality pass on the
            # 70k--160k node whole-city graphs.
            k = min(n, 10 if n > 50000 else 20 if n > 10000 else 300)
            bc = nx.betweenness_centrality(H, k=k, seed=1)
            cc = None
        g["node_metrics"] = {
            "degree": deg,
            "betweenness": np.array([bc[i] for i in range(n)]),
            "closeness": np.array([cc[i] for i in range(n)]) if cc else None,
        }
    return g["node_metrics"]


def spearman(a, b):
    from scipy.stats import spearmanr
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() < 5 or np.ptp(a[ok]) == 0 or np.ptp(b[ok]) == 0:
        return None
    r = spearmanr(a[ok], b[ok]).statistic
    return None if np.isnan(r) else round(float(r), 3)


def analyse(results, g, cfg):
    n, T = g["H"].number_of_nodes(), cfg["max_rounds"]
    N = len(results)
    wins = np.array([r["x_wins"] for r in results])
    rounds = np.array([r["rounds"] for r in results])
    x_vis = np.bincount(np.concatenate([r["x_path"] for r in results]), minlength=n).astype(float)
    d_vis = np.bincount(np.concatenate([r["det_visits"] or [0] for r in results]), minlength=n).astype(float)
    caught = np.zeros(n)
    for r in results:
        if not r["x_wins"] and r["caught_node"] >= 0:
            caught[r["caught_node"]] += 1
    rate = np.where(x_vis >= 20, caught / np.maximum(x_vis, 1), np.nan)
    p = wins.mean()
    lo, hi = sy.wilson(p, N)
    cap_hist = [int(((rounds == t) & (wins == 0)).sum()) for t in range(1, T + 1)]
    surv = [float(((wins == 1) | (rounds > t)).mean()) for t in range(0, T + 1)]
    nm = node_metrics(g)
    corr = {}
    for mname, arr in nm.items():
        if arr is None:
            continue
        corr[mname] = {"x_visits": spearman(arr, x_vis), "captures": spearman(arr, caught),
                       "capture_rate": spearman(arr, rate), "det_visits": spearman(arr, d_vis)}
    deg_start = {}
    deg = nm["degree"]
    for r in results:
        k = int(deg[r["x_start"]])
        a = deg_start.setdefault(k, [0, 0])
        a[0] += 1
        a[1] += r["x_wins"]
    start_by_degree = [{"degree": k, "games": v[0], "x_win_rate": v[1] / v[0]}
                       for k, v in sorted(deg_start.items())]
    order = np.argsort(-caught)[:8]
    top = [{"label": g["labels"][i], "captures": int(caught[i]), "degree": int(deg[i]),
            "x_visits": int(x_vis[i])} for i in order if caught[i] > 0]
    return {
        "games": N, "x_win_rate": float(p), "ci": [lo, hi], "mean_rounds": float(rounds.mean()),
        "median_rounds": float(np.median(rounds)), "config": cfg, "captured_by_round": cap_hist,
        "escaped": int(wins.sum()), "survival": surv, "x_visits": x_vis, "det_visits": d_vis,
        "captures": caught, "capture_rate": [None if np.isnan(v) else float(v) for v in rate],
        "correlations": corr, "start_by_degree": start_by_degree, "top_capture_nodes": top,
    }


def run_batch(games, cfg, g, batch):
    try:
        H = g["H"]
        edges = list(H.edges())
        cfg = dict(cfg, collect=True)
        base = random.randrange(1, 10**6)
        seeds = [base * 100003 + i for i in range(games)]
        n = H.number_of_nodes()
        available = max(1, (os.cpu_count() or 2) - 1)
        # Each worker reconstructs the road graph and sparse distance matrix.
        # More processes can exhaust RAM on full-city maps, so cap concurrency.
        workers = 1 if n > 50000 else min(2, available) if n > 10000 else available
        workers = min(workers, max(1, games))
        results = []
        with Pool(workers, initializer=sy._init_worker,
                  initargs=(edges, H.number_of_nodes(), cfg, g["starts"])) as pool:
            for r in pool.imap_unordered(sy._run_one, seeds, chunksize=max(1, games // (workers * 20))):
                results.append(r)
                batch["done"] = len(results)
        batch["rows"] = [(i, r["x_wins"], r["rounds"], r["reason"], g["labels"][r["x_start"]],
                          g["labels"][r["caught_node"]] if r["caught_node"] >= 0 else "")
                         for i, r in enumerate(results)]
        batch["result"] = analyse(results, g, cfg)
        batch["status"] = "done"
    except Exception as e:                      # noqa
        batch["status"], batch["error"] = "error", repr(e)


@app.route("/api/batch/start", methods=["POST"])
def api_batch_start():
    g = STATE["g"]
    if g is None:
        return jsonify({"error": "no graph loaded"}), 400
    d = request.get_json(force=True)
    cfg = cfg_from(d, g)
    games = max(10, min(50000, int(d.get("games", 1000))))
    batch = {"status": "running", "done": 0, "total": games, "result": None}
    STATE["batch"] = batch
    threading.Thread(target=run_batch, args=(games, cfg, g, batch), daemon=True).start()
    return jsonify({"ok": True})


@app.route("/api/batch/status")
def api_batch_status():
    b = STATE["batch"]
    if b is None:
        return jsonify({"status": "none"})
    out = {"status": b["status"], "done": b["done"], "total": b["total"]}
    if b["status"] == "done":
        out["result"] = b["result"]
        out["graph_name"] = STATE["g"]["name"]
    if b["status"] == "error":
        out["error"] = b.get("error")
    return jsonify(out)


@app.route("/api/batch/games.csv")
def api_batch_csv():
    b = STATE["batch"]
    if not b or b["status"] != "done":
        return "no finished batch", 404
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["game", "x_wins", "rounds", "reason", "x_start", "caught_node"])
    w.writerows(b["rows"])
    return Response(buf.getvalue(), mimetype="text/csv",
                    headers={"Content-Disposition": "attachment; filename=games.csv"})


def free_port(start=5000):
    for p in range(start, start + 50):
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", p)) != 0:
                return p
    return start


if __name__ == "__main__":
    freeze_support()
    port = free_port()
    url = f"http://127.0.0.1:{port}"
    print("Open", url, "(Ctrl+C to quit)")
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(host="127.0.0.1", port=port, threaded=True, use_reloader=False)
