"""
Export all analysis outputs as plain CSVs + README into share/, then zip it.

Produces share/ with:
  nodes_metrics.csv        osmid, lat, lon, degree, betweenness, eccentricity_hops
  edges_metrics.csv        u, v, key, length_m, street_name, is_bridge, betweenness_mean
  percolation_curves.csv   q + giant-component share per removal strategy
  degree_distributions.csv degree, count_bengaluru, count_er, count_ws
  model_comparison.csv     summary metrics table (copy)
  top_betweenness.csv      top-20 intersections (copy)
  reports/*.txt            all text reports (copies)
  README.md                what each file is, methods, provenance
Then zips everything to bengaluru_network_analysis.zip.

Usage:  python export_share_data.py
"""
import csv
import gc
import os
import shutil
import time
import zipfile

import networkx as nx
import numpy as np
import osmnx as ox

GRAPH = "data/bengaluru_city_bengaluru_drive.graphml"
TAG = "bengaluru_city_bengaluru_drive"
OUT = "share"
os.makedirs(f"{OUT}/reports", exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
U = ox.convert.to_undirected(G)
del G; gc.collect()
S = nx.Graph(U)
S.remove_edges_from(nx.selfloop_edges(S))

z = np.load(f"results/{TAG}_node_bc.npz", allow_pickle=True)
bc = dict(zip(z["nodes"].tolist(), z["bc"]))
z = np.load(f"results/{TAG}_eccentricity.npz", allow_pickle=True)
ecc = dict(zip(z["nodes"].tolist(), z["ecc_lb"]))

# ---------------------------------------------------------------- nodes
log("writing nodes_metrics.csv ...")
with open(f"{OUT}/nodes_metrics.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["osmid", "lat", "lon", "degree", "betweenness", "eccentricity_hops"])
    for u, d in U.nodes(data=True):
        b = bc.get(u)
        e = ecc.get(u)
        w.writerow([u, f"{d['y']:.7f}", f"{d['x']:.7f}", S.degree(u),
                    f"{b:.8f}" if b is not None else "", e if e is not None else ""])

# ---------------------------------------------------------------- edges
log("writing edges_metrics.csv (recomputing bridges) ...")
multi_uv = {frozenset((u, v)) for u, v, k in U.edges(keys=True) if k > 0}
bridges = {frozenset(e) for e in nx.bridges(S)} - multi_uv
with open(f"{OUT}/edges_metrics.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["osmid_u", "osmid_v", "key", "length_m", "street_name",
                "is_bridge", "betweenness_mean"])
    for u, v, k, d in U.edges(keys=True, data=True):
        nm = d.get("name", "")
        if isinstance(nm, list):
            nm = "; ".join(map(str, nm))
        bm = (bc.get(u, 0.0) + bc.get(v, 0.0)) / 2
        w.writerow([u, v, k, f"{float(d.get('length', 0)):.1f}", nm,
                    int(frozenset((u, v)) in bridges), f"{bm:.8f}"])
del U, S; gc.collect()

# ---------------------------------------------------------------- percolation
log("writing percolation_curves.csv ...")
z = np.load(f"results/{TAG}_percolation.npz")
keys = [k for k in z.files if k != "q"]
with open(f"{OUT}/percolation_curves.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["q_fraction_removed"] + [f"S_{k}" for k in keys])
    q = z["q"]
    cols = [z[k] for k in keys]
    for i in range(len(q)):
        w.writerow([f"{q[i]:.6f}"] + [f"{c[i]:.6f}" for c in cols])

# ---------------------------------------------------------------- degree hists
log("writing degree_distributions.csv ...")
z = np.load("results/model_degrees.npz")
names = list(z.files)
dmax = max(len(z[g]) for g in names) - 1
with open(f"{OUT}/degree_distributions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["degree"] + [f"count_{g.lower()}" for g in names])
    for d in range(dmax + 1):
        w.writerow([d] + [int(z[g][d]) if d < len(z[g]) else 0 for g in names])

# ---------------------------------------------------------------- copies
log("copying tables and reports ...")
shutil.copy("results/model_comparison.csv", f"{OUT}/model_comparison.csv")
shutil.copy(f"results/{TAG}_top_betweenness.csv", f"{OUT}/top_betweenness.csv")
for fn in os.listdir("results"):
    if fn.endswith(".txt"):
        shutil.copy(f"results/{fn}", f"{OUT}/reports/{fn}")
shutil.copy("ANALYSIS.md", f"{OUT}/ANALYSIS.md")

# ---------------------------------------------------------------- README
with open(f"{OUT}/README.md", "w") as f:
    f.write("""# Bengaluru road network — analysis data (ID5080 term project)

Network: OpenStreetMap drivable streets, BBMP boundary, fetched 2026-10-05 with
osmnx 2.0.7 (`fetch_namma.py bengaluru city`). 155,779 nodes / 205,640 undirected
segments / 13,001 km. The 145 MB graphml itself is NOT in this bundle — regenerate
it with the script (OSM data is cached) or ask Kanak.

## Files
- **nodes_metrics.csv** — one row per intersection/dead-end. `degree` is the simple
  undirected degree. `betweenness` is normalized node betweenness, length-weighted,
  approximated with 300 sampled sources (seed 42) on the largest SCC (empty = node
  outside the SCC). `eccentricity_hops` = max BFS depth from the node (64-landmark
  lower bound; it matched the exact diameter of 283, so values are essentially exact).
- **edges_metrics.csv** — one row per street segment (parallel edges share u,v with
  key 0,1,...). `is_bridge`=1 means the segment lies on no cycle (removing it
  disconnects something). `betweenness_mean` = mean of endpoint node betweenness.
- **percolation_curves.csv** — giant-component share S(q) vs fraction q of nodes
  removed, exact at every q (reverse Newman–Ziff union–find). Strategies: random
  (5-seed average), static degree, adaptive degree, static betweenness.
- **degree_distributions.csv** — degree histograms for Bengaluru and the matched
  ER G(n,m) and WS (k=4, p=0.1, seed 42) models (counts; normalize by column sum).
- **model_comparison.csv** — summary metrics: clustering (20k-trial approximation),
  sampled mean path length and diameter lower bound (48 BFS sources), spectral gap
  lambda_2 of the normalized Laplacian, assortativity, bridge fraction, cyclomatic
  number, alpha index.
- **top_betweenness.csv** — the 20 highest-betweenness intersections with
  coordinates and street names.
- **reports/** — plain-text reports as produced by each analysis script, including
  exact diameter (283 hops), Cheeger/Ramanujan verdict, min-cut sizes (N–S: 29,
  E–W: 24), gamma (power-law) check, and percolation thresholds.
- **ANALYSIS.md** — the write-up tying everything together (figures referenced
  there live in the project repo, not in this bundle).

## Key caveats
- Betweenness is a k=300 sampled approximation: top ranks are stable, exact values
  are not; rerun with python-igraph for exact numbers if needed.
- Model WS uses k=4 (integer constraint), so it has ~52% more edges than Bengaluru.
- Everything derives from OSM data (ODbL license — attribute "© OpenStreetMap
  contributors" in anything published).
""")

# ---------------------------------------------------------------- zip
log("zipping ...")
zpath = "bengaluru_network_analysis.zip"
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as zf:
    for root, _, files in os.walk(OUT):
        for fn in files:
            p = os.path.join(root, fn)
            zf.write(p, os.path.relpath(p, "."))
log(f"done: {zpath} ({os.path.getsize(zpath)/1e6:.1f} MB) + {OUT}/ folder")
