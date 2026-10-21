# Namma Graph — street network analysis (ID5080 term project)

Complex-networks analysis of city street networks from OpenStreetMap.
This repo contains the **Bengaluru** analysis end-to-end; every script is
city-agnostic (it takes a GraphML path), so the same pipeline runs unchanged
for **Chennai** — see [Running it for Chennai](#running-it-for-chennai).

Results write-up with all figures: **[ANALYSIS.md](ANALYSIS.md)**.

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt    # osmnx 2.x, networkx, scipy, matplotlib
```

## Pipeline

Step 0 downloads the network (needs internet; Overpass responses are cached in
`cache/`, so re-runs are free). Everything after runs offline on the saved graph.

| # | Command | What it does | Outputs |
|---|---|---|---|
| 0 | `python fetch_namma.py <city> city` | Fetch drivable network inside the municipal boundary | `data/<tag>.graphml` |
| 1 | `python analyze_bengaluru.py data/<tag>.graphml 300` | Betweenness centrality (length-weighted, 300 sampled sources) + cyclic structure: cyclomatic number, α/β/γ, bridges, 2-core, biconnected blocks | report, top-20 CSV, 2 maps |
| 2 | `python expander_bottleneck.py data/<tag>.graphml` | Spectral gap λ₂, Cheeger bounds, Ramanujan/expander verdict, Fiedler sweep cut, N–S and E–W max-flow min-cuts | report, 3 dark maps |
| 3 | `python model_comparison.py data/<tag>.graphml` | Same-size ER G(n,m) and WS (k=4, p=0.1) null models: clustering, paths, λ₂, assortativity, bridges | `results/model_comparison.csv`, degree histograms, metric panels |
| 4 | `python depth_diameter.py data/<tag>.graphml` | Eccentricity of every node (64-landmark BFS), exact diameter (extrema bounding — the slow last step, ~25 min), power-law γ check | report, eccentricity map + npz |
| 5 | `python percolation.py data/<tag>.graphml` | Robustness: S(q) under random failure vs degree/betweenness attacks (exact curves via reverse union–find) | report, npz, chart |
| 6 | `python percolation_models.py` | Percolation curves vs ER and WS (random failure) | dark chart |
| 7 | `python degree_distribution.py` | Degree histograms, raw + normalized (from step 3's cache) | dark chart |
| 8 | `python betweenness_edges.py data/<tag>.graphml` | Edge-painted betweenness map (light theme) | map |
| 9 | `python render_dark_maps.py data/<tag>.graphml` | Dark-theme betweenness / cycles / eccentricity maps | 3 maps |
| 10 | `python export_share_data.py` | Export everything as plain CSVs + README and zip it | `share/`, zip |

Notes on order: 1 must run before 5, 8, 9 (they reuse the cached betweenness
`results/<tag>_node_bc.npz`); 4 before 9 (eccentricity cache); 3 before 6 and 7.
Scripts 6, 7, 10 have Bengaluru names/sizes hardcoded near the top — edit the
constants for Chennai.

Timings on a 155k-node city (i7-class machine): step 0 ~5–15 min (Overpass),
step 1 ~12 min (betweenness is the slow part), step 2 ~3 min, step 4 ~5 min +
~25 min for the exact-diameter confirmation, the rest seconds-to-a-couple-minutes.

## Running it for Chennai

```bash
.venv/bin/python fetch_namma.py chennai city
.venv/bin/python analyze_bengaluru.py data/chennai_city_chennai_drive.graphml 300
# ...then steps 2, 4, 5, 8, 9 with that same path; edit constants in 6, 7, 10
```

## Gotchas we already hit (so you don't)

- **Nominatim neighbourhood queries fail** for some places ("did not geocode to
  (Multi)Polygon") — those OSM entries are points. Use
  `ox.graph_from_address(..., dist=...)` instead of `graph_from_place` for
  neighbourhood-scale graphs.
- **Self-loops**: city graphs contain them and `nx.k_core`/`nx.bridges` choke —
  strip with `G.remove_edges_from(nx.selfloop_edges(G))` on the simple graph
  (already handled in the scripts).
- **Betweenness is approximate** (k=300 sampled sources, seed 42): rankings are
  stable, exact values are not. For exact numbers use python-igraph (~50–100×
  faster than networkx).
- **WS density**: Watts–Strogatz needs an even integer k, so k=4 makes it ~52%
  denser than the street network (d̄ ≈ 2.6) — flagged wherever it matters.
- **λ₂ computation**: shift-invert Lanczos works great on near-planar graphs
  (cheap LU); for the random models use Lanczos on 2I − L instead (no
  factorization) — `expander_bottleneck.py`/`model_comparison.py` show both.
- **Charts/maps**: maps are dark-themed by preference; the colors come from a
  CVD-validated palette — keep series colors/orders as in the scripts.

## Repo layout

- `data/`, `cache/` — git-ignored (the Bengaluru graphml is 145 MB; regenerate
  with step 0, Overpass responses cache locally).
- `results/` — committed: reports (txt), tables (csv) and small npz caches,
  including node betweenness so you can re-render maps without the 7-min recompute.
- `figures/` — committed: all rendered maps and charts (`*_dark` = dark theme).
- `share/`, `*.zip` — git-ignored; `export_share_data.py` rebuilds the
  plain-CSV hand-off bundle in ~1 min.

## Attribution

Map data © OpenStreetMap contributors (ODbL) — keep this credit on anything published.
