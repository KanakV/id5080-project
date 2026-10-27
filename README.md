# Namma Graph — street network analysis (ID5080 term project)

Complex-networks analysis of city street networks. This repo contains the
**Bengaluru** analysis end-to-end; every script is city-agnostic, so the same
pipeline runs unchanged for **Chennai** (or any city OSM knows).

- Results write-up with all figures: **[ANALYSIS.md](ANALYSIS.md)**
- How to add a city or a new metric: **[CONTRIBUTIONS.md](CONTRIBUTIONS.md)**

## Where the data comes from

All networks are from **OpenStreetMap** (© OpenStreetMap contributors, ODbL —
keep this credit on anything published). `fetch_namma.py` uses **osmnx 2.x**,
which geocodes the municipal boundary polygon through **Nominatim** (e.g.
"Bengaluru, Karnataka, India") and downloads the drivable street network inside
it from the **Overpass API** (`network_type="drive"`, topologically simplified:
nodes = intersections/dead-ends, edges = street segments with length in metres).

To download (one-time internet need; responses are cached in `cache/`, so
re-runs are offline and free):

```bash
.venv/bin/python fetch_namma.py bengaluru city     # or: chennai
# -> data/bengaluru_city_bengaluru_drive.graphml  (~145 MB, 5–15 min; cache/ ~44 MB)
```

Heads-up: OSM is live data — a fetch today differs slightly from a fetch last
month. The team should analyse **one shared graphml per city** (the Bengaluru
one here was fetched 2026-10-05); don't re-fetch mid-analysis. `data/` and
`cache/` are git-ignored, so pass the graphml around directly if you want
bit-identical inputs.

## Quick start (plug and play)

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt     # osmnx 2.x, networkx, scipy, matplotlib

.venv/bin/python run_all.py chennai           # fetch (if needed) + every analysis step
.venv/bin/python run_all.py chennai --steps core,maps   # just some steps
.venv/bin/python run_all.py --list            # what the steps are
```

All outputs are prefixed with the graph's tag (`chennai_city_chennai_drive_...`),
so cities never overwrite each other. Expect ~45–60 min total for a 150k-node
city; the slow parts are betweenness (~12 min) and the exact-diameter
confirmation (~25 min, runs last in the `depth` step).

## The steps (also runnable by hand)

| Step | Script | What it computes | Needs |
|---|---|---|---|
| `core` | `analyze_network.py <graphml> 300` | Betweenness centrality (length-weighted, 300 sampled sources) + cyclic structure: cyclomatic number, α/β/γ, bridges, 2-core, biconnected blocks | — |
| `expander` | `expander_bottleneck.py` | Spectral gap λ₂, Cheeger bounds, Ramanujan/expander verdict, Fiedler sweep cut, N–S and E–W max-flow min-cuts | — |
| `models` | `model_comparison.py` | Same-size ER G(n,m) and WS (k=4, p=0.1) null models: clustering, paths, λ₂, assortativity, bridges | — |
| `depth` | `depth_diameter.py` | Per-node eccentricity (64-landmark BFS), exact diameter, power-law γ check | — |
| `percolation` | `percolation.py` | Robustness S(q): random failure vs degree/betweenness attacks (exact reverse union–find curves) | core |
| `perc-models` | `percolation_models.py` | Percolation vs ER and WS, dark chart | models, percolation |
| `degrees` | `degree_distribution.py` | Degree histograms (raw + normalized), dark chart | models |
| `maps` | `render_dark_maps.py` | Dark-theme betweenness / cycles / eccentricity maps | core, depth |
| `edges-map` | `betweenness_edges.py` | Light-theme edge-painted betweenness map | core |
| `assort` | `degree_correlations.py` | Assortativity r, e_jk, k_nn(k), local assortativity maps, road-type breakdown, boundary clusters (+ cache npz); compare cities with `assortativity_compare.py <graphml> <graphml> ...`, land use at clusters with `cluster_landuse.py` | core |
| `hubs` | `hub_bottlenecks.py` | Named hub pairs: effective resistance, corridor Fiedler λ₂ / conductance, min cut between hub areas | — |
| `edits` | `edit_distance.py` | Edges to add / remove to reach Barcelona, ER, scale-free and ring–spoke robustness & flow (exact degree lower bound + greedy upper bound); needs `data/barcelona_city_barcelona_drive.graphml` (`fetch_namma.py barcelona city`) | assort |
| `export` | `export_share_data.py` | Plain-CSV hand-off bundle (`share_<city>/` + zip) for non-Python teammates | all |

Every script takes the graphml path as its first argument and caches expensive
arrays in `results/*.npz`, so re-rendering figures costs seconds, not minutes.

## Gotchas we already hit (so you don't)

- **Nominatim neighbourhood queries fail** for some places ("did not geocode to
  (Multi)Polygon") — those OSM entries are points, not polygons. Use
  `ox.graph_from_address(..., dist=...)` for neighbourhood-scale graphs.
- **Self-loops**: real city graphs have them and `nx.k_core` / `nx.bridges`
  refuse them — strip with `G.remove_edges_from(nx.selfloop_edges(G))` on the
  simple graph (already handled in these scripts).
- **Betweenness is approximate** (k=300 sampled sources, seed 42): rankings are
  stable, exact values are not. For exact numbers use python-igraph
  (~50–100× faster than networkx).
- **WS density**: Watts–Strogatz needs an even integer k, so k=4 is ~52% denser
  than a street network (d̄ ≈ 2.6) — flagged wherever it matters.
- **λ₂**: shift-invert Lanczos is cheap on near-planar graphs; for random
  models use Lanczos on 2I − L instead (no factorization). Both patterns are in
  `expander_bottleneck.py` / `model_comparison.py`.

## Repo layout

- `data/`, `cache/` — git-ignored (graphml ~145 MB; regenerate via fetch, see above)
- `results/` — committed: reports (txt), tables (csv), small npz caches
- `figures/` — committed: rendered maps and charts (`*_dark` = dark theme)
- `share_*/`, `*.zip` — git-ignored hand-off bundles (`export` step rebuilds them)
