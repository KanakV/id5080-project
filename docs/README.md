# Aakash's contributions: Bengaluru vs Chennai street networks

Four analyses comparing the drive networks of Bengaluru and Chennai (OpenStreetMap). Each
one has a script in `tools/`, results in `results/` and figures in `figures/`.

| Analysis | Question | Script | Write-up |
|---|---|---|---|
| Degree assortativity | Do busy junctions connect to busy junctions, and where? | `degree_correlations.py`, `assortativity_compare.py`, `cluster_landuse.py` | [METHODS.md §2](METHODS.md) |
| Hub-to-hub bottlenecks | How well connected are the city's main activity centres to each other? | `hub_bottlenecks.py`, `hub_compare.py` | [METHODS.md §3](METHODS.md) |
| Edit distance | How many new local roads would make the centre robust? | `edit_spatial.py` | [EDIT_DISTANCE.md](EDIT_DISTANCE.md) |
| Conflict points | How much must drivers coordinate at a typical intersection? | `conflict_points.py` | [CONFLICT_POINTS.md](CONFLICT_POINTS.md) |

Discarded attempts are listed in [TRIED.md](TRIED.md).

## Results

### Degree assortativity

| | Bengaluru | Chennai | Barcelona (reference) |
|---|---|---|---|
| Assortativity r | −0.015 | −0.013 | +0.215 |
| Dead ends (degree 1) | 24.6% | 22.8% | 7.3% |
| 4-way junctions | 13.3% | 11.2% | 28.9% |

- **Both Indian cities are almost neutral (r ≈ 0).** Barcelona is clearly assortative.
- **Bengaluru shifts with distance from the centre.** Its core is assortative (local ρ·N
  +0.47 within 3 km) and its outskirts disassortative (−0.34 beyond 12 km). In the outskirts,
  arterials feed dead-end streets.
- **Chennai is flat at every distance.**

Figures: `figures/compare_assortativity_dark.png` and `figures/<city>_*assortativity*_dark.png`.

### Hub-to-hub bottlenecks (28 hub pairs in Bengaluru, 19 in Chennai, 4–30 km by road)

| Median over hub pairs | Bengaluru | Chennai |
|---|---|---|
| Resistance per road km (lower = more parallel routes) | 0.051 | 0.063 |
| Min cut between hub areas (higher = harder to cut off) | 21 | 9 |
| Corridor conductance φ ×10⁻³ (higher = no narrow neck) | 0.84 | 0.91 |

- **Chennai is easier to cut apart.** Its hubs can be cut off from each other by fewer than
  half as many road closures.
- **Chennai has fewer parallel routes.** Its resistance per road km is 25% higher.
- **Corridor necks are similar** in the two cities.
- **Hubs that are choke points in their own catchment:** Electronic City and Whitefield
  (Bengaluru), Sholinganallur (Chennai).

Figure: `figures/compare_hubs_dark.png`.

### Edit distance to a robust network

New local roads per 1,000 existing roads before the central 9,000 junctions reach 80% of a
perfect grid's robustness, under both random failure and targeted attack:

| City | Roads needed |
|---|---|
| Bengaluru | **239** |
| Chennai | **349** (1.5× more) |

The ranking holds when roads are added in random order (334 vs 407). Figure:
`figures/edit_spatial_dark.png`.

### Conflict points

| | Bengaluru | Chennai |
|---|---|---|
| Mean conflict points per intersection | 11.8 | 11.0 |
| Median | 9 (mostly T-junctions) | 9 |

Bengaluru's extra conflicts are crossing conflicts at bigger junctions.
Figure: `figures/conflict_points.png`.

## How to run

```bash
pip install osmnx networkx numpy scipy matplotlib
python tools/fetch_namma.py bengaluru city      # also: chennai city, barcelona city -> data/*.graphml
python tools/analyze_network.py data/bengaluru_city_bengaluru_drive.graphml 300   # top-betweenness CSV (defines the CBD)
python tools/degree_correlations.py data/bengaluru_city_bengaluru_drive.graphml
python tools/hub_bottlenecks.py data/bengaluru_city_bengaluru_drive.graphml
#   ...repeat the per-city steps for chennai (and barcelona for the assortativity reference)
python tools/assortativity_compare.py data/bengaluru_city_bengaluru_drive.graphml data/chennai_city_chennai_drive.graphml data/barcelona_city_barcelona_drive.graphml
python tools/cluster_landuse.py data/bengaluru_city_bengaluru_drive.graphml
python tools/hub_compare.py
python tools/edit_spatial.py
python tools/conflict_points.py
```

The `results/<city>_top_betweenness.csv` files are included, so the slow `analyze_network.py`
step can be skipped.

`tools/fetch_namma.py` and `tools/analyze_network.py` are copies of Kanak Varma's team
pipeline scripts, included only so this branch runs on its own.
