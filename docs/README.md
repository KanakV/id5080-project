# Aakash's contributions: Bengaluru vs Chennai street networks

Two analyses compare the drive networks of Bengaluru and Chennai (OpenStreetMap) directly.
Each has its scripts in `tools/`, its results in `results/` and its figures in `figures/`.

| Analysis | Question | Scripts | Write-up |
|---|---|---|---|
| Degree assortativity | Do busy junctions connect to busy junctions, and where in the city? | `degree_correlations.py`, `assortativity_compare.py`, `cluster_landuse.py` | [METHODS.md §2](METHODS.md) |
| Hub-to-hub bottlenecks | How well are the city's main activity centres connected to each other? | `hub_bottlenecks.py`, `hub_compare.py` | [METHODS.md §3](METHODS.md) |

**Read with care.** The two cities are compared directly, with no reference city. Chennai is
coastal and cut by rivers; Bengaluru is inland. They also differ in size and boundary. So a
difference can come from geography as well as from planning; see
[METHODS.md §4](METHODS.md). Dropped attempts are listed in [TRIED.md](TRIED.md).

## Results

### Degree assortativity

| | Bengaluru | Chennai |
|---|---|---|
| Assortativity r | −0.015 | −0.013 |
| Dead ends (degree 1) | 24.6% | 22.8% |
| 4-way junctions | 13.3% | 11.2% |
| Mean local assortativity ρ·N, < 3 km from the CBD | +0.47 | +0.11 |
| Mean local assortativity ρ·N, > 12 km from the CBD | −0.34 | −0.02 |

- **Both cities are almost neutral overall (r ≈ 0).**
- **They differ in where mixing happens.** Bengaluru has an assortative core and
  disassortative outskirts, where arterials feed dead-end streets; the two cancel out in r.
- **Chennai is flat at every distance from the centre.**

Figure: `figures/compare_assortativity_dark.png`.

### Hub-to-hub bottlenecks

28 hub pairs in Bengaluru and 19 in Chennai, all 4–30 km apart by road:

| Median over hub pairs | Bengaluru | Chennai |
|---|---|---|
| Resistance per road km (lower = more parallel routes) | 0.051 | 0.063 |
| Min cut between hub areas (higher = harder to cut off) | 21 | 9 |
| Corridor conductance φ ×10⁻³ (higher = no narrow neck) | 0.84 | 0.91 |

- **Chennai is easier to cut apart.** Fewer than half as many road closures separate its hubs.
- **Chennai has fewer parallel routes.** Its resistance per road km is about 25% higher.
- **Corridor necks are similar** in the two cities.
- **Choke hubs**, where the bottleneck is in the hub's own catchment: Electronic City and
  Whitefield (Bengaluru), Sholinganallur (Chennai).

Figure: `figures/compare_hubs_dark.png`.

## How to run

```bash
pip install osmnx networkx numpy scipy matplotlib
python tools/fetch_namma.py bengaluru city      # and: chennai city  -> data/*.graphml
python tools/analyze_network.py data/bengaluru_city_bengaluru_drive.graphml 300   # top-betweenness CSV (defines the CBD); slow
python tools/degree_correlations.py data/bengaluru_city_bengaluru_drive.graphml
python tools/cluster_landuse.py data/bengaluru_city_bengaluru_drive.graphml
python tools/hub_bottlenecks.py data/bengaluru_city_bengaluru_drive.graphml
#   ...repeat the per-city steps for chennai
python tools/assortativity_compare.py   # Bengaluru vs Chennai
python tools/hub_compare.py
```

The `results/<city>_top_betweenness.csv` files are included, so the slow `analyze_network.py`
step can be skipped.

`tools/fetch_namma.py` and `tools/analyze_network.py` are copies of Kanak Varma's team
pipeline scripts, included only so this branch runs on its own.
