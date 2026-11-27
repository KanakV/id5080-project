# Robustness and structure of the Chennai and Bengaluru road networks

Author: Omkar. Course project in complex networks (ID5080).

## The question

Which road network, Chennai or Bengaluru, is more robust when junctions fail,
and how does each one compare with the textbook networks (Erdős–Rényi random
graphs and scale-free networks)? The other scripts explain *why*, by describing
how the networks are structured, laid out in space, and how they grew.

## What is in this branch

```
road_robustness.py          main experiment: random failure vs targeted attacks
road_path_Struct2.py        how far apart junctions are (path length, diameter, clustering)
road_spatial_analysis2.py   street layout: orientation, circuity, density, road classes
road_growth_over_time.py    historical snapshots: growth, preferential attachment, robustness over time
Plots/                      all figures and summary tables produced by the scripts
```

## Setup

```
pip install osmnx networkx numpy scipy matplotlib geopandas
```

- Python 3.10 or newer. Every script downloads the drivable network from
  OpenStreetMap through `osmnx`, so it needs an internet connection.
- The first download of a city can take a few minutes. After that, `osmnx`
  caches the responses.
- For a zoomable plot window instead of the small VS Code preview, run
  `pip install pyqt5` and uncomment the two `matplotlib.use("QtAgg")` lines near
  the top of the script.

Run each script on its own, for example `python road_robustness.py`. They do not
depend on each other, except that `road_growth_over_time.py` is the slow one.

## Conventions

- Chennai is orange and Bengaluru is green. ER is blue. Scale-free networks are
  pink (γ = 2.5) and purple (γ = 3.5).
- Each road network is turned into a simple undirected graph (junctions as
  nodes, road segments as links) and reduced to its largest connected piece.
- Each comparison network (ER, scale-free) has the same number of nodes and about
  the same number of links as the city it is compared with.

---

## 1. `road_robustness.py`

**Purpose and intuition.** Take junctions away and watch the network fall
apart. The size of the giant component (the part of the city that is still
connected) tells you how robust the network is. Three ways of choosing which
junctions to remove show *why* a network is fragile:

| Removal | Models | Removes first |
|---|---|---|
| Random failure | a flood or accident anywhere | random junctions |
| Betweenness attack | key crossings and arterials fail | junctions on the most shortest routes |

If a betweenness attack hurts much more than a degree attack, the weakness
comes from position (corridors), not from hubs.

**What it does.**
1. Downloads each city and builds the ER and scale-free comparison networks.
2. For each network, computes the giant component after removing a growing
   fraction *f* of junctions, for all three removal types. This uses an exact
   one-pass method (add nodes back in reverse order with union-find).
3. Measures three numbers per curve: `f_c` (where the giant component falls
   below 1%), `f_50` (where it falls below 50%), and `R` (the area under the
   curve; larger means more robust).
4. Prints a diagnostics table (share of nodes in the initial giant component,
   share of degree ≤ 1 nodes, ⟨k²⟩/⟨k⟩, the Molloy–Reed prediction of `f_c`) and
   a verdict that says "no clear winner" when the metrics disagree or R differs
   by less than 5%.

**Things to know.**
- Betweenness is *estimated* from 500 sampled source nodes (`BC_SAMPLES`), because
  exact betweenness on 70,000 to 150,000 nodes is very slow. The ranking is
  therefore approximate.
- Random failure is averaged over `N_RANDOM_RUNS` (10) runs. Raise it to smooth
  the curve. Targeted curves use one fixed ranking.
- `GAMMAS = [2.5, 3.5]`. The labels are built from these values.
- `f_c` (the 1% rule) reads the long thin tail of the curve and is sensitive to
  small details. Prefer `R` and `f_50` when comparing.
- The road curves are jagged because removing one corridor junction can cut off a
  whole region at once. This is a finding, not noise.
- The scale-free comparison graphs are built from a degree sequence with a fixed
  mean degree. Their shape at low mean degree (few degree-1 nodes) affects how
  they compare with ER, so read the diagnostics table before drawing conclusions.
- Optional switches: `COMPARE_ON_LCC` (reduce comparison graphs to their largest
  component), `ADAPTIVE` (recompute betweenness during the attack; much slower),
  `BC_STABILITY_CHECK` (check the sampled ranking is stable).

**Outputs.** `{city}_robustness.png` (three panels), `chennai_vs_bengaluru.png`,
and the printed tables.

---

## 2. `road_path_Struct2.py`

**Purpose and intuition.** Measure how far apart junctions are. In a random
network a link can join any two places, so everything is a few hops from
everything. On roads, links only join nearby places, so the number of hops between
two points follows the physical distance between them. Long paths with no
shortcuts mean long trips rely on a few corridors, which connects to the
robustness results.

**What it does.**
1. Starts a breadth-first search from 200 random junctions to every other
   junction, giving the average path length (hops and km) and a lower bound on
   the diameter.
2. Computes the clustering coefficient and repeats the path measurements on an ER
   network of the same size.
3. Draws two iterations on the map: a random source, every road coloured by its
   hop distance from that source, and the longest path in red.
4. Plots the distribution of path lengths for road vs ER.

**Things to know.**
- Averages use 200 sampled sources (`N_SOURCES`) with all targets, which is
  millions of paths, so the mean is reliable. The diameter is a lower bound.
- Hop counts depend on how finely the roads were mapped into junctions. Compare
  the two cities mainly in km, and compare road with ER in hops (same size).
- The table also prints a small-world index (sigma). It is inflated by the very
  small clustering of the ER graph, so do not use it to draw conclusions. The
  path-length comparison is the meaningful result.
- Map options: `EDGE_ONLY`, `LINE_W`, `DPI` and `ZOOM`.

**Outputs.** `{city}_hop_maps.png`, `path_length_distributions.png`, summary table.

---

## 3. `road_spatial_analysis2.py`

**Purpose and intuition.** Roads are physical objects. Their geometry decides how
many real alternative routes exist when a road is blocked. A dense, direct mesh
offers many alternatives. This script describes the layout:

- **Orientation entropy (H):** how mixed the street directions are. A perfect grid
  gives about 1.4 and fully mixed directions give 3.58.
- **Circuity:** route distance divided by straight-line distance (1.0 means a
  straight route).
- **Intersection and street density:** junctions and km of road per km².
- **Road-class mix:** the share of road length by type (residential, primary, ...).

**What it does.**
1. Computes the bearing and length of every road and the entropy of their
   distribution.
2. Samples random origin-destination pairs, finds the shortest route and the
   straight-line distance, and computes circuity for about 300 pairs.
3. Divides junction count and road length by the city's administrative area.
4. Draws, per city, the orientation map, a polar histogram, and two circuity
   iterations (random origin and destination about 5 km apart).

**Things to know.**
- Area comes from the administrative boundary, so it includes lakes and open land.
  The comparison between cities is fair, the absolute density is approximate.
- Circuity is estimated from a sample of pairs (`N_SOURCES`, `TARGETS_PER_SOURCE`),
  and pairs closer than 500 m are ignored. The Chennai-Bengaluru gap is small, so
  compare medians as well as means.
- One-way streets and turns are ignored (all roads are treated as two-way).
- Map options: `LINE_W`, `DPI`, `ZOOM`.

**Outputs.** `{city}_spatial.png`, `circuity_comparison.png`, summary table.

---

## 4. `road_growth_over_time.py`

**Purpose and intuition.** Real networks grow. The Barabási–Albert model says
growth by preferential attachment (new links go to already well-connected nodes)
creates hubs. Historical OpenStreetMap snapshots (2018, 2020, 2022) let the script
test that growth rule on real roads and check whether growth changes robustness.

**What it does.**
1. Downloads each snapshot from the Overpass API as it was on the given date.
2. Growth and structure per year: nodes, road length, degree distribution, mean
   degree, path length and clustering against ER.
3. Preferential-attachment test: for old nodes, how many new links they gained as a
   function of their old degree, with a fitted exponent α (α = 1 would mean linear
   preferential attachment, α = 0 no preference).
4. Runs the robustness attacks on every snapshot.
5. Maps the new roads, growth by distance from the centre, and how the most
   central roads change from year to year.

**Things to know.**
- **Part of the "growth" is better mapping, not construction.** The script prints a
  warning when edges grow by more than 40% and reports how much new road lies
  within 300 m of the old network (often mapping completion).
- Overpass history only goes back to September 2012.
- Cities are downloaded in small tiles (`TILE_DEG`) with retries (`MAX_RETRIES`)
  and waiting passes (`PASSES`, `ROUND_WAIT`), because big historical queries make
  the public server drop the connection. Each finished tile is saved in
  `snapshot_cache/tiles/`, so an interrupted run resumes where it stopped. Do not
  delete that folder between runs. If you change `TILE_DEG`, delete it first.
- The first run is slow (possibly an hour or more per city). Reruns are fast.
- Road degrees only reach about 4 to 6, so the α fit rests on few points. Report it
  cautiously.

**Outputs (in `growth_plots/`).** `_1_structure_over_time`,
`_2_preferential_attachment`, `_3_robustness_over_time`, `_4_new_roads`,
`_5_corridors` for each city, and `cities_over_time.png`.

---

## Plots folder

| File | Produced by | What it shows |
|---|---|---|
| `Chennai_robustness.png`, `Bengaluru_robustness.png` | `road_robustness.py` | giant component vs fraction removed: random, betweenness, degree; city vs ER and scale-free |
| `chennai_vs_bengaluru.png` | `road_robustness.py` | the two cities on the same axes |
| `path_length_distributions.png` | `road_path_Struct2.py` | path-length histograms, road vs ER |
| `Chennai_hop_maps.png`, `Bengaluru_hop_maps.png` | `road_path_Struct2.py` | hop distance from a random source, two iterations |
| `Chennai_spatial.png`, `Bengaluru_spatial.png` | `road_spatial_analysis2.py` | orientation map, polar histogram, two circuity examples |
| `circuity_comparison.png` | `road_spatial_analysis2.py` | circuity histogram of sampled trips |
| table screenshots | all scripts | printed summary tables |

For how to read each figure, see the separate plot guide.

## Key findings

Values below come from the runs in `Plots/`. Robustness values for the roads and
ER come from the first run with random failure and betweenness attack. Curve
positions for the degree attack are read from the plots by eye, so treat them as
approximate.

| | Chennai | Bengaluru | ER (same size) |
|---|---|---|---|
| Nodes / links | 69,389 / 91,678 | 155,843 / 205,130 | same |
| Mean degree ⟨k⟩ | 2.64 | 2.63 | 2.6 |
| Random failure, `f_c` | 0.32 | 0.34 | 0.62 |
| Betweenness attack, `f_c` | 0.20 | 0.19 | 0.35 |
| Betweenness attack, `R` | 0.038 | 0.044 | 0.197 |
| Average path (hops) | 107.0 | 101.7 | 11.2 / 12.1 |
| Average path (km) | 15.1 | 15.6 | n/a |
| Orientation entropy H | 3.16 | 3.28 | n/a |
| Median circuity | 1.24 | 1.20 | n/a |
| Junctions per km² | 121 | 162 | n/a |

1. **Both road networks are far less robust than ER**: roughly half the random-failure
   threshold, and a much lower attack `R`.
2. **Random failure: a tie** between the cities (`f_c` 0.32 vs 0.34, `R` 0.172 vs 0.169).
3. **Betweenness attack: Bengaluru is slightly ahead** (`R` 0.044 vs 0.038, holds longer
   at roughly f = 0.02 to 0.06). The metrics disagree on `f_c`, so call the gap narrow.
4. **Corridors, not hubs.** A betweenness attack breaks the roads at about
   f ≈ 0.07, close to scale-free γ = 2.5. A degree attack only breaks them at
   about f ≈ 0.2, close to ER. Importance comes from position in the city.
5. **Roads are not small-world.** Paths are about 9 times longer than ER's, and
   distance follows geography (the colour rings in the hop maps).
6. **Both layouts are organic** (H close to the fully mixed 3.58). Bengaluru has the
   denser mesh and slightly lower circuity, which fits its narrow edge under
   targeted attack. This is an association, not proof.

## Limitations

- Betweenness is sampled, so the exact attack order is approximate.
- Differences between the cities are small, so avoid declaring a decisive winner.
- Only junction removal is modelled. Road capacity, traffic and one-way streets are
  not included in these scripts.
- OpenStreetMap detail differs between cities and years, and this affects hop
  counts and any growth measure.
- Administrative areas used for densities include non-urban land.

## Troubleshooting

| Problem | Fix |
|---|---|
| `ResponseStatusCodeError` or `ConnectTimeout` while downloading | The public Overpass server is busy or unreachable. Rerun later. Finished tiles are saved. Check `https://overpass-api.de/api/status`. |
| `Overpass kept failing for this tile` | Rerun, make `TILE_DEG` smaller (and delete `snapshot_cache/tiles/`), or try a quieter time of day. |
| Boundary check skipped in the spatial script | Install `geopandas`. The rest of the script still runs. |
| Plots too small to read | Use `ZOOM`, raise `DPI`, or switch to the Qt backend (see Setup). |
