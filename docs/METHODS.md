# Methods: Aakash's analyses

How every number in Aakash's part of the project is computed. Each section gives the
definition, the procedure the code follows, the constants, and the script. Scripts live in
`tools/` and are run from the repo root (`python tools/<script>.py ...`).

---

## 1. The network

| Step | What is done | Where |
|---|---|---|
| Source | OpenStreetMap drive network inside the municipal boundary, via osmnx (`network_type="drive"`) | `fetch_namma.py` |
| Simplification | osmnx default: degree-2 shape points merged, so a node is a junction or a dead end, and an edge is a road segment between junctions | osmnx |
| Undirected, simple | One-way directions dropped; parallel edges merged keeping the **shortest** `length`; self-loops dropped | `graph_utils.load_simple` |
| Giant component | Only the largest connected component is kept | `graph_utils.load_simple` |
| Coordinates | lon/lat → km with an equirectangular projection around the city mean (1° = 111.2 km, x scaled by cos lat) | `graph_utils.to_km` |

Exception: `degree_correlations.py` uses the whole simple graph, not only the giant component.

**Weights.** Unless a section says otherwise, every edge has weight 1 (topology only).
Edge length is used only in the following places:

- R_len: conductance = 1 / length;
- road distance between hubs.

Lanes, speed and capacity are never used.

---

## 2. Degree assortativity (`degree_correlations.py`, `assortativity_compare.py`, `cluster_landuse.py`)

### 2a. Global r (Newman 2003, *Networks* eq. 7.82)

Over all M edges, with degrees k_i, k_j at the two ends:

```
r = [ Σ_e k_i k_j / M  −  ( Σ_e (k_i + k_j) / 2M )² ]
    / [ Σ_e (k_i² + k_j²) / 2M  −  ( Σ_e (k_i + k_j) / 2M )² ]
```

- r > 0: busy junctions link to busy junctions (assortative).
- r < 0: hubs link to dead ends (disassortative).

The code checks its own value against `networkx.degree_assortativity_coefficient`, to 1e-9.


### 2b. Local assortativity ρ_v (Piraveenan et al. 2010)

Uses excess degree q = k − 1. Its mean μ_q and variance σ_q² are taken over edge ends.

```
ρ_v = q_v · Σ_{w ~ v} (q_w − μ_q) / (2M σ_q²)        and   Σ_v ρ_v = r   (asserted to 1e-6)
```

- **Plotted as ρ·N,** so the mean over all nodes equals r.
- **Road-category means:** 95% bootstrap CI.
- **Boundary clusters:**
  - Take the nodes more than 10 km from the CBD whose ρ·N is at least the city's 90th percentile.
  - Link nodes within 150 m of each other (single linkage).
  - Keep clusters of 40 nodes or more.
- **Script:** `degree_correlations.py`. It also produces e_jk (the joint degree distribution
  over edge ends) and k_nn(k) (the mean neighbour degree), compared with the uncorrelated
  baseline ⟨k²⟩/⟨k⟩.

### 2c. CBD

The CBD is the centroid of the top-8 betweenness junctions, read from
`results/<tag>_top_betweenness.csv`. That file comes from the team script `analyze_network.py`
(length-weighted betweenness on the directed graph, 300 sampled sources, seed 42). Every
"distance from CBD" and every central window is measured from this point.

---

## 3. Hub-to-hub bottlenecks (`hub_bottlenecks.py`, `hub_compare.py`)

### 3a. Hubs

- **Choice:** eight named activity centres per city, from approximate centre coordinates
  (`graph_utils.HUBS`).
- **Snapping:** each is moved to the nearest junction.
- **Dropped hubs:** any more than 1.5 km from the network. Only Tambaram (9.4 km outside the
  Chennai boundary) is dropped.
- **Hub area:** every node within 1 km of the snapped hub.

### 3b. Effective resistance (resistance distance)

Treat the road network as an electrical circuit: each segment is a resistor, and a unit
current is pushed in at hub a and out at hub b. The voltage needed is R_ab.

```
R_ab = (e_a − e_b)ᵀ L⁺ (e_a − e_b) = X_aa + X_bb − 2 X_ab
```

- **Computing L⁺:** never formed directly. Ground one non-hub node: delete its row and column
  from the Laplacian, giving L_g, which is invertible.
  1. Factorise L_g once with sparse LU.
  2. Solve for the 8 hub columns.
  3. X is the hub × hub block of L_g⁻¹.
- **Exactness:** this is exact, and the choice of ground node does not change R_ab.

The two weightings:

| Version | Conductance of a segment | Units of R | Meaning |
|---|---|---|---|
| R_unit | 1 | "segments" | Topology only |
| R_len | 1 / length_km (length floored at 1 m) | km | A road of length L is a resistor of resistance L |

**Normalised by road distance.** This is the quantity plotted in `figures/compare_hubs_dark.png`.

```
resistance per road km = R_len / d_road
```

- **d_road:** the shortest road distance between the two hubs (Dijkstra on length).
- **Reading it:**
  - A single road of length L gives R_len = L, so the ratio is 1.
  - k equally long parallel roads give R = L/k, so the ratio is **1/k**.
- **Example:** a value of 0.05 means the network carries current between the hubs like about 20
  parallel roads of the shortest-route length. Lower is better.
- **Why divide:** R_len grows with distance. Dividing by d_road takes that growth out, so pairs
  at different distances, and cities of different size, can be compared.

### 3c. Corridor conductance φ

The corridor is the set of nodes p inside the ellipse

```
|p − a| + |p − b| ≤ 1.3 d_ab + 2 km
```

- d_ab is the straight-line distance between the hubs.
- The corridor's giant component is taken, **unweighted**.
- Its Fiedler sweep cut (`graph_utils.sweep_cut`) gives the corridor's worst bottleneck:
  1. Sort nodes by the Fiedler vector of the normalised Laplacian, divided by √degree.
  2. For every prefix S, compute φ(S) = cut(S) / min(vol S, vol V∖S), where vol is the sum of degrees.
  3. Keep the lowest φ with at least 10% of the volume on each side, so the cut cannot just clip
     off one cul-de-sac.
- It reports:
  - λ₂;
  - φ, the conductance of the corridor's worst cut;
  - `splits`, which records whether that cut separates a from b.

φ is reported × 10³. Higher is better: there is no narrow neck between the hubs.

### 3d. Min cut between hub areas and the choke-hub rule

- **Min cut:** unit-capacity max-flow between the two 1 km hub areas, with nodes shared by both
  areas removed. The result is the fewest road segments whose loss separates the two areas.
  - It can never exceed the number of roads leaving the smaller hub's area. The script logs that number.
- **Choke-hub rule** (`graph_utils.classify_cuts`):
  - A hub is **choked at c** if its min cut to every other hub (allowing one exception, and needing
    at least 3 pairs) equals the same value c.
  - The bottleneck then sits in that hub's own catchment, not in the corridor.
  - These pairs are drawn hollow in panel B.

### 3e. Cross-city comparison (`hub_compare.py`)

- **Panels A–C:** each metric plotted against d_road, with lines at the median in each of
  7 geometric distance bins.
- **Matched range:** the road distances every compared city has, from max(min d_road) to
  min(max d_road). Bengaluru vs Chennai gives 4.2–30.4 km.
- **Panel D:** the ratio of medians inside that range, oriented so that higher = better.
  For resistance this is baseline ÷ city.
- **Map colours:** one shared scale, from the 2nd to the 98th percentile of all pairs.
  - **Caveat:** pair distance confounds the colours. Panel D is the fair comparison.

---

## 4. Edit distance to a robust street network (`edit_spatial.py`)

The full write-up is in [EDIT_DISTANCE.md](EDIT_DISTANCE.md).

- **Window:** the 9,000 junctions nearest the CBD.
- **Metric:** Schneider R under random failure and adaptive degree attack.
- **Target:** 60 / 70 / 80% of a 95×95 square grid's R, the same for both cities.
- **Edits:** new local roads only, between Delaunay-neighbour junctions ≤ 300 m apart, both
  ends with fewer than 4 roads. Each road counts 1.
- **Order:** greedy by robustness gain, with exact R after every batch, plus a random-order baseline.

---

## 5. Conflict points per intersection (`conflict_points.py`)

The full write-up is in [CONFLICT_POINTS.md](CONFLICT_POINTS.md).

---

## 6. Constants at a glance

| Constant | Value | Used in |
|---|---|---|
| Random seed | 42 | Everywhere |
| CBD | centroid of the top-8 betweenness junctions | §2c, §4 |
| Assortativity clusters | > 10 km from the CBD, ≥ 90th percentile ρ·N, 150 m link, ≥ 40 nodes | §2b |
| Hub snap limit / hub area radius | 1.5 km / 1 km | §3 |
| Corridor ellipse | 1.3 d + 2 km | §3c |
| Sweep-cut balance | ≥ 10% volume per side | §3c |
| Length floor in R_len | 1 m | §3b |
| Edit distance: window / new road / target | 9,000 junctions / Delaunay neighbour ≤ 300 m, both degree < 4 / 60–80% of a 95×95 grid's R | §4 |
