# Edit distance to a robust street network: Bengaluru vs Chennai

**Question.** How many new roads would each city centre need before its street network
becomes robust? This is measured against the same absolute target for both cities.

**Answer.** To reach 80% of a perfect grid's robustness under both random failure and
targeted attack:

| City | New local roads per 1,000 existing roads |
|---|---|
| Bengaluru | **239** |
| Chennai | **349**, about 1.5× more |

Bengaluru's centre needs fewer new roads at every target level. This holds both when roads
are chosen well and when they are chosen at random.

Script: `tools/edit_spatial.py`. Outputs: `results/edit_spatial_report.txt`,
`results/edit_spatial.csv` and `figures/edit_spatial_dark.png`.

---

## 1. What "edit distance" means here

The graph edit distance between two networks is the fewest insertions and deletions of nodes
and edges that turn one into the other [Sanfeliu & Fu 1983]. Computing it exactly is NP-hard
[Zeng et al. 2009], and impossible at about 9,000 junctions. We therefore use the restricted
version that answers the planning question:

> **Edit distance = the number of new road segments that must be added before the network
> reaches a target robustness.**

- **Only insertions.** Roads are added, never removed. Demolishing roads is not a realistic
  way to make a city more robust.
- **Every road costs 1.** No weights, so the number means exactly what it says.
- **The target is a level of a metric, not another network.** An earlier version used
  Barcelona and model networks as targets. That rewarded a city for already resembling the
  target's shape (§7).

## 2. The steps

| Step | What we do | Why |
|---|---|---|
| **Window** | The 9,000 junctions nearest the CBD (centroid of the top-betweenness junctions), largest connected piece | The same size for both cities, so counts are comparable. The whole cities differ 2.2× in size. |
| **Robustness metric** | Schneider R: remove junctions one by one, track the share of junctions still in the largest connected piece, average over the whole removal. Maximum 0.5. | Standard single-number robustness [Schneider et al. 2011]. It rewards networks that stay connected long, not just ones that survive the first removal. |
| **Failure modes** | (a) Random failure: random order, mean of 3 fixed-seed orders. (b) Degree attack: always remove the currently most-connected junction, recomputed after every removal. | The two classic tests [Albert et al. 2000]. The recomputed attack is the stronger, standard version [Holme et al. 2002]. Computed exactly with the Newman–Ziff union-find sweep [Newman & Ziff 2000]. |
| **Absolute target** | A % of the R of a perfect 95 × 95 square grid (same number of junctions): R = 0.327 under random failure, 0.258 under degree attack. Levels: 60 / 70 / 80%. | The square grid is the textbook planned street layout and the usual reference for street networks [Barthélemy 2011]. Both cities face the same finish line. |
| **Which roads may be built** | A new road joins two **neighbouring** junctions: Delaunay neighbours ≤ 300 m apart, both with fewer than 4 roads | Real roads join nearby junctions without cutting through blocks. Delaunay triangulation is a common model of "who can connect to whom" in planar street networks [Masucci et al. 2009; Barthélemy 2011]. The 4-road cap keeps junctions ordinary crossroads. |
| **Order of building (main rule)** | Greedy: build first the roads that keep the most junctions attached to the largest piece during the collapse (§3). Recompute R exactly after every batch. | Greedy edge addition is the standard approach to improving robustness [Beygelzimer et al. 2005]. |
| **Baseline rule** | The same candidate roads, in random order | Shows that the city comparison does not depend on how cleverly roads are chosen (§5) |
| **Edit distance** | New roads at the first checkpoint where **both** failure modes reach the level, per 1,000 existing roads | Checkpoints are every 5% of the buildable roads, plus finer steps at the start |

## 3. The greedy rule in one paragraph

R is the area under the "largest piece vs junctions removed" curve.

1. Take 10 snapshots of each collapse.
2. At each snapshot, see which pieces the network has broken into.
3. A new road between junctions *a* and *b* helps at that snapshot only if both are still
   standing and lie in different pieces. The road then joins them, and the largest piece
   grows (by the smaller piece's size, when one side is the largest piece).
4. Add this up over the snapshots to estimate how much the road raises R. Build the roads
   with the largest estimates first.

Within one batch, a junction gets at most one new road, and a stranded piece is reconnected
at most once. Once no road helps any more, the remaining buildable roads are added in random
order, so both rules can use every road.

## 4. Results

**Starting robustness** (R, and % of the grid):

| | Random failure | Degree attack |
|---|---|---|
| Perfect grid (ceiling) | 0.327 (100%) | 0.258 (100%) |
| **Bengaluru centre** | 0.214 (66%) | 0.125 (48%) |
| **Chennai centre** | 0.179 (55%) | 0.105 (41%) |
| Barcelona, whole city (reference only) | 0.251 (77%) | 0.167 (65%) |

**Edit distance**: new roads per 1,000 existing roads until **both** modes reach the level.

| Level | Bengaluru, greedy | Chennai, greedy | Chennai ÷ Bengaluru | Bengaluru, random | Chennai, random |
|---|---|---|---|---|---|
| 60% | **78** | **139** | 1.8× | 143 | 233 |
| 70% | **143** | **233** | 1.6× | 239 | 325 |
| 80% | **239** | **349** | 1.5× | 334 | 407 |

**With every buildable road built** (about 340 per 1,000 roads in Bengaluru, 425 in Chennai),
both cities reach about 99% under random failure and about 90% under degree attack. Both
centres *can* be made robust with local roads alone; Chennai simply needs more of them.

**Figure** (`figures/edit_spatial_dark.png`):

- **A:** the edit distance at each level.
- **B:** robustness, the weaker of the two modes, against roads added. The greedy rule is
  solid and the random baseline is dashed.
- **Maps:** where the new roads go.

## 5. Why the result can be trusted

1. **It does not depend on the greedy rule.** With roads added at random, Chennai still
   needs more at every level: 233 vs 143, 325 vs 239, 407 vs 334. The ranking is the same
   under both rules.
2. **The greedy rule does what it claims.** At the same number of roads, it raises robustness
   faster than random roads in both cities (Figure B: the solid lines are above the dashed
   lines). For example, it reaches 60% with 78 roads per 1,000 vs 143 for random in
   Bengaluru, and 139 vs 233 in Chennai.
3. **Robustness is computed exactly, never estimated.** The greedy rule's estimate only
   decides the *order* of building. Every reported R is a full recomputation.
4. **Deterministic.** Fixed seeds throughout. Rerunning gives identical numbers.
5. **Resolution.** Checkpoints are every 5% of the buildable roads, about ±20 roads per
   1,000. The gap between the cities (61–110 per 1,000 with the greedy rule, 73–90 with
   random roads) is several checkpoints wide.

## 6. What this does and does not show

**It shows:**
- Chennai's centre starts less robust than Bengaluru's under both failure modes. Its
  starting values (55% and 41%) sit below Bengaluru's (66% and 48%).
- Chennai's centre needs about 1.5–1.8× more new local roads to reach the same absolute
  robustness.
- Random failure is cheap to fix. Most of the early gain comes from closing dead ends into
  loops.
- Withstanding targeted attack on major junctions is what drives the cost.

**It does not show:**
- **The real cost of construction.** Every road counts 1, whatever its length, land price or
  the demolition it needs.
- **The minimum number of roads.** The greedy order is good, but not proven optimal. The
  true minimum could be lower for both cities.
- **Anything about the whole city.** The window covers the central 9,000 junctions,
  3.5 km around Bengaluru's CBD and 4.4 km around Chennai's.
- **Traffic.** Robustness here means staying connected, not travel time or congestion.

**Known simplifications:**
- **Edge of the window.** The window's edge creates artificial dead ends. Both cities are cut
  the same way, but this lowers both starting values.
- **Undirected network.** Roads are treated as two-way, and lanes and capacity are ignored,
  as in the rest of the project.
- **Static betweenness attack not used.** It is not a usable target: a perfect grid survives
  it almost completely (R = 0.49 of a possible 0.5), so no finish line can be set from it.

## 7. Earlier versions, and why they were replaced

These are not in the repository; see [TRIED.md](TRIED.md) for where to find them.

| Version | Problem | Replaced by |
|---|---|---|
| Edit distance to Barcelona, Erdős–Rényi, scale-free and ring–spoke networks (`edit_distance.py`) | The target was a structure, so a city whose layout already resembled it looked cheap. Random-graph targets need links several km long that can't be built. | An absolute robustness target |
| Distance-constrained version with flyovers (cost 10) | The flyover cost was arbitrary, and road choice ignored the target (it only shortened detours) | Local roads only, chosen by robustness gain |
| Robustness target including a static betweenness attack | That mode's grid ceiling cannot be reached, so the run spent its budget on long links | Random failure and degree attack only |

## 8. Reproduce

```
python tools/edit_spatial.py        # about 3 minutes; Bengaluru and Chennai by default
```

It needs `data/<city>_city_<city>_drive.graphml` and `data/barcelona_city_barcelona_drive.graphml`
(reference line only; download with `tools/fetch_namma.py <city> city`), and
`results/<tag>_assortativity.npz`, which holds the CBD (from `degree_correlations.py`).

## References

- Albert, R., Jeong, H. & Barabási, A.-L. (2000). Error and attack tolerance of complex networks. *Nature* 406, 378–382.
- Barthélemy, M. (2011). Spatial networks. *Physics Reports* 499, 1–101.
- Beygelzimer, A., Grinstein, G., Linsker, R. & Rish, I. (2005). Improving network robustness by edge modification. *Physica A* 357, 593–612.
- Holme, P., Kim, B. J., Yoon, C. N. & Han, S. K. (2002). Attack vulnerability of complex networks. *Physical Review E* 65, 056109.
- Masucci, A. P., Smith, D., Crooks, A. & Batty, M. (2009). Random planar graphs and the London street network. *European Physical Journal B* 71, 259–271.
- Newman, M. E. J. & Ziff, R. M. (2000). Efficient Monte Carlo algorithm and high-precision results for percolation. *Physical Review Letters* 85, 4104–4107.
- Sanfeliu, A. & Fu, K.-S. (1983). A distance measure between attributed relational graphs for pattern recognition. *IEEE Transactions on Systems, Man, and Cybernetics* 13(3), 353–362.
- Schneider, C. M., Moreira, A. A., Andrade, J. S., Havlin, S. & Herrmann, H. J. (2011). Mitigation of malicious attacks on networks. *PNAS* 108(10), 3838–3841.
- Zeng, Z., Tung, A. K. H., Wang, J., Feng, J. & Zhou, L. (2009). Comparing stars: on approximating graph edit distance. *Proceedings of the VLDB Endowment* 2(1), 25–36.
