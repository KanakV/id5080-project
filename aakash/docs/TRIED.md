# Things tried and dropped

Attempts that were replaced or abandoned, kept out of the final files so the branch holds only
the two analyses in the report: degree assortativity and hub-to-hub bottlenecks. Each line
says why it was dropped and where to find it.

| Attempt | Why it was dropped | Where it is |
|---|---|---|
| **Barcelona as a reference city** (assortativity, hubs, edit distance) | Barcelona's Cerdà grid comes from a different planning history, continent and era. Any Indian city differs from it, so the comparison says nothing about Bengaluru vs Chennai. | Commit `db3c5bb` (repo history); commit `c46b744` on `main` |
| **Indian peer cities as the reference** (Pune and Hyderabad: Census 2011 urban agglomeration growth +34% / +33%, density 10.1k / 6.3k per km², vs Bengaluru +49% / 11.4k and Chennai +29% / 9.3k) | Chosen to separate planning from geography and size, but not finished: the city downloads from Overpass were slow and unreliable. Kept as a stated limitation of the direct comparison (METHODS.md §4). | Not committed |
| **Edit distance to a robust network**: new local roads (Delaunay neighbours ≤ 300 m) needed to reach 60–80% of a square grid's robustness under random failure and degree attack | Removed from the final scope of the report. Result at removal: Bengaluru 239 vs Chennai 349 new roads per 1,000 existing roads for 80%. | Commit `db3c5bb` (repo history) (`tools/edit_spatial.py`, `docs/EDIT_DISTANCE.md`) |
| Earlier edit-distance versions: target networks (Barcelona, Erdős–Rényi, scale-free, ring–spoke); flyovers at cost 10; a betweenness-attack target | Targets that were structures rewarded a city for already resembling them. The flyover cost was arbitrary. The betweenness target cannot be reached (a grid scores R = 0.49 of 0.5). | Commit `c46b744` on `main`; the others were not committed |
| **Conflict points per intersection** (diverging, merging and crossing points; Bengaluru 11.8 vs Chennai 11.0 per intersection) | Removed from the final scope of the report | Commit `db3c5bb` (repo history) (`tools/conflict_points.py`, `docs/CONFLICT_POINTS.md`) |
| **Global assortativity in a 10 km disc** (Bengaluru, Chennai, Hyderabad; Seoul as the benchmark) | The Seoul download failed repeatedly, so there was no benchmark. Three cities gave no pattern by growth speed. | Not committed |
| **Hub plot of "effective parallel routes"** (N_eff = road distance / R_len) | Replaced by its inverse, resistance per road km (lower = better), which reads directly as "1 = a single road" | Not committed |
| **Hub plots against straight-line distance**, with an a + b·ln d fit | The apparent "5–7× steeper" growth came from comparing different distance ranges. Plots now use road distance over a matched range. | Not committed |
