# Things tried and dropped

Attempts that were replaced or abandoned, kept out of the repository so it holds only
final results. Each line says why it was dropped and where to find it, if anywhere.

| Attempt | Why it was dropped | Where it is |
|---|---|---|
| **Edit distance to target networks** (Barcelona, Erdős–Rényi, scale-free, ring–spoke): degree-sequence lower bound + greedy upper bound on 6 metrics | The target was a structure, so a city already shaped like it looked cheap. Random-graph targets need links several km long that can't be built. The two bounds measured different things. | Commit `c46b744` on `main` (`tools/edit_distance.py`, `tools/edit_compare.py`) |
| **Edit distance with flyovers** (local road = 1, flyover up to 2 km = 10), targets Barcelona / ring–spoke | The flyover cost of 10 was arbitrary. Roads were chosen to shorten detours, which ignored the target. | Not committed |
| **Edit distance with a betweenness-attack target** | A perfect grid survives a static betweenness attack almost intact (R = 0.49 of 0.5), so the target can't be reached, and the run spent its budget on long links | Not committed |
| **Global assortativity in a 10 km disc** (Bengaluru, Chennai, Hyderabad; Seoul as the benchmark) | The Seoul download from Overpass failed repeatedly, so there was no benchmark. Three cities gave no growth-speed pattern. | Not committed |
| **Hub comparison including Barcelona** | Barcelona is about 10 km across, so its hub pairs only overlap the Indian cities at 4–9 km; the comparison now covers Bengaluru vs Chennai over 4–30 km | Commits `c46b744` / `1ca9ea9` on `main` |
| **Hub plot of "effective parallel routes"** (N_eff = road distance / R_len) | Replaced by its inverse, resistance per road km (lower = better), which reads directly as "1 = a single road" | Not committed |
| **Hub plots against straight-line distance**, with an a + b·ln d fit | The apparent "5–7× steeper" growth came from comparing different distance ranges; plots now use road distance over a matched range | Not committed |
