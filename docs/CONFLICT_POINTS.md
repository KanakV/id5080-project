# Coordination friendliness: conflict points per intersection

Script: `tools/conflict_points.py` -> `figures/conflict_points.png`, `results/conflict_points_report.txt`

## Idea
Every point where two vehicle streams must merge, split or cross is a place where drivers
(or signals) have to coordinate. Fewer such **conflict points** at an intersection means
less to coordinate. This is the standard way traffic engineers compare intersection
layouts (e.g. a 4-leg junction has 32 conflict points, a roundabout 8).

## Definition (one number per intersection)
Intersection = graph node with >= 3 distinct neighbouring nodes ("legs").
Each leg has one entry point and one exit point (left-hand traffic, as in India).
A movement is entry(a) -> exit(b), a != b (no U-turns), allowed only if the directed
graph has both edges, so one-way streets are respected.

| Type | Count |
|---|---|
| Diverging | a stream entering at a leg and splitting into m exits gives m-1 points |
| Merging | m streams joining one exit give m-1 points |
| Crossing | two movements whose paths intersect inside the junction (chord intersection on the circle of leg points). Right turns across opposing traffic fall here |

Conflict points = diverging + merging + crossing. Reported: mean over all intersections
of a city.

## Validation (asserted in the script)
For a fully two-way junction this reproduces the textbook counts:
- 3 legs (T): 3 diverge + 3 merge + 3 cross = **9**
- 4 legs: 8 diverge + 8 merge + 16 cross = **32**
Closed form for n two-way legs: n(n-2) + n(n-2) + n^2(n-1)(n-2)/6 (3 legs: 9, 4 legs: 32, 5 legs: 80, with 5-leg diverge/merge 15 each, cross 50).

## References (checked from memory, please verify page numbers before citing)
- FHWA, *Roundabouts: An Informational Guide* (FHWA-RD-00-067, 2000) - conflict-point
  diagrams: 32 for a conventional 4-leg intersection vs 8 for a single-lane roundabout;
  also T-junction 9.
- Garber & Hoel, *Traffic and Highway Engineering* - conflict points (crossing, merging,
  diverging) at at-grade intersections.
- Conflict points as a proxy for safety/complexity: Hauer, *Observational Before-After Studies in Road Safety* / traffic-conflict-technique literature (Hyden 1987).

## Assumptions / limitations
- Potential conflicts, not observed: no signals, turn bans (beyond one-way), traffic volume,
  or lane counts. Signals reduce simultaneous conflicts but are not in OSM drive graphs reliably.
- Legs are neighbouring graph nodes; bearing from node to neighbour sets the clockwise order.
  Parallel edges between the same pair collapse to one leg.
- Dual carriageways modelled as separate nodes (not consolidated) appear as several small
  intersections, which lowers the per-intersection mean in cities that have many of them.
- A per-intersection mean compares layout complexity, not total network load; a city
  with more intersections per km still has more points to coordinate overall.

## Result
See `results/conflict_points_report.txt`. Bengaluru ~11.8 vs Chennai ~11.0 conflict points
per intersection (median 9 in both, i.e. mostly T-junctions): Bengaluru's difference comes
from more crossing conflicts (4.7 vs 4.3) at bigger junctions.
