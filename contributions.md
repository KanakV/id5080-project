# Contributions

## Aakash (me24b072)

- Degree correlations (Newman, Networks sec. 7.13): assortativity r (eq. 7.82, cross-checked vs networkx), joint degree distribution e_jk, neighbour degree k_nn(k), and a map of local assortativity with the top-betweenness intersections labelled. `degree_correlations.py` -> `figures/bengaluru_degree_correlations_dark.png`, `results/bengaluru_degree_correlations_report.txt`.
- Local assortativity by road type: nodes classed by OSM highway tag (expressway/trunk, primary, secondary, tertiary, local), flyovers (bridge/viaduct), ring roads, NH/SH refs, and distance from the CBD; full-city and CBD-zoom maps of local assortativity plus a per-category comparison with bootstrap CIs (`figures/bengaluru_local_assortativity_*_dark.png`, `figures/bengaluru_assortativity_by_road_type_dark.png`).
- Extremes of local assortativity: most disassortative nodes/edges (hub junctions with dead-end spokes), radial centre-vs-boundary profile, and the high-assortativity clusters beyond 10 km from the CBD (planned grid layouts), in `figures/bengaluru_assortativity_extremes_dark.png` and the report.
