# Contributions

## Who did what

| Area | Person |
|---|---|
| Bengaluru analysis, pipeline scripts, write-up (ANALYSIS.md) | Kanak Varma |
| Chennai and Bangalore, Pagerank Centrality, Detour Index comparison (WIP)| S.Sathyanand |
| Growth-over-time modeling | *(open)* |
| Scotland Yard Pursuit and Evasion Simulation and Analysis | Jayvardhan More |
| Degree correlations / assortativity (Bengaluru, Chennai, Barcelona), edit distance to optimal and model networks, hub-to-hub bottlenecks | Aakash (me24b072) |

Add a row (or extend one) in the same PR/commit as your work.

## Adding your city

Everything is keyed by the graphml filename (its "tag"). Run:

```bash
.venv/bin/python run_all.py chennai
```

That fetches `data/chennai_city_chennai_drive.graphml` once (internet needed;
Overpass responses cache in `cache/`), then runs every analysis step. All
outputs land in `results/` and `figures/` prefixed with the tag, so two cities
never overwrite each other. Commit `results/` + `figures/` for your city;
`data/` and `cache/` stay out of git (see `.gitignore`).

Write your findings into `ANALYSIS.md` (new section or a sibling
`ANALYSIS_chennai.md` — your call), embedding figures by relative path.

## Adding a new metric / analysis step

Conventions every script in this repo follows — keep them so the runner and
your teammates' assumptions keep working:

1. **Interface**: take the graphml path as `sys.argv[1]`, derive
   `tag = os.path.splitext(os.path.basename(GRAPH))[0]`, and prefix every output
   file with `{tag}_`.
2. **Outputs**: numbers → `results/{tag}_<name>_report.txt` (+ CSV if tabular);
   expensive arrays → `results/{tag}_<name>.npz` so re-renders don't recompute;
   figures → `figures/{tag}_<name>.png`.
3. **Graph handling**: analyses run on the largest component; strip self-loops
   before `nx.bridges` / `nx.k_core`-style algorithms; use `weight="length"`
   for anything path-based.
4. **Plot style**: maps on the dark surface `#1a1a19` (see `render_dark_maps.py`
   for the ramps); charts use the validated categorical palette — city
   `#3987e5`, ER `#d95926`, WS `#199e70` on dark (light variants in
   `model_comparison.py`). Keep series colors/order fixed; always include a
   legend; put titles inside the axes on osmnx maps (see the comment in
   `render_dark_maps.py` for why).
5. **Reproducibility**: seed randomness (`SEED = 42` everywhere here), log
   progress with timestamps, and make re-runs cheap (check for cached npz
   before recomputing).
6. **Register it**: add the script to the `STEPS` list in `run_all.py` (in
   dependency order) and one row to the pipeline table in `README.md`.

## Workflow

- Work on a branch per feature/city (`chennai-core`, `growth-model`, ...),
  merge to `main` when the step runs end-to-end.
- Don't commit `data/`, `cache/`, venvs, or share bundles — they're
  regenerable and large.
- Anything published (report, slides) must credit
  "© OpenStreetMap contributors" (ODbL) — see README.

## Aakash (me24b072) — details

Outputs are per city, prefixed with the graph tag (`<tag>` = e.g. `bengaluru_city_bengaluru_drive`).

- Degree correlations (Newman, Networks sec. 7.13): assortativity r (eq. 7.82, cross-checked vs networkx), joint degree distribution e_jk, neighbour degree k_nn(k), and local assortativity (Piraveenan et al. 2010; node contributions sum exactly to r) mapped with the top-betweenness intersections labelled. `degree_correlations.py` -> `figures/<tag>_degree_correlations_dark.png`, `figures/<tag>_local_assortativity_{map,cbd}_dark.png`, `results/<tag>_degree_correlations_report.txt`, cache `results/<tag>_assortativity.npz`.
- Local assortativity by road type (OSM highway class, flyovers, ring roads, NH/SH, distance from CBD, bootstrap CIs) -> `figures/<tag>_assortativity_by_road_type_dark.png`; extremes (most disassortative nodes/edges, centre-vs-boundary profile, high-assortativity clusters beyond 10 km) -> `figures/<tag>_assortativity_extremes_dark.png`. `cluster_landuse.py` checks OSM land use at those clusters (Bengaluru: industrial estates and planned BDA layouts, not IT corridors) -> `results/<tag>_cluster_landuse.txt`.
- City comparison (Bengaluru, Chennai, Barcelona): `assortativity_compare.py` -> `figures/compare_assortativity_dark.png`, `results/compare_assortativity_report.txt`.
- Edit distance to target networks (Barcelona as the reference optimum; Erdős–Rényi; scale-free Chung–Lu; ring–spoke): exact degree-sequence lower bound on the full city + greedy constructive upper bound on a Barcelona-sized central window, counting edges added and removed, with robustness (Schneider R, random / degree attack), bridges, n·λ₂, min cut / √n, efficiency and detour tracked against the number of edits. `edit_distance.py` -> `results/<tag>_edit_distance_report.txt`, `figures/<tag>_edit_distance_dark.png`, `figures/<tag>_edit_map_dark.png`; comparison `edit_compare.py` -> `figures/compare_edit_{metrics,distance}_dark.png`.
- Hub-to-hub bottlenecks for ~8 named hubs per city: effective resistance (unit and length-weighted, exact via sparse LU), corridor Fiedler λ₂ and sweep-cut conductance, max-flow min cut between hub areas and where that cut sits. `hub_bottlenecks.py` -> `results/<tag>_hub_bottlenecks.{csv,_report.txt}`, `figures/<tag>_hub_bottlenecks_dark.png`; comparison `hub_compare.py` -> `figures/compare_hubs_dark.png`.
- Shared helpers in `graph_utils.py` (copies of λ₂ / sweep cut / percolation routines so base scripts stay untouched). Barcelona added to `fetch_namma.py`; steps `assort`, `hubs`, `edits` added to `run_all.py`.
