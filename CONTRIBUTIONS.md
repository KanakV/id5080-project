# Contributions

## Who did what

| Area | Person |
|---|---|
| Bengaluru analysis, pipeline scripts, write-up (ANALYSIS.md) | Kanak Varma |
| Chennai analysis | *(add yourself here)* |
| Growth-over-time modeling | *(open)* |

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
