# Upgrade: play on complete city networks

This patch changes the app to play on the **entire converted connected road network** for Chennai and Bengaluru. It removes the Area and Size controls, increases the detective limit, sets sensible city-map defaults, and adds performance safeguards for large graphs.

## Install

1. Stop the running game in PowerShell with `Ctrl+C`.
2. Extract this patch ZIP directly into the existing inner `scotland_yard_roads` project folder (the folder containing `app.py`).
3. Allow these files to be overwritten: `app.py`, `sy_sim.py`, `README.md`, and `static/index.html`.
4. **Do not delete or replace your `maps` folder.** The patch ZIP intentionally contains no map CSVs.
5. Start the app again: `python app.py`.
6. Choose `Chennai (whole city)` or `Bengaluru (whole city)`. Both should display the complete converted connected graph.

## Defaults and performance

- City maps default to 30 detectives and 80 rounds; both controls remain editable (up to 100 detectives and 500 rounds).
- The app starts city batches at 25 games to avoid unexpectedly long waits. Increase the batch size after checking runtime.
- Large-map node/edge drawing is cached, large belief overlays are omitted from the animation, city batch multiprocessing is limited to protect RAM, and expensive betweenness estimates use sampling. The full graph topology is still used for gameplay.
- The strategic policy is made scalable for large maps: Mr. X maximizes distance from the nearest detective, while detectives pursue the last reported location after a reveal and favor high-degree junctions before the first reveal. This is the same policy on each map, but differs from the previous prototype's belief-average heuristic.
- On large maps, graph path-length and betweenness metrics are sampled estimates.
