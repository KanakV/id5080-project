# Scotland Yard on Networks (road-only)

One Mr. X vs. several detectives on any road network. Every edge is a road, every player moves one
edge per turn, no tickets or transport types. Yellow junctions, red roads.

## Setup and run
    pip install numpy scipy networkx flask
    python app.py                 # opens http://127.0.0.1:5000

**Main menu:** pick a map, set the rules, then
- **Watch one game:** animated, with fog of war. Mr. X only shows on reveal rounds. On smaller maps, the
  cyan glow shows where detectives think he could be; that overlay is hidden for full-city maps to keep the
  animation responsive. On the original board, the game runs on the real board picture (toggle it and change
  its brightness in the side panel).
- **Run simulations and analyse:** 25 to 5,000 games, then win rates with a 95% interval,
  game-length histogram, survival curve, heatmaps on the map, correlations between node degree /
  betweenness / closeness and outcomes, deadliest nodes, and a `games.csv` download.

Play styles (both sides): **Strategic** (Mr. X moves away from the nearest detective; after a sighting,
detectives pursue the last reported location, and before the first sighting they favor high-degree junctions)
or **Random**. The strategic policy is deliberately simple and scales to city-sized graphs.

## Maps
- **Original board:** the 199 stations at their true positions on the board picture, taxi edges only.
  Positions were read off the picture (`data/board_nodes.csv`), and `data/board_edges.csv` uses the
  real board numbers. See `data/LABEL_FIX.txt`: the original repo's `taxi_map.txt` had five shifted
  station numbers; with them corrected, no road crosses another when drawn on the board.
- **Erdős–Rényi, Watts–Strogatz, Barabási–Albert:** random graphs (node slider, "new random graph").
- **Chennai and Bengaluru (real road networks):** the game uses the entire largest connected component
  of each city road graph, drawn with intersections at their true map coordinates. There is no area selector
  or subgraph sampling: Bengaluru has about 155,000 intersections and Chennai about 69,000 in the supplied
  graphs. City maps default to 30 detectives and 80 rounds; adjust those values based on your experiment.

### Getting Chennai and Bengaluru into the game
The project zip did not contain the networks (the `data/*.graphml` files are git-ignored), so create
them once on your machine:

    # in your id5080-project folder (skip any file you already have)
    python fetch_namma.py chennai city
    python fetch_namma.py bengaluru city
    # then, in this folder:
    python prepare_city.py path/to/id5080-project/data
    python app.py

`prepare_city.py` converts every `*.graphml` in that folder to `maps/<name>_edges.csv` and
`maps/<name>_nodes.csv` (lon, lat). It streams the GraphML parser, so the 145 MB Bengaluru graph is fine.
Any converted map is played in full; if you also convert neighborhood GraphML files, they appear as separate
smaller maps in the selector.

## Rules
- Round t: Mr. X moves, then detectives move one at a time. Nobody may enter a node held by a
  detective; a detective with no free neighbour stays put.
- Landing on Mr. X, or boxing him in, wins for the detectives. Surviving all rounds wins for Mr. X.
- Reveal rounds scale from the board game's 3, 8, 13, 18, 24 when you change the round limit.

## Command line (no browser)
    python sy_sim.py --graph taxi --games 5000
    python sy_sim.py --graph er --n 199 --games 1000 --instances 20 --out results
See `python sy_sim.py --help`.

## Files
`app.py` (server), `static/index.html` (menu, animation, charts), `static/boards/scotland_yard.jpg`
(board picture), `sy_sim.py` (engine and CLI), `prepare_city.py`, `fetch_city.py`, `data/`, `maps/`.
The board artwork belongs to Ravensburger; it is included for your course project only.

## Caveats
- Whole-city networks have far more nodes than the original board, so their results should not be compared
  directly with the 199-node board without accounting for graph scale and the chosen detective/round settings.
- ER with mean degree 3.31 is often disconnected, so its giant component (a bit under 199 nodes) is used.
- Whole-city maps are large. Start with 25–100 games and increase the batch after checking runtime.
  The app limits multiprocessing on city-sized graphs to avoid excessive RAM use, caches the rendered road
  layer, avoids sending the enormous belief set to the animation, and samples expensive centrality estimates.
- City graph diameter/path-length and betweenness metrics are sampled estimates on large graphs; report them
  as estimates rather than exact values. The full road topology is still used for gameplay.
- City support has not been validated against every regional GraphML variant; if your maps render incorrectly,
  check that they were converted from OSMnx GraphML and include valid coordinates.
