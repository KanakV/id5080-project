# Scotland Yard on Networks (road-only)

One Mr. X vs. several detectives on any road network. Every edge is a road, every player moves one
edge per turn, no tickets or transport types. Yellow junctions, red roads.

## Setup and run
    pip install numpy scipy networkx flask
    python app.py                 # opens http://127.0.0.1:5000

**Main menu:** pick a map, set the rules, then
- **Watch one game:** animated, with fog of war. Mr. X only shows on reveal rounds; the cyan glow is
  where the detectives think he could be. On the original board the game runs on the real board picture
  (toggle it and change its brightness in the side panel).
- **Run simulations and analyse:** 100 to 10,000 games in parallel, then win rates with a 95% interval,
  game-length histogram, survival curve, heatmaps on the map, correlations between node degree /
  betweenness / closeness and outcomes, deadliest nodes, and a `games.csv` download.

Play styles (both sides): **Strategic** (Mr. X runs away from the nearest detective; detectives track the
set of nodes where he could be) or **Random**.

## Maps
- **Original board:** the 199 stations at their true positions on the board picture, taxi edges only.
  Positions were read off the picture (`data/board_nodes.csv`), and `data/board_edges.csv` uses the
  real board numbers. See `data/LABEL_FIX.txt`: the original repo's `taxi_map.txt` had five shifted
  station numbers; with them corrected, no road crosses another when drawn on the board.
- **Erdős–Rényi, Watts–Strogatz, Barabási–Albert:** random graphs (node slider, "new random graph").
- **Chennai and Bengaluru (real road networks):** drawn with intersections at their true map
  coordinates, no satellite image. Whole cities are far too big for one game (Bengaluru has about 155,000
  intersections), so the game uses a piece of 30 to 3,000 nodes: pick an area (T. Nagar, MG Road, ...) or
  a random one.

### Getting Chennai and Bengaluru into the game
The project zip did not contain the networks (the `data/*.graphml` files are git-ignored), so create
them once on your machine:

    # in your id5080-project folder (skip any file you already have)
    python fetch_namma.py chennai city
    python fetch_namma.py bengaluru city
    # then, in this folder:
    python prepare_city.py path/to/id5080-project/data
    python app.py

`prepare_city.py` converts every `*.graphml` in that folder (city, neighbourhood or metro scale) to
`maps/<name>_edges.csv` and `maps/<name>_nodes.csv` (lon, lat). It streams the file, so the 145 MB
Bengaluru graphml is fine. Neighbourhood graphs (T. Nagar, Mylapore, Jayanagar, Chickpet) smaller than the
size setting are played whole. Alternatively `pip install osmnx` and use `fetch_city.py` to download straight
into `maps/`.

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
- Cutting a piece out of a city adds artificial dead ends at its edge; compare pieces of the same size.
- ER with mean degree 3.31 is often disconnected, so its giant component (a bit under 199 nodes) is used.
- On whole-city scale Mr. X would nearly always win a 24-round game; scale rounds and detectives with the
  graph and report it.
- City support was tested on a synthetic 62,500-node graphml in osmnx's format, not on your real files.
