#!/usr/bin/env python3
"""
Download a real road network and save it for the game (needs internet + osmnx).

  pip install osmnx
  python fetch_city.py "Adyar, Chennai, India" --name chennai_adyar
  python fetch_city.py "Chennai, India" --name chennai --subgraph 300
  python fetch_city.py "Bangalore, India" --point 12.9716 77.5946 --dist 1500 --name bangalore_centre

Writes maps/<name>_edges.csv (u,v) and maps/<name>_nodes.csv (id,lon,lat); the game then
lists it on the menu and draws it with real coordinates.
"""
import argparse
import csv
import os
import random

import networkx as nx

import sy_sim as sy

ap = argparse.ArgumentParser()
ap.add_argument("place", nargs="?", help='place name, e.g. "Adyar, Chennai, India"')
ap.add_argument("--name", required=True)
ap.add_argument("--point", type=float, nargs=2, metavar=("LAT", "LON"))
ap.add_argument("--dist", type=int, default=1500, help="metres around --point")
ap.add_argument("--subgraph", type=int, default=0, help="keep a connected ~N-node piece")
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

import osmnx as ox

if a.point:
    D = ox.graph_from_point(tuple(a.point), dist=a.dist, network_type="drive")
elif a.place:
    D = ox.graph_from_place(a.place, network_type="drive")
else:
    ap.error("give a place name or --point")
U = nx.Graph()
for u, v in nx.Graph(D.to_undirected()).edges():
    if u != v:
        U.add_edge(u, v)
U = sy.largest_component(U)
if a.subgraph and U.number_of_nodes() > a.subgraph:
    U = sy.ball_subgraph(U, a.subgraph, random.Random(a.seed))
new = {old: i for i, old in enumerate(U.nodes())}
os.makedirs("maps", exist_ok=True)
with open(f"maps/{a.name}_edges.csv", "w", newline="") as f:
    csv.writer(f).writerows((new[u], new[v]) for u, v in U.edges())
with open(f"maps/{a.name}_nodes.csv", "w", newline="") as f:
    csv.writer(f).writerows((new[o], D.nodes[o]["x"], D.nodes[o]["y"]) for o in U.nodes())
print(f"saved maps/{a.name}_edges.csv and _nodes.csv: {U.number_of_nodes()} nodes, {U.number_of_edges()} edges")
