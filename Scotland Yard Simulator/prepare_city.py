#!/usr/bin/env python3
"""
Turn the road networks your team already downloaded (osmnx .graphml files, e.g. from
fetch_namma.py) into the compact map files the game reads.  Needs only networkx.

  python prepare_city.py path/to/bengaluru_city_bengaluru_drive.graphml
  python prepare_city.py path/to/chennai_city_chennai_drive.graphml
  python prepare_city.py path/to/data            # converts every *.graphml in that folder

Writes into ./maps/ :
  <name>_edges.csv   one line per street segment: node,node  (undirected, no self-loops)
  <name>_nodes.csv   node,lon,lat                (so the game draws the real map layout)

<name> defaults to the file name without "_drive", e.g. bengaluru_city_bengaluru,
chennai_city_chennai, chennai_nbhd_tnagar.  Only the largest connected component is kept.
Parsing is streamed, so the 145 MB Bengaluru file does not need much memory (a minute or two).
"""
import argparse
import csv
import glob
import os
import sys
import xml.etree.ElementTree as ET

import networkx as nx


def parse_graphml(path):
    keys, xs, ys, edges = {}, {}, {}, []
    cur, graph_el, seen = None, None, 0
    for ev, el in ET.iterparse(path, events=("start", "end")):
        tag = el.tag.rsplit("}", 1)[-1]
        if ev == "start":
            if tag == "graph":
                graph_el = el
            elif tag == "node":
                cur = el.get("id")
            continue
        if tag == "key":
            keys[el.get("id")] = (el.get("for"), el.get("attr.name"))
        elif tag == "data" and cur is not None:
            kind = keys.get(el.get("key"), (None, None))
            if kind[0] == "node" and kind[1] in ("x", "y") and el.text:
                (xs if kind[1] == "x" else ys)[cur] = float(el.text)
        elif tag == "node":
            cur = None
            seen += 1
        elif tag == "edge":
            edges.append((el.get("source"), el.get("target")))
            seen += 1
        if tag in ("node", "edge") and seen % 50000 == 0 and graph_el is not None:
            graph_el.clear()                     # free memory, nothing else needs it
    return xs, ys, edges


def convert(path, name, out_dir):
    print(f"reading {path} ...", flush=True)
    xs, ys, edges = parse_graphml(path)
    if not xs:
        sys.exit("no node coordinates (x, y) found; is this an osmnx graphml file?")
    G = nx.Graph()
    for u, v in edges:
        if u != v and u in xs and v in xs and u in ys and v in ys:
            G.add_edge(u, v)
    comp = max(nx.connected_components(G), key=len)
    G = G.subgraph(comp).copy()
    ids = {old: i for i, old in enumerate(G.nodes())}
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, f"{name}_edges.csv"), "w", newline="") as f:
        csv.writer(f).writerows((ids[u], ids[v]) for u, v in G.edges())
    with open(os.path.join(out_dir, f"{name}_nodes.csv"), "w", newline="") as f:
        csv.writer(f).writerows((ids[o], f"{xs[o]:.6f}", f"{ys[o]:.6f}") for o in G.nodes())
    print(f"  saved {name}: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges "
          f"(largest component of {len(xs):,} nodes)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help=".graphml file, or a folder containing .graphml files")
    ap.add_argument("--name", help="output name (single file only)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "maps"))
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.path, "*.graphml"))) if os.path.isdir(a.path) else [a.path]
    if not files:
        sys.exit("no .graphml files found")
    for f in files:
        stem = os.path.splitext(os.path.basename(f))[0]
        name = a.name if (a.name and len(files) == 1) else (stem[:-6] if stem.endswith("_drive") else stem)
        convert(f, name, a.out)
    print("done. Restart the game (python app.py) and pick the map from the menu.")


if __name__ == "__main__":
    main()
