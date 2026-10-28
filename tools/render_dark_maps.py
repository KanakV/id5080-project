"""
Re-render the main network maps on a dark theme (user feedback: light-gray-on-white
maps were hard to tell from the background).

Produces (suffix _dark):
  <tag>_betweenness_edges_dark.png — edges colored by mean endpoint betweenness (cached)
  <tag>_cycles_dark.png            — bridges vs cycle edges (bridges recomputed, cheap)
  <tag>_eccentricity_dark.png      — nodes colored by eccentricity (cached, if present)

Dark style: surface #1a1a19, base marks #3a3a37, sequential ramps dim -> bright.

Usage:  python render_dark_maps.py [path/to/graph.graphml]
"""
import gc
import os
import sys
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
from matplotlib.lines import Line2D

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
tag = os.path.splitext(os.path.basename(GRAPH))[0]
os.makedirs("figures", exist_ok=True)

SURFACE, BASE, INKD = "#1a1a19", "#3a3a37", "#f0efec"
SUBD = "#c3c2b7"
BLUE_RAMP = mpl.colors.LinearSegmentedColormap.from_list(
    "dkblue", ["#232b36", "#3987e5", "#bcd9ff"])
ECC_RAMP = mpl.colors.LinearSegmentedColormap.from_list(
    "dkecc", ["#2e5f98", "#3f8ef0", "#d6ecff"])
ORANGE = "#ff8a3d"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def style(fig, ax, title):
    # draw the title INSIDE the axes: osmnx axes fill the canvas, so a normal
    # title sits on the bbox_inches extension strip, which saves as transparent
    ax.text(0.5, 0.995, title, transform=ax.transAxes, ha="center", va="top",
            color=INKD, fontsize=12)
    fig.patch.set_facecolor(SURFACE)


log(f"loading {GRAPH}")
G = ox.load_graphml(GRAPH)
Gs = ox.truncate.largest_component(G, strongly=True)

# ---------------------------------------------------------------- betweenness (cached)
bc_file = f"results/{tag}_node_bc.npz"
if os.path.exists(bc_file):
    log("map 1: betweenness, edge-painted, dark ...")
    z = np.load(bc_file, allow_pickle=True)
    bc = dict(zip(z["nodes"].tolist(), z["bc"]))
    ev = np.array([(bc.get(u, 0.0) + bc.get(v, 0.0)) / 2 for u, v, _ in Gs.edges(keys=True)])
    s = np.sqrt(ev / max(ev.max(), 1e-12))
    ec = [BLUE_RAMP(x) for x in s]
    ew = list(0.3 + 2.4 * s)
    fig, ax = ox.plot_graph(Gs, bgcolor=SURFACE, node_size=0, edge_color=ec,
                            edge_linewidth=ew, show=False, close=False, figsize=(12, 12))
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, ev.max()), cmap=BLUE_RAMP)
    cb = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.01)
    cb.set_label("betweenness centrality (mean of edge endpoints)", color=SUBD)
    cb.ax.tick_params(colors=SUBD, labelsize=7)
    cb.outline.set_edgecolor(BASE)
    style(fig, ax, f"Betweenness centrality, edge-painted — {tag}")
    fig.savefig(f"figures/{tag}_betweenness_edges_dark.png", dpi=220,
                bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
else:
    log(f"skipping betweenness map: {bc_file} missing")

# ---------------------------------------------------------------- cycles / bridges
log("map 2: cyclic vs tree-like, dark (recomputing bridges) ...")
U = ox.convert.to_undirected(G)
del G, Gs; gc.collect()
Ssimple = nx.Graph(U)
Ssimple.remove_edges_from(nx.selfloop_edges(Ssimple))
multi_uv = {frozenset((u, v)) for u, v, k in U.edges(keys=True) if k > 0}
bridges = {frozenset(e) for e in nx.bridges(Ssimple)} - multi_uv
ec, ew = [], []
for u, v, k in U.edges(keys=True):
    if frozenset((u, v)) in bridges:
        ec.append(ORANGE); ew.append(0.9)
    else:
        ec.append("#6e6e68"); ew.append(0.35)
fig, ax = ox.plot_graph(U, bgcolor=SURFACE, node_size=0, edge_color=ec,
                        edge_linewidth=ew, show=False, close=False, figsize=(12, 12))
ax.legend(handles=[Line2D([], [], color="#6e6e68", lw=2, label="on a cycle (loops/grid)"),
                   Line2D([], [], color=ORANGE, lw=2, label="bridge / tree-like (no cycle)")],
          loc="lower left", fontsize=9, frameon=False, labelcolor=INKD)
style(fig, ax, f"Cyclic vs tree-like street segments — {tag}")
fig.savefig(f"figures/{tag}_cycles_dark.png", dpi=220, bbox_inches="tight", facecolor=SURFACE)
plt.close(fig)

# ---------------------------------------------------------------- eccentricity (cached)
ecc_file = f"results/{tag}_eccentricity.npz"
if os.path.exists(ecc_file):
    log("map 3: eccentricity, dark ...")
    z = np.load(ecc_file, allow_pickle=True)
    enodes, ecc = z["nodes"].tolist(), z["ecc_lb"]
    pos = {u: (U.nodes[u]["x"], U.nodes[u]["y"]) for u in enodes}
    x = np.array([pos[u][0] for u in enodes]); y = np.array([pos[u][1] for u in enodes])
    from matplotlib.collections import LineCollection
    segs = np.array([[[U.nodes[u]["x"], U.nodes[u]["y"]],
                      [U.nodes[v]["x"], U.nodes[v]["y"]]] for u, v in Ssimple.edges])
    fig, ax = plt.subplots(figsize=(11, 11), dpi=220)
    ax.set_facecolor(SURFACE)
    ax.add_collection(LineCollection(segs, colors=BASE, linewidths=0.25, zorder=1))
    sc = ax.scatter(x, y, c=ecc, cmap=ECC_RAMP, s=0.8, lw=0, zorder=2)
    ax.set_aspect(1 / np.cos(np.deg2rad(y.mean())))
    ax.autoscale(); ax.set_axis_off()
    cb = fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.01)
    cb.set_label("eccentricity (hops) — max BFS depth from the node", color=SUBD)
    cb.ax.tick_params(colors=SUBD, labelsize=8)
    cb.outline.set_edgecolor(BASE)
    style(fig, ax, f"Node eccentricity (graph depth) — {tag} — diameter ≥ {int(ecc.max())} hops")
    fig.savefig(f"figures/{tag}_eccentricity_dark.png", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
else:
    log(f"skipping eccentricity map: {ecc_file} missing")

log("done (dark maps)")
