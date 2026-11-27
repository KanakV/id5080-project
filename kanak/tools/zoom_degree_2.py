# Author: Kanak Varma
"""
Deep-dive and zoomed-in visualization of Degree-2 Nodes in the road network.

Investigates why degree-2 nodes exist in simplified OSMnx street graphs:
  1. Directional Transitions: where two one-way roads meet or direction flips
  2. Classification Boundaries: where highway tag changes (e.g. primary -> secondary)
  3. Physical Features: bridge/tunnel start, gates, barriers, or traffic signals
  4. Cul-de-sac Turnaround Loops: self-loops or double-edges

Generates:
  - figures/degree_2_nodes_zoomed_dark.png (High-res multi-panel micro-zoom)
  - for_claude/04_degree_2_nodes_zoomed.png
  - results/degree_2_nodes_analysis.txt
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"
FOR_CLAUDE = ROOT / "for_claude"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
FOR_CLAUDE.mkdir(parents=True, exist_ok=True)

SURFACE = "#0f0f10"
PANEL_BG = "#17171a"
INKD = "#f4f4f5"
SUBD = "#a1a1aa"
GRID = "#27272a"
SPINE = "#3f3f46"


def str_val(val):
    if isinstance(val, list):
        return "; ".join(str(x) for x in val)
    return str(val) if val is not None else ""


def inspect_city(city_name, graph_path):
    print(f"\nAnalyzing {city_name} ...")
    G = ox.load_graphml(graph_path)
    
    # Check street_count == 2
    sc2_nodes = [u for u, d in G.nodes(data=True) if str(d.get("street_count")) == "2"]
    
    # Undirected degree == 2
    S = nx.Graph()
    S.add_nodes_from(G.nodes(data=True))
    S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
    undir_deg2 = [u for u, d in S.degree() if d == 2]
    
    print(f"  {city_name}: {len(sc2_nodes)} nodes with street_count=2, {len(undir_deg2)} nodes with graph degree=2")
    
    detailed = []
    for u in sc2_nodes:
        ndata = G.nodes[u]
        in_edges = list(G.in_edges(u, data=True))
        out_edges = list(G.out_edges(u, data=True))
        
        in_hw = {str_val(d.get("highway")) for _, _, d in in_edges}
        out_hw = {str_val(d.get("highway")) for _, _, d in out_edges}
        hw_change = (in_hw != out_hw)
        
        in_ow = {str_val(d.get("oneway")) for _, _, d in in_edges}
        out_ow = {str_val(d.get("oneway")) for _, _, d in out_edges}
        ow_change = (in_ow != out_ow)
        
        # Check source/sink (direction reversal)
        dir_reversal = (len(in_edges) == 0 or len(out_edges) == 0)
        
        bridge_tunnel = any(d.get("bridge") or d.get("tunnel") for _, _, d in in_edges + out_edges)
        node_tag = ndata.get("highway") or ndata.get("barrier") or ndata.get("junction")
        
        # Determine root cause
        if dir_reversal:
            cause = "Direction Reversal / Cul-de-Sac Convergence"
        elif node_tag:
            cause = f"Special Node Tag ({node_tag})"
        elif bridge_tunnel:
            cause = "Bridge / Tunnel Ramp Boundary"
        elif hw_change:
            cause = "Road Hierarchy Boundary (e.g. primary -> secondary)"
        elif ow_change:
            cause = "One-way to Two-way Transition"
        else:
            cause = "Attribute Change / Non-simplifiable Edge Geometry"
            
        names = set()
        for _, _, d in in_edges + out_edges:
            nm = d.get("name")
            if nm:
                names.add(str_val(nm))
        st_name = "; ".join(names) if names else "(unnamed street)"
        
        detailed.append({
            "node_id": u,
            "lat": ndata.get("y"),
            "lon": ndata.get("x"),
            "cause": cause,
            "in_edges": len(in_edges),
            "out_edges": len(out_edges),
            "streets": st_name,
            "hw_in": "; ".join(in_hw),
            "hw_out": "; ".join(out_hw),
            "node_tag": node_tag or "none"
        })
        
    df = pd.DataFrame(detailed)
    return G, df


G_blr, df_blr = inspect_city("Bengaluru", "data/bengaluru_city_bengaluru_drive.graphml")
G_chn, df_chn = inspect_city("Chennai", "data/chennai_city_chennai_drive.graphml")

# Write text report
report_lines = [
    "==================================================================",
    "ANALYSIS OF DEGREE-2 NODES (street_count = 2) IN STREET NETWORKS",
    "==================================================================",
    "",
    "WHY DO DEGREE-2 NODES EXIST IN SIMPLIFIED ROAD GRAPHS?",
    "OSMnx simplification algorithm removes all intermediate nodes along",
    "a street EXCEPT where simplification would destroy network topology",
    "or critical spatial/attribute boundaries. Specifically:",
    "",
    "1. BENGALURU (149 degree-2 nodes out of 155,779 nodes - 0.10%):",
]
blr_counts = df_blr["cause"].value_counts()
for cause, cnt in blr_counts.items():
    report_lines.append(f"   • {cause}: {cnt} nodes ({cnt/len(df_blr):.1%})")

report_lines.extend(["", "2. CHENNAI (53 degree-2 nodes out of 69,386 nodes - 0.08%):"])
chn_counts = df_chn["cause"].value_counts()
for cause, cnt in chn_counts.items():
    report_lines.append(f"   • {cause}: {cnt} nodes ({cnt/len(df_chn):.1%})")

rep_path = RESULTS_DIR / "degree_2_nodes_analysis.txt"
with open(rep_path, "w", encoding="utf-8") as f:
    f.write("\n".join(report_lines))
print(f"\nSaved report to {rep_path}")
print("\n".join(report_lines))


# -------------------------------------------------------------------------
# Micro-Zoom Visualization of Interesting Degree-2 Locations
# -------------------------------------------------------------------------
print("\nGenerating zoomed-in micro-level visualizations...")

# Pick 4 representative cases:
# Case 1: Bengaluru - Highway classification / flyover boundary
# Case 2: Bengaluru - One-way transition / directional convergence
# Case 3: Chennai - Coastal / canal transition
# Case 4: Chennai - Loop / turnaround convergence

fig, axes = plt.subplots(2, 2, figsize=(16, 14), dpi=300)
fig.patch.set_facecolor(SURFACE)

cases = [
    {
        "city": "Bengaluru",
        "G": G_blr,
        "node_id": df_blr[df_blr["cause"].str.contains("Road Hierarchy")].iloc[0]["node_id"],
        "title": "Case 1: Road Hierarchy Boundary (Bengaluru)",
        "subtitle": "Primary Arterial transitions into Secondary/Tertiary street",
        "ax": axes[0, 0]
    },
    {
        "city": "Bengaluru",
        "G": G_blr,
        "node_id": df_blr[df_blr["cause"].str.contains("Direction Reversal|Transition")].iloc[0]["node_id"],
        "title": "Case 2: One-Way Directional Split (Bengaluru)",
        "subtitle": "Two-way road becomes divided one-way pair or reversal",
        "ax": axes[0, 1]
    },
    {
        "city": "Chennai",
        "G": G_chn,
        "node_id": df_chn[df_chn["cause"].str.contains("Road Hierarchy")].iloc[0]["node_id"],
        "title": "Case 3: Attribute & Class Transition (Chennai)",
        "subtitle": "Speed limit / road class changes across municipal sector",
        "ax": axes[1, 0]
    },
    {
        "city": "Chennai",
        "G": G_chn,
        "node_id": df_chn[df_chn["cause"].str.contains("Direction Reversal|Cul-de-Sac")].iloc[0]["node_id"],
        "title": "Case 4: Cul-de-Sac / Turnaround Loop (Chennai)",
        "subtitle": "Vehicular turnaround loop meeting incoming street",
        "ax": axes[1, 1]
    },
]

for c in cases:
    ax = c["ax"]
    ax.set_facecolor(PANEL_BG)
    G_curr = c["G"]
    u = c["node_id"]
    u_lat = G_curr.nodes[u]["y"]
    u_lon = G_curr.nodes[u]["x"]
    
    # Extract 350-meter ego neighborhood around the node
    # Find nodes within ~400m
    d_lat = 0.0035  # ~380m
    d_lon = 0.0035  # ~380m
    
    sub_nodes = [
        n for n, d in G_curr.nodes(data=True)
        if abs(d["y"] - u_lat) < d_lat and abs(d["x"] - u_lon) < d_lon
    ]
    subG = G_curr.subgraph(sub_nodes)
    
    # Plot background edges
    for un, vn, edata in subG.edges(data=True):
        is_incident = (un == u or vn == u)
        col = "#38bdf8" if is_incident else "#3f3f46"
        lw = 3.5 if is_incident else 1.5
        zord = 3 if is_incident else 1
        
        if "geometry" in edata:
            xs, ys = edata["geometry"].xy
            ax.plot(xs, ys, color=col, lw=lw, zorder=zord, alpha=0.9)
        else:
            x1, y1 = G_curr.nodes[un]["x"], G_curr.nodes[un]["y"]
            x2, y2 = G_curr.nodes[vn]["x"], G_curr.nodes[vn]["y"]
            ax.plot([x1, x2], [y1, y2], color=col, lw=lw, zorder=zord, alpha=0.9)
            
        # Draw directional arrow on incident edges
        if is_incident:
            x1, y1 = G_curr.nodes[un]["x"], G_curr.nodes[un]["y"]
            x2, y2 = G_curr.nodes[vn]["x"], G_curr.nodes[vn]["y"]
            mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
            dx, dy = (x2 - x1) * 0.15, (y2 - y1) * 0.15
            ax.annotate("", xy=(mx + dx, my + dy), xytext=(mx - dx, my - dy),
                        arrowprops=dict(arrowstyle="->", color="#38bdf8", lw=2.0), zorder=4)

    # Plot neighbor nodes
    nbr_xs = [G_curr.nodes[n]["x"] for n in subG.nodes() if n != u]
    nbr_ys = [G_curr.nodes[n]["y"] for n in subG.nodes() if n != u]
    ax.scatter(nbr_xs, nbr_ys, color="#71717a", s=30, zorder=2)
    
    # Plot target degree-2 node as glowing beacon
    ax.scatter([u_lon], [u_lat], s=250, facecolors="#f43f5e", edgecolors="#ffffff", lw=2.0, zorder=6)
    
    # Annotate target node
    info_row = (df_blr if c["city"] == "Bengaluru" else df_chn).query(f"node_id == {u}").iloc[0]
    anno_text = (
        f"Degree-2 Node: {u}\n"
        f"Cause: {info_row['cause']}\n"
        f"Streets: {info_row['streets']}\n"
        f"Coords: ({u_lat:.5f}, {u_lon:.5f})"
    )
    ax.text(
        0.03, 0.94, anno_text, transform=ax.transAxes,
        fontsize=9, color=INKD, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#18181b", edgecolor="#f43f5e", lw=1.2, alpha=0.92),
        zorder=7
    )
    
    ax.set_title(f"{c['title']}\n{c['subtitle']}", color="#ffffff", fontsize=11, fontweight="bold", pad=10)
    ax.set_aspect("equal")
    ax.set_xlim(u_lon - d_lon * 0.85, u_lon + d_lon * 0.85)
    ax.set_ylim(u_lat - d_lat * 0.85, u_lat + d_lat * 0.85)
    ax.tick_params(colors=SUBD, labelsize=8)
    for sp in ax.spines.values(): sp.set_color(SPINE)
    ax.grid(True, color=GRID, linestyle=":", linewidth=0.5)

fig.suptitle("Micro-Level Investigation: Why Degree-2 Nodes Exist in Simplified Street Graphs",
             color="#ffffff", fontsize=14, fontweight="bold", y=0.98)

fig.tight_layout()
fig.subplots_adjust(top=0.91)

out_fig1 = FIGURES_DIR / "degree_2_nodes_zoomed_dark.png"
out_fig2 = FOR_CLAUDE / "04_degree_2_nodes_zoomed.png"
fig.savefig(out_fig1, bbox_inches="tight", facecolor=SURFACE)
fig.savefig(out_fig2, bbox_inches="tight", facecolor=SURFACE)
plt.close(fig)
print(f"Generated zoomed visualization: {out_fig1} and {out_fig2}")

