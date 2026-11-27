# Author: Kanak Varma
"""
Comparative Visualization: Full Street Network vs Arterial Backbone (Without Residential)
Shows the morphological transformation when removing 80%+ residential streets:
- Bengaluru: Capillary sponge -> Concentric Radial-Ring Skeleton
- Chennai: Capillary sponge -> Linear Coastal-Radial Funnel
Uses Matplotlib LineCollection for ultra-fast C-level vector rendering.
"""

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import networkx as nx
from pathlib import Path
from shapely import wkt

PROJECT_ROOT = Path("/home/kanak/Documents/University/ID5080-complex_networks/project")
FIGURES_DIR = PROJECT_ROOT / "figures"
FOR_CLAUDE_DIR = PROJECT_ROOT / "for_claude"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
FOR_CLAUDE_DIR.mkdir(parents=True, exist_ok=True)

BG_COLOR = "#0B0D11"
PANEL_BG = "#11141C"
TEXT_MAIN = "#F1F5F9"
TEXT_MUTED = "#94A3B8"

def plot_graph_fast(ax, G, is_backbone_only=False):
    ax.set_facecolor(PANEL_BG)
    
    # Layer definitions:
    # 0: residential/unclassified
    # 1: tertiary
    # 2: secondary
    # 3: primary
    # 4: trunk/motorway
    segments = {0: [], 1: [], 2: [], 3: [], 4: []}
    
    node_coords = {n: (float(d["x"]), float(d["y"])) for n, d in G.nodes(data=True)}
    
    min_x, max_x = float("inf"), float("-inf")
    min_y, max_y = float("inf"), float("-inf")
    
    for u, v, d in G.edges(data=True):
        hw = str(d.get("highway", "")).lower()
        is_res = ("residential" in hw or "living_street" in hw or "unclassified" in hw)
        
        if is_backbone_only and is_res:
            continue
            
        geom_str = d.get("geometry")
        if geom_str:
            line = wkt.loads(geom_str)
            coords = list(line.coords)
        else:
            coords = [node_coords[u], node_coords[v]]
            
        for x, y in coords:
            if x < min_x: min_x = x
            if x > max_x: max_x = x
            if y < min_y: min_y = y
            if y > max_y: max_y = y
            
        if any(t in hw for t in ["trunk", "motorway"]):
            segments[4].append(coords)
        elif "primary" in hw:
            segments[3].append(coords)
        elif "secondary" in hw:
            segments[2].append(coords)
        elif "tertiary" in hw:
            segments[1].append(coords)
        else:
            segments[0].append(coords)
            
    # Render layers via LineCollection
    styles = {
        0: {"color": "#1E293B", "alpha": 0.22, "lw": 0.4, "z": 1},
        1: {"color": "#64748B", "alpha": 0.65, "lw": 0.75, "z": 2},
        2: {"color": "#38BDF8", "alpha": 0.85, "lw": 1.15, "z": 3},
        3: {"color": "#F59E0B", "alpha": 0.95, "lw": 1.5, "z": 4},
        4: {"color": "#EC4899", "alpha": 1.0, "lw": 1.9, "z": 5},
    }
    
    for layer_id in range(5):
        segs = segments[layer_id]
        if not segs:
            continue
        st = styles[layer_id]
        alpha = st["alpha"] if not is_backbone_only else min(1.0, st["alpha"] + 0.15)
        lw = st["lw"] if not is_backbone_only else st["lw"] * 1.3
        
        lc = LineCollection(segs, colors=st["color"], alpha=alpha, linewidths=lw, zorder=st["z"])
        ax.add_collection(lc)
        
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_aspect("equal")
    ax.tick_params(colors=TEXT_MUTED, labelsize=7)
    for spine in ax.spines.values():
        spine.set_color("#1E293B")
    ax.grid(True, color="#1E293B", linestyle=":", alpha=0.35)

def main():
    print("Loading graphs...")
    G_blr = nx.read_graphml(PROJECT_ROOT / "data" / "bengaluru_city_bengaluru_drive.graphml")
    G_chn = nx.read_graphml(PROJECT_ROOT / "data" / "chennai_city_chennai_drive.graphml")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 15), dpi=300)
    fig.patch.set_facecolor(BG_COLOR)
    
    # -------------------------------------------------------------
    # ROW 1: BENGALURU
    # -------------------------------------------------------------
    print("Plotting Bengaluru full network...")
    plot_graph_fast(axes[0, 0], G_blr, is_backbone_only=False)
    axes[0, 0].set_title("Bengaluru: Full Network (All Roads)\n155,779 Nodes | 394,401 Edges | 24,261 km",
                         color=TEXT_MAIN, fontsize=11, fontweight="bold", pad=8)
    
    print("Plotting Bengaluru arterial backbone...")
    plot_graph_fast(axes[0, 1], G_blr, is_backbone_only=True)
    axes[0, 1].set_title("Bengaluru: Arterial Backbone (Excl. Residential)\n34,022 Nodes (21.8%) | 61,314 Edges (15.5%) | 4,169 km",
                         color=TEXT_MAIN, fontsize=11, fontweight="bold", pad=8)
    
    # -------------------------------------------------------------
    # ROW 2: CHENNAI
    # -------------------------------------------------------------
    print("Plotting Chennai full network...")
    plot_graph_fast(axes[1, 0], G_chn, is_backbone_only=False)
    axes[1, 0].set_title("Chennai: Full Network (All Roads)\n69,386 Nodes | 175,689 Edges | 12,410 km",
                         color=TEXT_MAIN, fontsize=11, fontweight="bold", pad=8)
    
    print("Plotting Chennai arterial backbone...")
    plot_graph_fast(axes[1, 1], G_chn, is_backbone_only=True)
    axes[1, 1].set_title("Chennai: Arterial Backbone (Excl. Residential)\n15,624 Nodes (22.5%) | 27,394 Edges (15.6%) | 2,224 km",
                         color=TEXT_MAIN, fontsize=11, fontweight="bold", pad=8)
    
    fig.suptitle("Impact of Removing Residential Roads: Capillary Mesh vs Arterial Backbone",
                 color=TEXT_MAIN, fontsize=14, fontweight="bold", y=0.98)
    
    plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.04, wspace=0.10, hspace=0.15)
    
    out1 = FIGURES_DIR / "network_backbone_comparison_dark.png"
    out2 = FOR_CLAUDE_DIR / "06_network_backbone_comparison.png"
    
    fig.savefig(out1, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(out2, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close()
    
    print(f"Generated comparison figures:\n  1. {out1}\n  2. {out2}")

if __name__ == "__main__":
    main()
