# Author: Kanak Varma
"""
Publication-Quality Side-by-Side Comparison for a Degree-2 Node in Dense Urban Context:
Left: Rich OSM Street Network Grid showing surrounding topology and the degree-2 breakpoint.
Right: High-Resolution Satellite Ground Truth showing matching urban blocks and street layout.
"""

import math
import urllib.request
from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt
from PIL import Image
from shapely import wkt
from shapely.geometry import box, LineString

PROJECT_ROOT = Path("/home/kanak/Documents/University/ID5080-complex_networks/project")
FIGURES_DIR = PROJECT_ROOT / "figures"
FOR_CLAUDE_DIR = PROJECT_ROOT / "for_claude"
CACHE_DIR = PROJECT_ROOT / "data" / "satellite_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

TARGET_NODE = "314698438"
GRAPH_PATH = PROJECT_ROOT / "data" / "bengaluru_city_bengaluru_drive.graphml"

def deg2num(lat_deg, lon_deg, zoom):
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
    return xtile, ytile

def num2deg(xtile, ytile, zoom):
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return lat_deg, lon_deg

def fetch_tile(xtile, ytile, zoom):
    tile_file = CACHE_DIR / f"{zoom}_{xtile}_{ytile}.jpg"
    if tile_file.exists():
        return Image.open(tile_file).convert("RGB")
    url = f"https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{zoom}/{ytile}/{xtile}.jpg"
    headers = {"User-Agent": "Mozilla/5.0"}
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=10) as resp:
        content = resp.read()
        with open(tile_file, "wb") as f:
            f.write(content)
    return Image.open(tile_file).convert("RGB")

def get_satellite_extent(center_lat, center_lon, zoom=17, radius_tiles=1):
    cx, cy = deg2num(center_lat, center_lon, zoom)
    tiles_x = range(cx - radius_tiles, cx + radius_tiles + 1)
    tiles_y = range(cy - radius_tiles, cy + radius_tiles + 1)
    
    w = len(tiles_x) * 256
    h = len(tiles_y) * 256
    stitched = Image.new("RGB", (w, h))
    
    for i, tx in enumerate(tiles_x):
        for j, ty in enumerate(tiles_y):
            tile = fetch_tile(tx, ty, zoom)
            stitched.paste(tile, (i * 256, j * 256))
            
    nw_lat, nw_lon = num2deg(tiles_x[0], tiles_y[0], zoom)
    se_lat, se_lon = num2deg(tiles_x[-1] + 1, tiles_y[-1] + 1, zoom)
    return stitched, [nw_lon, se_lon, se_lat, nw_lat]

def main():
    print(f"Loading graph from {GRAPH_PATH}...")
    G = nx.read_graphml(GRAPH_PATH)
    
    nd = G.nodes[TARGET_NODE]
    t_lat = float(nd["y"])
    t_lon = float(nd["x"])
    print(f"Target Degree-2 Node: {TARGET_NODE} at ({t_lat:.6f}, {t_lon:.6f})")
    
    # Urban field of view: ~350m radius (700m × 700m urban block context)
    d_lat = 0.0032
    d_lon = 0.0032
    bbox = [t_lon - d_lon, t_lon + d_lon, t_lat - d_lat, t_lat + d_lat]
    bbox_geom = box(bbox[0], bbox[2], bbox[1], bbox[3])
    
    # 1. Fetch satellite tiles
    print("Fetching satellite imagery...")
    sat_img, sat_extent = get_satellite_extent(t_lat, t_lon, zoom=17, radius_tiles=1)
    
    # 2. Extract intersecting edges & nodes within the urban grid
    intersecting_edges = []
    seen_nodes = set()
    for u, v, d in G.edges(data=True):
        geom_str = d.get("geometry")
        if geom_str:
            line = wkt.loads(geom_str)
        else:
            p1 = (float(G.nodes[u]["x"]), float(G.nodes[u]["y"]))
            p2 = (float(G.nodes[v]["x"]), float(G.nodes[v]["y"]))
            line = LineString([p1, p2])
            
        if line.intersects(bbox_geom):
            intersecting_edges.append((u, v, d, line))
            if bbox[0] <= float(G.nodes[u]["x"]) <= bbox[1] and bbox[2] <= float(G.nodes[u]["y"]) <= bbox[3]:
                seen_nodes.add(u)
            if bbox[0] <= float(G.nodes[v]["x"]) <= bbox[1] and bbox[2] <= float(G.nodes[v]["y"]) <= bbox[3]:
                seen_nodes.add(v)
                
    print(f"Urban street network context: {len(intersecting_edges)} edges, {len(seen_nodes)} nodes")
    
    # Color palette: Minimal, high contrast
    BG_COLOR = "#0B0D11"
    PANEL_BG = "#13171F"
    TEXT_MAIN = "#F1F5F9"
    TEXT_MUTED = "#94A3B8"
    COLOR_PRIMARY = "#38BDF8"    # Sky Blue: Arterial / Dhanvanthari Road
    COLOR_SECONDARY = "#818CF8"  # Indigo: Collector / Connecting streets
    COLOR_LOCAL = "#475569"      # Slate: Urban grid / residential
    COLOR_TARGET = "#EF4444"     # Red: Target Degree-2 Node
    
    fig, (ax_graph, ax_sat) = plt.subplots(1, 2, figsize=(15, 7.2), dpi=300)
    fig.patch.set_facecolor(BG_COLOR)
    
    # -------------------------------------------------------------
    # PANEL 1: URBAN STREET NETWORK GRAPH
    # -------------------------------------------------------------
    ax_graph.set_facecolor(PANEL_BG)
    
    for u, v, d, line in intersecting_edges:
        xs, ys = line.xy
        hw = str(d.get("highway"))
        name = str(d.get("name") or "")
        is_incident = (u == TARGET_NODE or v == TARGET_NODE)
        
        if is_incident:
            col = "#F59E0B" if (d.get("oneway") == "True" or d.get("oneway") == True) else COLOR_PRIMARY
            lw = 4.2
            zorder = 7
        elif "primary" in hw or "Dhanvanthari" in name:
            col = COLOR_PRIMARY
            lw = 2.8
            zorder = 5
        elif "secondary" in hw or "tertiary" in hw:
            col = COLOR_SECONDARY
            lw = 2.0
            zorder = 4
        else:
            col = COLOR_LOCAL
            lw = 1.2
            zorder = 2
            
        ax_graph.plot(xs, ys, color=col, linewidth=lw, alpha=0.9, zorder=zorder)
        
        # Flow arrows on key incident edges
        if is_incident:
            mx = (xs[0] + xs[-1]) / 2.0
            my = (ys[0] + ys[-1]) / 2.0
            dx = (xs[-1] - xs[0]) * 0.15
            dy = (ys[-1] - ys[0]) * 0.15
            ax_graph.annotate(
                "", xy=(mx + dx, my + dy), xytext=(mx - dx, my - dy),
                arrowprops=dict(arrowstyle="->", color=col, lw=2.2, mutation_scale=14),
                zorder=8
            )
            
    # Regular urban junctions
    other_xs = [float(G.nodes[n]["x"]) for n in seen_nodes if n != TARGET_NODE]
    other_ys = [float(G.nodes[n]["y"]) for n in seen_nodes if n != TARGET_NODE]
    ax_graph.scatter(other_xs, other_ys, c="#64748B", s=22, zorder=4, alpha=0.75, edgecolors="none")
    
    # Target Degree-2 Node highlight
    ax_graph.scatter(t_lon, t_lat, c=COLOR_TARGET, s=260, zorder=10, edgecolors="white", linewidth=2.2)
    ax_graph.scatter(t_lon, t_lat, c=COLOR_TARGET, s=650, zorder=9, alpha=0.35)
    
    # Clean, concise annotations
    ax_graph.annotate(
        f"Degree-2 Node ({TARGET_NODE})\n"
        "1-Way $\\rightarrow$ 2-Way Transition\n"
        "(Simplification blocked)",
        xy=(t_lon, t_lat), xytext=(t_lon - 0.0018, t_lat + 0.0010),
        color=TEXT_MAIN, fontsize=9.5, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#1E293B", edgecolor=COLOR_TARGET, lw=1.5, alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="white", lw=1.5),
        zorder=12
    )
    
    ax_graph.annotate(
        "1-Way Inflow (Amber)",
        xy=(t_lon + 0.0006, t_lat + 0.0001),
        color="#F59E0B", fontsize=9, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#0F172A", edgecolor="#F59E0B", lw=1),
        zorder=12
    )
    
    ax_graph.annotate(
        "2-Way Road (Cyan)",
        xy=(t_lon - 0.0007, t_lat - 0.0005),
        color=COLOR_PRIMARY, fontsize=9, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#0F172A", edgecolor=COLOR_PRIMARY, lw=1),
        zorder=12
    )

    ax_graph.set_xlim(bbox[0], bbox[1])
    ax_graph.set_ylim(bbox[2], bbox[3])
    ax_graph.set_title("OSM Street Network Graph", color=TEXT_MAIN, fontsize=12, fontweight="bold", pad=10)
    ax_graph.tick_params(colors=TEXT_MUTED, labelsize=8)
    for spine in ax_graph.spines.values():
        spine.set_color("#2D3748")
    ax_graph.grid(True, color="#1E293B", linestyle="--", alpha=0.5)

    # -------------------------------------------------------------
    # PANEL 2: SATELLITE IMAGERY
    # -------------------------------------------------------------
    ax_sat.set_facecolor(PANEL_BG)
    ax_sat.imshow(sat_img, extent=sat_extent, aspect="auto")
    
    # Overlay target reticle
    ax_sat.scatter(t_lon, t_lat, c="none", s=340, edgecolors=COLOR_TARGET, linewidth=2.8, zorder=10)
    ax_sat.scatter(t_lon, t_lat, c=COLOR_TARGET, s=60, edgecolors="white", linewidth=1.5, zorder=11)
    
    cross_d = 0.0002
    ax_sat.plot([t_lon - cross_d, t_lon + cross_d], [t_lat, t_lat], color="white", lw=1.2, zorder=10)
    ax_sat.plot([t_lon, t_lon], [t_lat - cross_d, t_lat + cross_d], color="white", lw=1.2, zorder=10)

    # Clean satellite ground-truth callout
    ax_sat.annotate(
        "Dhanvanthari Road Corridor\n(Urban Grid Alignment)",
        xy=(t_lon, t_lat), xytext=(t_lon - 0.0019, t_lat + 0.0012),
        color="white", fontsize=9.5, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="black", edgecolor=COLOR_TARGET, lw=1.5, alpha=0.85),
        arrowprops=dict(arrowstyle="->", color=COLOR_TARGET, lw=1.8),
        zorder=12
    )

    ax_sat.set_xlim(bbox[0], bbox[1])
    ax_sat.set_ylim(bbox[2], bbox[3])
    ax_sat.set_title("Satellite Imagery", color=TEXT_MAIN, fontsize=12, fontweight="bold", pad=10)
    ax_sat.tick_params(colors=TEXT_MUTED, labelsize=8)
    for spine in ax_sat.spines.values():
        spine.set_color("#2D3748")

    # Short single-line title
    fig.suptitle("Degree-2 Node Investigation: Dhanvanthari Road, Gandhinagar Grid",
                 color=TEXT_MAIN, fontsize=13, fontweight="bold", y=0.97)
    
    plt.subplots_adjust(left=0.06, right=0.96, top=0.90, bottom=0.08, wspace=0.12)
    
    out1 = FIGURES_DIR / "degree_2_sample_comparison_dark.png"
    out2 = FOR_CLAUDE_DIR / "05_degree_2_sample_ground_truth.png"
    
    fig.savefig(out1, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(out2, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close()
    
    print(f"Generated comparison figures:\n  1. {out1}\n  2. {out2}")

if __name__ == "__main__":
    main()
