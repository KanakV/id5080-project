# Author: Kanak Varma
"""
Network Voronoi Analysis for Police and Fire Stations in Bengaluru & Chennai.

Partitions the drivable road network into catchment areas (Network Voronoi cells)
for each station using length-weighted multi-source shortest paths on the road graph.

Produces:
  - Euclidean Voronoi boundaries vs Network road catchments
  - Colored road segments assigned to their closest emergency station (by road distance)
  - Distance-to-nearest service heatmaps
  - Service coverage statistics (road km covered, max response distance per station)

Outputs:
  - data_analysis/figures/<city>_police_network_voronoi_dark.png
  - data_analysis/figures/<city>_fire_station_network_voronoi_dark.png
  - data_analysis/figures/<city>_emergency_combined_voronoi_dark.png
  - data_analysis/results/<city>_emergency_catchment_stats.csv
"""

import csv
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from scipy.spatial import KDTree, Voronoi, Delaunay

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "data_analysis" / "results"
FIGURES_DIR = ROOT / "data_analysis" / "figures"
ROOT_FIGURES_DIR = ROOT / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Deep obsidian dark background palette
SURFACE = "#0f0f10"
INKD = "#f5f5f3"
SUBD = "#b0b0a8"
GRID = "#262629"
ROAD_BG = "#1c1c1f"

LIGHT_SURFACE = "#ffffff"
LIGHT_INK = "#0f172a"
LIGHT_SUB = "#475569"
LIGHT_GRID = "#cbd5e1"
LIGHT_ROAD_BG = "#e2e8f0"

# Curated harmonious neon/jewel-tone palette for adjacent Voronoi catchments
# Designed for maximum visual contrast on deep dark backgrounds (#0f0f10)
HARMONIC_PALETTE = [
    "#38bdf8",  # Sky Blue
    "#f43f5e",  # Rose Red
    "#10b981",  # Vivid Emerald
    "#fbbf24",  # Warm Amber
    "#a855f7",  # Bright Purple
    "#06b6d4",  # Cyan
    "#f97316",  # Tangerine Orange
    "#ec4899",  # Hot Pink
    "#22c55e",  # Vibrant Green
    "#818cf8",  # Soft Indigo
    "#eab308",  # Deep Gold
    "#14b8a6",  # Teal
]

LIGHT_HARMONIC_PALETTE = [
    "#0284c7",  # Cerulean Blue
    "#e11d48",  # Crimson Rose
    "#059669",  # Emerald
    "#d97706",  # Rich Amber
    "#7c3aed",  # Violet
    "#0891b2",  # Deep Cyan
    "#ea580c",  # Burnt Orange
    "#db2777",  # Magenta
    "#16a34a",  # Forest Green
    "#4f46e5",  # Indigo
    "#b45309",  # Ochre
    "#0d9488",  # Teal
]


def assign_adjacent_cell_colors(coords, palette=HARMONIC_PALETTE):
    """
    Uses Delaunay triangulation (the planar dual of Voronoi diagrams) to construct
    the exact adjacency graph of cells, then applies greedy graph coloring so that
    no two adjacent catchments share the same color.
    """
    n = len(coords)
    if n <= 3:
        return [palette[i % len(palette)] for i in range(n)]

    try:
        tri = Delaunay(coords)
        adj_graph = nx.Graph()
        for i in range(n):
            adj_graph.add_node(i)
        for simplex in tri.simplices:
            for i in range(3):
                for j in range(i + 1, 3):
                    adj_graph.add_edge(simplex[i], simplex[j])

        coloring = nx.coloring.greedy_color(adj_graph, strategy="largest_first")
        assigned_colors = [palette[coloring[i] % len(palette)] for i in range(n)]
        return assigned_colors
    except Exception:
        # Fallback if triangulation degenerates
        return [palette[i % len(palette)] for i in range(n)]

CITIES = {
    "bengaluru": {
        "graphml": DATA_DIR / "bengaluru_city_bengaluru_drive.graphml",
        "csv": DATA_DIR / "bengaluru_emergency_services.csv",
        "title": "Bengaluru",
    },
    "chennai": {
        "graphml": DATA_DIR / "chennai_city_chennai_drive.graphml",
        "csv": DATA_DIR / "chennai_emergency_services.csv",
        "title": "Chennai",
    },
}


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_stations(csv_path: Path, category: str):
    stations = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["category"] == category:
                try:
                    lat = float(row["lat"])
                    lon = float(row["lon"])
                    stations.append({
                        "name": row["name"],
                        "lat": lat,
                        "lon": lon,
                        "id": row["osm_id"],
                    })
                except (ValueError, TypeError):
                    continue
    return stations


def compute_network_voronoi(G, stations):
    """
    Snap stations to nearest graph nodes using scipy KDTree, then run multi-source Dijkstra.
    Returns:
      node_to_station: dict mapping each node in G -> station index
      node_distance: dict mapping each node in G -> shortest road distance (meters)
      station_to_node: list of snapped node ids
    """
    graph_nodes = list(G.nodes())
    node_coords = np.array([(G.nodes[n]["x"], G.nodes[n]["y"]) for n in graph_nodes])
    kdtree = KDTree(node_coords)

    st_coords = np.array([(s["lon"], s["lat"]) for s in stations])
    _, nearest_idx = kdtree.query(st_coords)
    snapped_nodes = [graph_nodes[i] for i in nearest_idx]

    node_to_station_idx = {}
    valid_sources = []
    for idx, node in enumerate(snapped_nodes):
        if node in G:
            valid_sources.append(node)
            node_to_station_idx[node] = idx

    # Run multi-source Dijkstra length-weighted
    dist, paths = nx.multi_source_dijkstra(G, sources=valid_sources, weight="length")

    node_owner = {}
    for node, path in paths.items():
        source_node = path[0]
        node_owner[node] = node_to_station_idx[source_node]

    return node_owner, dist, snapped_nodes


def plot_network_voronoi(city_key, category, G, stations, node_owner, node_dist, snapped_nodes, theme="dark"):
    title_cat = "Police Stations" if category == "police" else "Fire Stations"
    log(f"Plotting Network Voronoi ({theme}) for {CITIES[city_key]['title']} {title_cat} ...")

    if theme == "light":
        surface, ink, sub, grid, road_bg = LIGHT_SURFACE, LIGHT_INK, LIGHT_SUB, LIGHT_GRID, LIGHT_ROAD_BG
        palette = LIGHT_HARMONIC_PALETTE
        ridge_col, ridge_alpha = "#334155", 0.45
        halo_col, pin_col, pin_edge = "#ffffff", "#0f172a", "#ffffff"
        leg_bg = "#f8fafc"
    else:
        surface, ink, sub, grid, road_bg = SURFACE, INKD, SUBD, GRID, ROAD_BG
        palette = HARMONIC_PALETTE
        ridge_col, ridge_alpha = "#ffffff", 0.40
        halo_col, pin_col, pin_edge = "#ffffff", "#ffffff", "#09090b"
        leg_bg = "#242422"

    fig, ax = plt.subplots(figsize=(14, 14), facecolor=surface)
    ax.set_facecolor(surface)

    num_stations = len(stations)
    coords = np.array([[s["lon"], s["lat"]] for s in stations])
    active_station_ids = set(node_owner.values())
    active_count = len(active_station_ids)

    # Assign harmonious colors via Delaunay dual planar graph coloring
    colors = assign_adjacent_cell_colors(coords, palette=palette)

    # Segment lines grouped by owner
    lines_by_station = {i: [] for i in range(num_stations)}
    unreached_lines = []

    for u, v, data in G.edges(data=True):
        if "geometry" in data:
            xs, ys = data["geometry"].xy
            pts = list(zip(xs, ys))
        else:
            pts = [(G.nodes[u]["x"], G.nodes[u]["y"]), (G.nodes[v]["x"], G.nodes[v]["y"])]

        owner_u = node_owner.get(u)
        owner_v = node_owner.get(v)

        if owner_u is not None:
            lines_by_station[owner_u].append(pts)
        elif owner_v is not None:
            lines_by_station[owner_v].append(pts)
        else:
            unreached_lines.append(pts)

    # Background road grid for unreached/peripheral streets
    if unreached_lines:
        lc_un = LineCollection(unreached_lines, colors=road_bg, linewidths=0.5, alpha=0.45, zorder=1)
        ax.add_collection(lc_un)

    # Add each station's road catchment collection with crisp line weighting and high alpha
    for st_idx, segs in lines_by_station.items():
        if segs:
            col = colors[st_idx]
            lc = LineCollection(segs, colors=[col], linewidths=0.9, alpha=0.92, zorder=2)
            ax.add_collection(lc)

    # Subtle Euclidean Voronoi ridge lines
    if len(coords) >= 4:
        try:
            vor = Voronoi(coords)
            for simplex in vor.ridge_vertices:
                simplex = np.asarray(simplex)
                if np.all(simplex >= 0):
                    ax.plot(
                        vor.vertices[simplex, 0],
                        vor.vertices[simplex, 1],
                        color=ridge_col,
                        linestyle=(0, (4, 4)),
                        linewidth=0.85,
                        alpha=ridge_alpha,
                        zorder=3,
                    )
        except Exception as e:
            log(f"Voronoi ridge plotting skipped: {e}")

    # Plot facility origin markers: outer halo + crisp center pin
    marker = "o" if category == "police" else "^"
    st_xs = coords[:, 0]
    st_ys = coords[:, 1]

    # Outer soft halo
    ax.scatter(
        st_xs,
        st_ys,
        s=150 if category == "police" else 230,
        c=halo_col,
        alpha=0.65 if theme == "light" else 0.25,
        marker=marker,
        zorder=5,
    )
    # Inner crisp pin
    ax.scatter(
        st_xs,
        st_ys,
        s=48 if category == "police" else 85,
        c=pin_col,
        edgecolors=pin_edge,
        linewidth=1.2,
        marker=marker,
        zorder=6,
        label=f"{title_cat} ({active_count} active)",
    )

    # Bounds
    all_xs = [d["x"] for _, d in G.nodes(data=True) if "x" in d]
    all_ys = [d["y"] for _, d in G.nodes(data=True) if "y" in d]
    px = (max(all_xs) - min(all_xs)) * 0.02
    py = (max(all_ys) - min(all_ys)) * 0.02
    ax.set_xlim(min(all_xs) - px, max(all_xs) + px)
    ax.set_ylim(min(all_ys) - py, max(all_ys) + py)
    ax.set_aspect("equal")

    # Title & Legend
    count_str = f"n = {active_count} active" if active_count == len(stations) else f"n = {active_count} active ({len(stations)} raw OSM)"
    ax.set_title(
        f"{CITIES[city_key]['title']} — {title_cat} Road Catchments (Network Voronoi)\n"
        f"Partitioned by shortest road distance along drivable streets ({count_str})",
        fontsize=13,
        color=ink,
        weight="bold",
        pad=16,
    )

    legend_elements = [
        Line2D([0], [0], marker=marker, color=pin_col, markeredgecolor=pin_edge,
               label=f"{title_cat} ({active_count} active)", markersize=8, linestyle="None"),
        Line2D([0], [0], color=palette[0], lw=2, label="Road Catchment (closest station via network)"),
        Line2D([0], [0], color=ridge_col, lw=1, linestyle="--", alpha=0.7, label="Euclidean Voronoi boundary"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", facecolor=leg_bg, edgecolor=grid, labelcolor=ink, fontsize=9.5)

    ax.tick_params(colors=sub, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(grid)
    ax.grid(True, color=grid, linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("Longitude", color=sub, fontsize=9)
    ax.set_ylabel("Latitude", color=sub, fontsize=9)

    fname = f"{city_key}_{category}_network_voronoi_{theme}.png"
    plt.tight_layout()
    for out_dir in [FIGURES_DIR, ROOT_FIGURES_DIR]:
        out_file = out_dir / fname
        plt.savefig(out_file, dpi=250, facecolor=surface)
        log(f"Saved: {out_file.relative_to(ROOT)}")
    plt.close()


def compute_and_save_stats(city_key, category, G, stations, node_owner, node_dist):
    stats = []
    # Total street length per station
    station_road_length = {i: 0.0 for i in range(len(stations))}
    station_node_count = {i: 0 for i in range(len(stations))}
    station_max_dist = {i: 0.0 for i in range(len(stations))}

    for node, owner in node_owner.items():
        station_node_count[owner] += 1
        d = node_dist.get(node, 0.0)
        if d > station_max_dist[owner]:
            station_max_dist[owner] = d

    # Edges length
    for u, v, data in G.edges(data=True):
        owner = node_owner.get(u) or node_owner.get(v)
        if owner is not None:
            station_road_length[owner] += data.get("length", 0.0)

    for idx, s in enumerate(stations):
        stats.append({
            "city": city_key,
            "category": category,
            "station_name": s["name"],
            "lat": s["lat"],
            "lon": s["lon"],
            "catchment_nodes": station_node_count[idx],
            "catchment_road_km": round(station_road_length[idx] / 1000.0, 2),
            "max_response_road_km": round(station_max_dist[idx] / 1000.0, 2),
        })
    return stats


def main():
    target = sys.argv[1].lower() if len(sys.argv) > 1 else "all"
    cities = [target] if target in CITIES else ["bengaluru", "chennai"]

    for city in cities:
        cfg = CITIES[city]
        log(f"--- Processing Network Voronoi for {cfg['title']} ---")
        log(f"Loading {cfg['graphml'].name} ...")
        G = ox.load_graphml(str(cfg["graphml"]))

        city_stats = []
        for cat in ["police", "fire_station"]:
            stations = load_stations(cfg["csv"], cat)
            log(f"Loaded {len(stations)} {cat} facilities.")
            if not stations:
                continue

            node_owner, node_dist, snapped_nodes = compute_network_voronoi(G, stations)
            plot_network_voronoi(city, cat, G, stations, node_owner, node_dist, snapped_nodes, theme="dark")
            plot_network_voronoi(city, cat, G, stations, node_owner, node_dist, snapped_nodes, theme="light")
            c_stats = compute_and_save_stats(city, cat, G, stations, node_owner, node_dist)
            city_stats.extend(c_stats)

        # Write CSV stats
        stats_csv = RESULTS_DIR / f"{city}_emergency_catchment_stats.csv"
        with open(stats_csv, "w", newline="", encoding="utf-8") as f:
            fieldnames = [
                "city", "category", "station_name", "lat", "lon",
                "catchment_nodes", "catchment_road_km", "max_response_road_km"
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(city_stats)
        log(f"Saved stats: {stats_csv.relative_to(ROOT)}")

    # Auto-generate 2x2 comparison grid if all 4 images exist
    p_blr = FIGURES_DIR / "bengaluru_police_network_voronoi_dark.png"
    f_blr = FIGURES_DIR / "bengaluru_fire_station_network_voronoi_dark.png"
    p_chn = FIGURES_DIR / "chennai_police_network_voronoi_dark.png"
    f_chn = FIGURES_DIR / "chennai_fire_station_network_voronoi_dark.png"
    if p_blr.exists() and f_blr.exists() and p_chn.exists() and f_chn.exists():
        try:
            import matplotlib.image as mpimg
            fig, axes = plt.subplots(2, 2, figsize=(22, 22), facecolor=SURFACE)
            axes[0, 0].imshow(mpimg.imread(str(p_blr))); axes[0, 0].axis("off")
            axes[0, 1].imshow(mpimg.imread(str(f_blr))); axes[0, 1].axis("off")
            axes[1, 0].imshow(mpimg.imread(str(p_chn))); axes[1, 0].axis("off")
            axes[1, 1].imshow(mpimg.imread(str(f_chn))); axes[1, 1].axis("off")
            plt.tight_layout()
            out_2x2 = FIGURES_DIR / "compare_network_voronoi_2x2_dark.png"
            plt.savefig(out_2x2, dpi=180, facecolor=SURFACE)
            plt.close()
            log(f"Saved 2x2 comparison panel: {out_2x2.relative_to(ROOT)}")
        except Exception as e:
            log(f"Could not build 2x2 comparison panel: {e}")


if __name__ == "__main__":
    main()
