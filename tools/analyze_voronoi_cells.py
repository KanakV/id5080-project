# Author: Kanak Varma
"""
Network Voronoi Numerical Analysis & Subgraph Extraction.

For each emergency service facility (Police & Fire Station):
  1. Identifies the cell catchment (Network Voronoi cell) via multi-source Dijkstra.
  2. Extracts each cell as an independent networkx.MultiDiGraph subgraph.
  3. Computes detailed road length metrics per cell (total length, mean, median, max response road km).
  4. Profiles road types (highway hierarchy: primary, secondary, tertiary, residential, trunk, etc.):
     - Across the entire Voronoi cell
     - Immediately around the cell origin (immediate incident edges and 500m ego network)
  5. Exports:
     - CSV with full numerical breakdown per cell/station
     - Summary report per city
     - Subgraphs saved into data/voronoi_subgraphs/<city>/<category>_cell_<idx>.graphml (optional or indexed bundle)
     - Distribution charts of cell road lengths and road-type compositions

Usage:
  python tools/analyze_voronoi_cells.py chennai
  python tools/analyze_voronoi_cells.py bengaluru
  python tools/analyze_voronoi_cells.py all
"""

import csv
import os
import sys
import time
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
from scipy.spatial import KDTree

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "data_analysis" / "results"
FIGURES_DIR = ROOT / "data_analysis" / "figures"
SUBGRAPHS_DIR = DATA_DIR / "voronoi_subgraphs"

DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
SUBGRAPHS_DIR.mkdir(parents=True, exist_ok=True)

SURFACE = "#1a1a19"
INKD = "#f0efec"
SUBD = "#c3c2b7"
GRID = "#33332f"

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

CATEGORIES = ["police", "fire_station"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_highway(val):
    """Normalize OSM highway tags (handles lists, e.g. ['residential', 'unclassified'])."""
    if isinstance(val, list):
        return val[0]
    if isinstance(val, str) and val.startswith("[") and val.endswith("]"):
        val = val.strip("[]'\" ")
        return val.split(",")[0].strip("'\" ")
    return str(val) if val else "unclassified"


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
                        "osm_id": row["osm_id"],
                    })
                except (ValueError, TypeError):
                    continue
    return stations


def partition_and_analyze(city_key: str, category: str, G, stations, save_graphmls: bool = False):
    cfg = CITIES[city_key]
    title_cat = "Police" if category == "police" else "Fire Station"
    log(f"Partitioning network for {cfg['title']} {title_cat} ({len(stations)} stations)...")

    # Snap stations to nearest road nodes using KDTree
    graph_nodes = list(G.nodes())
    node_coords = np.array([(G.nodes[n]["x"], G.nodes[n]["y"]) for n in graph_nodes])
    kdtree = KDTree(node_coords)

    st_coords = np.array([(s["lon"], s["lat"]) for s in stations])
    _, nearest_indices = kdtree.query(st_coords)

    # Handle co-located or duplicate snaps by grouping stations sharing the same snap node
    # Each station gets an index 0..N-1
    node_to_station_idx = {}
    valid_sources = []
    snapped_nodes = []

    for idx, node_idx in enumerate(nearest_indices):
        node = graph_nodes[node_idx]
        snapped_nodes.append(node)
        if node not in node_to_station_idx:
            node_to_station_idx[node] = idx
            valid_sources.append(node)

    # Multi-source Dijkstra (length-weighted)
    log("Running multi-source Dijkstra length-weighted ...")
    dist, paths = nx.multi_source_dijkstra(G, sources=valid_sources, weight="length")

    # Map each reachable node to the station index
    node_owner = {}
    for node, path in paths.items():
        src_node = path[0]
        node_owner[node] = node_to_station_idx[src_node]

    # Group graph nodes by station index
    station_nodes = {i: [] for i in range(len(stations))}
    for node, owner_idx in node_owner.items():
        station_nodes[owner_idx].append(node)

    # Output directory for subgraphs
    city_subgraph_dir = SUBGRAPHS_DIR / city_key / category
    if save_graphmls:
        city_subgraph_dir.mkdir(parents=True, exist_ok=True)

    cell_records = []

    log(f"Extracting subgraphs and computing road profile for each cell ...")
    for idx, st in enumerate(stations):
        nodes_in_cell = set(station_nodes[idx])
        origin_node = snapped_nodes[idx]

        # Extract cell subgraph
        subG = G.subgraph(nodes_in_cell).copy()

        # Road metrics in cell
        total_road_length_m = 0.0
        edge_count = subG.number_of_edges()
        node_count = subG.number_of_nodes()

        road_type_counter = Counter()
        for u, v, k, data in subG.edges(keys=True, data=True):
            length = float(data.get("length", 0.0))
            total_road_length_m += length
            hwy = parse_highway(data.get("highway", "unclassified"))
            road_type_counter[hwy] += length

        # Max road distance within the cell
        max_dist_m = 0.0
        for n in nodes_in_cell:
            d = dist.get(n, 0.0)
            if d > max_dist_m:
                max_dist_m = d

        # Immediate road environment at origin (incident edges to origin_node)
        origin_incident_types = Counter()
        for u, v, k, data in G.out_edges(origin_node, keys=True, data=True):
            hwy = parse_highway(data.get("highway", "unclassified"))
            origin_incident_types[hwy] += 1
        for u, v, k, data in G.in_edges(origin_node, keys=True, data=True):
            hwy = parse_highway(data.get("highway", "unclassified"))
            origin_incident_types[hwy] += 1

        # 500m Ego-network road types around origin
        origin_ego_types = Counter()
        try:
            ego_500 = nx.ego_graph(G, origin_node, radius=500.0, distance="length")
            for u, v, k, data in ego_500.edges(keys=True, data=True):
                hwy = parse_highway(data.get("highway", "unclassified"))
                origin_ego_types[hwy] += float(data.get("length", 0.0))
        except Exception:
            pass

        # Save subgraph to GraphML if requested
        subgraph_filename = f"{category}_cell_{idx:03d}_{st['osm_id']}.graphml"
        if save_graphmls and node_count > 0:
            ox.save_graphml(subG, filepath=city_subgraph_dir / subgraph_filename)

        # Top 3 road types in cell by length
        top_cell_types = ", ".join([f"{k}: {round(v/1000.0, 1)}km" for k, v in road_type_counter.most_common(3)])
        top_origin_types = ", ".join([f"{k}: {v}" for k, v in origin_incident_types.most_common(3)])
        top_ego_types = ", ".join([f"{k}: {round(v/1000.0, 1)}km" for k, v in origin_ego_types.most_common(3)])

        cell_records.append({
            "cell_id": idx,
            "city": city_key,
            "category": category,
            "station_name": st["name"],
            "osm_id": st["osm_id"],
            "lat": st["lat"],
            "lon": st["lon"],
            "snapped_node_id": origin_node,
            "cell_nodes": node_count,
            "cell_edges": edge_count,
            "total_road_km": round(total_road_length_m / 1000.0, 2),
            "max_response_road_km": round(max_dist_m / 1000.0, 2),
            "primary_km": round(road_type_counter.get("primary", 0.0) / 1000.0, 2),
            "secondary_km": round(road_type_counter.get("secondary", 0.0) / 1000.0, 2),
            "tertiary_km": round(road_type_counter.get("tertiary", 0.0) / 1000.0, 2),
            "residential_km": round(road_type_counter.get("residential", 0.0) / 1000.0, 2),
            "trunk_km": round(road_type_counter.get("trunk", 0.0) / 1000.0, 2),
            "unclassified_km": round(road_type_counter.get("unclassified", 0.0) / 1000.0, 2),
            "top_cell_road_hierarchy": top_cell_types,
            "origin_adjacent_road_types": top_origin_types,
            "origin_500m_road_hierarchy": top_ego_types,
        })

    return cell_records


def plot_numerical_distributions(city_key: str, all_records: list):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), facecolor=SURFACE)
    for row in axes:
        for ax in row:
            ax.set_facecolor(SURFACE)

    police = [r for r in all_records if r["category"] == "police" and r["cell_nodes"] > 0]
    fire = [r for r in all_records if r["category"] == "fire_station" and r["cell_nodes"] > 0]

    # 1. Road length distribution (Police)
    p_lengths = [r["total_road_km"] for r in police]
    axes[0, 0].hist(p_lengths, bins=25, color="#29b6f6", edgecolor="#111110", alpha=0.85)
    axes[0, 0].set_title(f"{CITIES[city_key]['title']} — Police Cell Road Length Distribution", color=INKD, fontsize=11, weight="bold")
    axes[0, 0].set_xlabel("Total Road Length in Cell (km)", color=SUBD, fontsize=9)
    axes[0, 0].set_ylabel("Number of Stations", color=SUBD, fontsize=9)
    mean_p = np.mean(p_lengths) if p_lengths else 0
    axes[0, 0].axvline(mean_p, color="#ffd54f", linestyle="--", linewidth=1.5, label=f"Mean: {mean_p:.1f} km")
    axes[0, 0].legend(facecolor="#242422", edgecolor=GRID, labelcolor=INKD, fontsize=9)

    # 2. Road length distribution (Fire)
    f_lengths = [r["total_road_km"] for r in fire]
    axes[0, 1].hist(f_lengths, bins=15, color="#ffa726", edgecolor="#111110", alpha=0.85)
    axes[0, 1].set_title(f"{CITIES[city_key]['title']} — Fire Station Cell Road Length Distribution", color=INKD, fontsize=11, weight="bold")
    axes[0, 1].set_xlabel("Total Road Length in Cell (km)", color=SUBD, fontsize=9)
    axes[0, 1].set_ylabel("Number of Stations", color=SUBD, fontsize=9)
    mean_f = np.mean(f_lengths) if f_lengths else 0
    axes[0, 1].axvline(mean_f, color="#ffd54f", linestyle="--", linewidth=1.5, label=f"Mean: {mean_f:.1f} km")
    axes[0, 1].legend(facecolor="#242422", edgecolor=GRID, labelcolor=INKD, fontsize=9)

    # 3. Max Response Road Distance Comparison
    axes[1, 0].boxplot(
        [[r["max_response_road_km"] for r in police], [r["max_response_road_km"] for r in fire]],
        tick_labels=["Police Stations", "Fire Stations"],
        patch_artist=True,
        boxprops=dict(facecolor="#29b6f6", color=GRID),
        medianprops=dict(color="#ffffff", linewidth=2),
    )
    axes[1, 0].set_title(f"{CITIES[city_key]['title']} — Max Response Road Distance (km)", color=INKD, fontsize=11, weight="bold")
    axes[1, 0].set_ylabel("Max Shortest Road Distance in Cell (km)", color=SUBD, fontsize=9)

    # 4. Aggregate Road Hierarchy per Category
    r_types = ["primary_km", "secondary_km", "tertiary_km", "residential_km", "trunk_km"]
    labels = ["Primary", "Secondary", "Tertiary", "Residential", "Trunk"]
    p_type_totals = [sum(r[k] for r in police) for k in r_types]
    p_pcts = [v / sum(p_type_totals) * 100 if sum(p_type_totals) > 0 else 0 for v in p_type_totals]

    bars = axes[1, 1].bar(labels, p_pcts, color="#4dd0e1", edgecolor="#111110", alpha=0.9)
    axes[1, 1].set_title(f"{CITIES[city_key]['title']} — Police Catchment Road Type Share (%)", color=INKD, fontsize=11, weight="bold")
    axes[1, 1].set_ylabel("Percentage of Total Road km", color=SUBD, fontsize=9)
    for bar in bars:
        h = bar.get_height()
        axes[1, 1].text(bar.get_x() + bar.get_width() / 2.0, h + 1.0, f"{h:.1f}%", ha="center", va="bottom", color=INKD, fontsize=8)

    for row in axes:
        for ax in row:
            ax.tick_params(colors=SUBD, labelsize=8)
            for sp in ax.spines.values():
                sp.set_color(GRID)
            ax.grid(True, color=GRID, linestyle=":", linewidth=0.4, alpha=0.4)

    plt.tight_layout()
    out_fig = FIGURES_DIR / f"{city_key}_voronoi_numerical_analysis_dark.png"
    plt.savefig(out_fig, dpi=200, facecolor=SURFACE)
    plt.close()
    log(f"Saved figure: {out_fig.relative_to(ROOT)}")


def main():
    target = sys.argv[1].lower() if len(sys.argv) > 1 else "all"
    cities = [target] if target in CITIES else ["bengaluru", "chennai"]

    for city in cities:
        cfg = CITIES[city]
        log(f"=== Running Voronoi Numerical Analysis for {cfg['title']} ===")
        log(f"Loading {cfg['graphml'].name} ...")
        G = ox.load_graphml(str(cfg["graphml"]))

        city_records = []
        for cat in CATEGORIES:
            stations = load_stations(cfg["csv"], cat)
            if not stations:
                continue
            records = partition_and_analyze(city, cat, G, stations, save_graphmls=True)
            city_records.extend(records)

        # Save comprehensive CSV
        out_csv = RESULTS_DIR / f"{city}_voronoi_cells_numerical_analysis.csv"
        fieldnames = [
            "cell_id", "city", "category", "station_name", "osm_id", "lat", "lon",
            "snapped_node_id", "cell_nodes", "cell_edges", "total_road_km", "max_response_road_km",
            "primary_km", "secondary_km", "tertiary_km", "residential_km", "trunk_km", "unclassified_km",
            "top_cell_road_hierarchy", "origin_adjacent_road_types", "origin_500m_road_hierarchy"
        ]
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(city_records)
        log(f"Saved numerical data to {out_csv.relative_to(ROOT)}")

        plot_numerical_distributions(city, city_records)


if __name__ == "__main__":
    main()

