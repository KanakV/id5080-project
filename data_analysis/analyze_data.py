# Author: Kanak Varma
"""
OSM GraphML Data Profile and Completeness Analysis.

Inspects all node, edge, and graph-level attributes in the street network GraphML:
  1. Attribute Completeness: known vs unknown/missing data coverage across all tags
  2. Road Hierarchy: distribution of street classifications by count and total km
  3. Node Topology & Features: junction degrees (T-junctions, crossroads, dead ends)
     and traffic infrastructure (traffic signals, roundabouts, crossings)
  4. Physical Road Properties: one-way vs two-way, lanes, flyovers/bridges, tunnels,
     speed limits, and street naming coverage

Outputs:
  data_analysis/results/
    - <tag>_data_profile_report.txt
    - <tag>_attribute_completeness.csv
    - <tag>_road_hierarchy.csv
    - <tag>_node_features.csv
    - <tag>_road_attributes_summary.csv
  data_analysis/figures/
    - <tag>_attribute_completeness_dark.png
    - <tag>_road_hierarchy_distribution_dark.png
    - <tag>_node_features_and_junctions_dark.png
    - <tag>_road_attributes_profile_dark.png

Usage:
  python data_analysis/analyze_data.py [path/to/graph.graphml]
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
import numpy as np
import osmnx as ox

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

GRAPH_PATH = sys.argv[1] if len(sys.argv) > 1 else str(PROJECT_ROOT / "data" / "bengaluru_city_bengaluru_drive.graphml")
TAG = Path(GRAPH_PATH).stem
CITY = TAG.split("_")[0].capitalize()

# ── Styling Palette ────────────────────────────────────────────────────────────
SURFACE = "#1a1a19"
INKD = "#f0efec"
SUBD = "#c3c2b7"
GRID = "#33332f"
SPINE = "#44443f"
BLUE = "#3987e5"
ORANGE = "#ff8a3d"
GREEN = "#199e70"
RED = "#e53935"
AMBER = "#ffb300"
PURPLE = "#9c27b0"
CYAN = "#00b4d8"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def style_ax(ax, title=None, xlabel=None, ylabel=None):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=9)
    for s in ax.spines.values():
        s.set_color(SPINE)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    if title:
        ax.set_title(title, color=INKD, fontsize=12, pad=8, fontweight="bold")
    if xlabel:
        ax.set_xlabel(xlabel, color=SUBD, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, color=SUBD, fontsize=10)


def normalize_val(val):
    if isinstance(val, list):
        return val[0]
    return val


def main():
    log(f"Loading {GRAPH_PATH} ...")
    t0 = time.time()
    G = ox.load_graphml(GRAPH_PATH)
    n_nodes = len(G.nodes)
    n_edges = len(G.edges)
    log(f"Loaded in {time.time() - t0:.1f}s: {n_nodes:,} nodes, {n_edges:,} directed edges")

    # ───────────────────────────────────────────────────────────────────────────
    # 1. Attribute Completeness
    # ───────────────────────────────────────────────────────────────────────────
    log("Computing attribute completeness ...")
    
    # Expected/standard OSMnx attributes
    node_keys_expected = ["x", "y", "street_count", "highway", "junction", "ref"]
    edge_keys_expected = [
        "osmid", "highway", "oneway", "reversed", "length", "name", "geometry",
        "lanes", "maxspeed", "access", "bridge", "ref", "width", "junction",
        "tunnel", "est_width", "area"
    ]

    node_counts = Counter()
    for _, d in G.nodes(data=True):
        for k in d:
            node_counts[k] += 1

    edge_counts = Counter()
    edge_length_by_key = Counter()
    total_network_len_m = 0.0

    for _, _, d in G.edges(data=True):
        length_m = float(d.get("length", 0.0))
        total_network_len_m += length_m
        for k in d:
            edge_counts[k] += 1
            edge_length_by_key[k] += length_m

    # Build completeness records
    node_completeness = []
    for k in sorted(set(list(node_counts.keys()) + node_keys_expected)):
        cnt = node_counts[k]
        pct = (cnt / n_nodes) * 100 if n_nodes else 0
        node_completeness.append({
            "target": "node", "attribute": k, "present": cnt,
            "missing": n_nodes - cnt, "pct_present": pct, "total": n_nodes
        })

    edge_completeness = []
    for k in sorted(set(list(edge_counts.keys()) + edge_keys_expected)):
        cnt = edge_counts[k]
        pct = (cnt / n_edges) * 100 if n_edges else 0
        km = edge_length_by_key[k] / 1000.0
        edge_completeness.append({
            "target": "edge", "attribute": k, "present": cnt,
            "missing": n_edges - cnt, "pct_present": pct, "length_km": km, "total": n_edges
        })

    # Save completeness CSV
    comp_csv = RESULTS_DIR / f"{TAG}_attribute_completeness.csv"
    with open(comp_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "attribute", "present_count", "missing_count", "pct_present", "length_covered_km"])
        for r in node_completeness:
            w.writerow([r["target"], r["attribute"], r["present"], r["missing"], f"{r['pct_present']:.2f}%", "N/A"])
        for r in edge_completeness:
            w.writerow([r["target"], r["attribute"], r["present"], r["missing"], f"{r['pct_present']:.2f}%", f"{r['length_km']:.1f}"])
    log(f"Saved: {comp_csv.relative_to(PROJECT_ROOT)}")

    # ───────────────────────────────────────────────────────────────────────────
    # 2. Road Hierarchy Distribution
    # ───────────────────────────────────────────────────────────────────────────
    log("Analyzing road hierarchy & classifications ...")
    hw_counts = Counter()
    hw_lengths = Counter()

    for _, _, d in G.edges(data=True):
        hw = normalize_val(d.get("highway", "unclassified"))
        length_km = float(d.get("length", 0.0)) / 1000.0
        hw_counts[hw] += 1
        hw_lengths[hw] += length_km

    hierarchy_csv = RESULTS_DIR / f"{TAG}_road_hierarchy.csv"
    with open(hierarchy_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["highway_class", "edge_count", "pct_edges", "total_length_km", "pct_length"])
        for hw, cnt in hw_counts.most_common():
            km = hw_lengths[hw]
            pct_e = (cnt / n_edges) * 100 if n_edges else 0
            pct_l = (km / (total_network_len_m / 1000.0)) * 100 if total_network_len_m else 0
            w.writerow([hw, cnt, f"{pct_e:.2f}%", f"{km:.2f}", f"{pct_l:.2f}%"])
    log(f"Saved: {hierarchy_csv.relative_to(PROJECT_ROOT)}")

    # ───────────────────────────────────────────────────────────────────────────
    # 3. Node Features & Junction Degree
    # ───────────────────────────────────────────────────────────────────────────
    log("Analyzing node features & junction degrees ...")
    street_counts = Counter()
    node_hw_counts = Counter()
    node_junction_counts = Counter()

    for _, d in G.nodes(data=True):
        sc = d.get("street_count", 0)
        try:
            sc = int(sc)
        except (ValueError, TypeError):
            pass
        street_counts[sc] += 1

        if "highway" in d:
            node_hw_counts[normalize_val(d["highway"])] += 1
        if "junction" in d:
            node_junction_counts[normalize_val(d["junction"])] += 1

    node_feat_csv = RESULTS_DIR / f"{TAG}_node_features.csv"
    with open(node_feat_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["feature_category", "feature_name", "count", "pct_of_nodes"])
        for sc, cnt in sorted(street_counts.items(), key=lambda kv: -kv[1]):
            w.writerow(["street_count", f"degree_{sc}", cnt, f"{(cnt/n_nodes)*100:.2f}%"])
        for feat, cnt in node_hw_counts.most_common():
            w.writerow(["node_highway", feat, cnt, f"{(cnt/n_nodes)*100:.4f}%"])
        for junc, cnt in node_junction_counts.most_common():
            w.writerow(["node_junction", junc, cnt, f"{(cnt/n_nodes)*100:.4f}%"])
    log(f"Saved: {node_feat_csv.relative_to(PROJECT_ROOT)}")

    # ───────────────────────────────────────────────────────────────────────────
    # 4. Physical Road Properties
    # ───────────────────────────────────────────────────────────────────────────
    log("Summarizing physical properties (lanes, speed, oneway, bridges, tunnels) ...")
    oneway_counts = Counter()
    oneway_lengths = Counter()
    lanes_counts = Counter()
    lanes_lengths = Counter()
    maxspeed_counts = Counter()
    bridge_count = 0
    bridge_len_km = 0.0
    tunnel_count = 0
    tunnel_len_km = 0.0
    named_count = 0
    named_len_km = 0.0

    for _, _, d in G.edges(data=True):
        length_km = float(d.get("length", 0.0)) / 1000.0
        # oneway
        ow = d.get("oneway")
        is_ow = str(ow).lower() in ("true", "1", "yes")
        oneway_counts[is_ow] += 1
        oneway_lengths[is_ow] += length_km

        # name
        if "name" in d and d["name"]:
            named_count += 1
            named_len_km += length_km

        # lanes
        if "lanes" in d and d["lanes"]:
            l_val = normalize_val(d["lanes"])
            lanes_counts[str(l_val)] += 1
            lanes_lengths[str(l_val)] += length_km

        # speed
        if "maxspeed" in d and d["maxspeed"]:
            s_val = normalize_val(d["maxspeed"])
            maxspeed_counts[str(s_val)] += 1

        # bridge / tunnel
        if "bridge" in d and str(d["bridge"]).lower() not in ("no", "false", "0"):
            bridge_count += 1
            bridge_len_km += length_km
        if "tunnel" in d and str(d["tunnel"]).lower() not in ("no", "false", "0"):
            tunnel_count += 1
            tunnel_len_km += length_km

    road_props_csv = RESULTS_DIR / f"{TAG}_road_attributes_summary.csv"
    with open(road_props_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "category", "edge_count", "length_km", "pct_edges", "pct_length"])
        total_km = total_network_len_m / 1000.0
        w.writerow(["oneway", "one_way", oneway_counts[True], f"{oneway_lengths[True]:.2f}",
                    f"{(oneway_counts[True]/n_edges)*100:.2f}%", f"{(oneway_lengths[True]/total_km)*100:.2f}%"])
        w.writerow(["oneway", "bidirectional", oneway_counts[False], f"{oneway_lengths[False]:.2f}",
                    f"{(oneway_counts[False]/n_edges)*100:.2f}%", f"{(oneway_lengths[False]/total_km)*100:.2f}%"])
        w.writerow(["naming", "named", named_count, f"{named_len_km:.2f}",
                    f"{(named_count/n_edges)*100:.2f}%", f"{(named_len_km/total_km)*100:.2f}%"])
        w.writerow(["naming", "unnamed", n_edges - named_count, f"{total_km - named_len_km:.2f}",
                    f"{((n_edges - named_count)/n_edges)*100:.2f}%", f"{((total_km - named_len_km)/total_km)*100:.2f}%"])
        w.writerow(["structures", "bridge_or_flyover", bridge_count, f"{bridge_len_km:.2f}",
                    f"{(bridge_count/n_edges)*100:.2f}%", f"{(bridge_len_km/total_km)*100:.2f}%"])
        w.writerow(["structures", "tunnel_or_underpass", tunnel_count, f"{tunnel_len_km:.2f}",
                    f"{(tunnel_count/n_edges)*100:.2f}%", f"{(tunnel_len_km/total_km)*100:.2f}%"])
    log(f"Saved: {road_props_csv.relative_to(PROJECT_ROOT)}")

    # ───────────────────────────────────────────────────────────────────────────
    # Comprehensive Text Report
    # ───────────────────────────────────────────────────────────────────────────
    total_km = total_network_len_m / 1000.0
    report_path = RESULTS_DIR / f"{TAG}_data_profile_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"""================================================================================
OSM ROAD NETWORK DATA PROFILE: {CITY} ({TAG})
================================================================================
Total Nodes (Intersections & Dead Ends) : {n_nodes:,}
Total Directed Edges (Street Segments)  : {n_edges:,}
Total Street Length (Directed)          : {total_km:,.1f} km

1. ATTRIBUTE COMPLETENESS SUMMARY
--------------------------------------------------------------------------------
Node Attributes (n = {n_nodes:,}):
  • Coordinates (x, y)          : {node_counts['x']:,} / {n_nodes:,} (100.0%)
  • Physical Degree (street_cnt): {node_counts['street_count']:,} / {n_nodes:,} (100.0%)
  • Traffic Infrastructure      : {node_counts['highway']:,} / {n_nodes:,} ({(node_counts['highway']/n_nodes)*100:.2f}%)
  • Roundabout/Junction Tagged  : {node_counts['junction']:,} / {n_nodes:,} ({(node_counts['junction']/n_nodes)*100:.3f}%)

Edge Attributes (m = {n_edges:,}, total length = {total_km:,.1f} km):
  • Essential Topology (osmid, highway, oneway, reversed, length): 100.0% complete
  • Street Name                 : {edge_counts['name']:,} ({(edge_counts['name']/n_edges)*100:.1f}% edges, {named_len_km:,.1f} km)
    └─ Unnamed Streets           : {n_edges - edge_counts['name']:,} ({((n_edges - edge_counts['name'])/n_edges)*100:.1f}% edges, {total_km - named_len_km:,.1f} km)
  • Curved Centerline Geometry  : {edge_counts['geometry']:,} ({(edge_counts['geometry']/n_edges)*100:.1f}%) — remainder are straight links
  • Lane Count (lanes)          : {edge_counts['lanes']:,} ({(edge_counts['lanes']/n_edges)*100:.1f}%)
  • Speed Limit (maxspeed)      : {edge_counts['maxspeed']:,} ({(edge_counts['maxspeed']/n_edges)*100:.1f}%)
  • Access Control              : {edge_counts['access']:,} ({(edge_counts['access']/n_edges)*100:.1f}%)
  • Elevated Bridge / Flyover   : {bridge_count:,} ({(bridge_count/n_edges)*100:.2f}%, {bridge_len_km:.1f} km)
  • Highway Route Code (ref)    : {edge_counts['ref']:,} ({(edge_counts['ref']/n_edges)*100:.2f}%)
  • Width                       : {edge_counts['width']:,} ({(edge_counts['width']/n_edges)*100:.2f}%)
  • Tunnel / Underpass          : {tunnel_count:,} ({(tunnel_count/n_edges)*100:.3f}%, {tunnel_len_km:.1f} km)

2. ROAD HIERARCHY BREAKDOWN
--------------------------------------------------------------------------------
""")
        for hw, cnt in hw_counts.most_common():
            km = hw_lengths[hw]
            f.write(f"  • {hw:22s}: {cnt:7,d} edges ({cnt/n_edges:5.1%}) | {km:8.1f} km ({km/total_km:5.1%})\n")

        f.write(f"""
3. JUNCTION DEGREES & TRAFFIC SIGNALS
--------------------------------------------------------------------------------
Intersection Types (by incident street count):
  • 3-way T-junctions (deg 3)   : {street_counts[3]:,d} ({(street_counts[3]/n_nodes)*100:.1f}%)
  • Dead ends / Cul-de-sacs (1) : {street_counts[1]:,d} ({(street_counts[1]/n_nodes)*100:.1f}%)
  • 4-way Crossroads (deg 4)    : {street_counts[4]:,d} ({(street_counts[4]/n_nodes)*100:.1f}%)
  • 5-way+ Complex junctions    : {street_counts[5] + street_counts[6]:,d} ({((street_counts[5] + street_counts[6])/n_nodes)*100:.2f}%)
  • Interstitial / Simpl. (2)   : {street_counts[2]:,d} ({(street_counts[2]/n_nodes)*100:.2f}%)

Specialized Intersection Infrastructure:
""")
        for feat, cnt in node_hw_counts.most_common():
            f.write(f"  • {feat:26s}: {cnt:5,d} nodes\n")
        for junc, cnt in node_junction_counts.most_common():
            f.write(f"  • junction={junc:17s}: {cnt:5,d} nodes\n")

        f.write(f"""
4. TRAFFIC FLOW & PHYSICAL CHARACTERISTICS
--------------------------------------------------------------------------------
  • One-Way Street Segments     : {oneway_counts[True]:,d} ({oneway_counts[True]/n_edges:.1%}) | {oneway_lengths[True]:,.1f} km ({oneway_lengths[True]/total_km:.1%})
  • Two-Way Street Segments     : {oneway_counts[False]:,d} ({oneway_counts[False]/n_edges:.1%}) | {oneway_lengths[False]:,.1f} km ({oneway_lengths[False]/total_km:.1%})
  • Bridges / Flyovers / Viaducts: {bridge_count:,d} ({bridge_len_km:.1f} km)
  • Tunnels / Underpasses       : {tunnel_count:,d} ({tunnel_len_km:.1f} km)
  • Named Streets               : {named_count:,d} ({named_len_km:,.1f} km)
  • Unnamed (Local/Gully) Streets: {n_edges - named_count:,d} ({total_km - named_len_km:,.1f} km)
""")
    log(f"Saved report: {report_path.relative_to(PROJECT_ROOT)}")

    # ───────────────────────────────────────────────────────────────────────────
    # Figures Generation
    # ───────────────────────────────────────────────────────────────────────────
    log("Rendering publication figures in figures/ ...")

    # Figure 1: Attribute Completeness (Known vs Missing)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5), facecolor=SURFACE, constrained_layout=True)
    
    # Sort edges by completeness
    edge_sorted = sorted(edge_completeness, key=lambda x: x["pct_present"])
    e_names = [e["attribute"] for e in edge_sorted]
    e_pcts = [e["pct_present"] for e in edge_sorted]
    e_colors = [GREEN if p >= 90 else AMBER if p >= 20 else RED for p in e_pcts]

    ax1.barh(e_names, e_pcts, color=e_colors, height=0.65, edgecolor=SURFACE)
    style_ax(ax1, title="Edge Attribute Completeness (%)", xlabel="Percentage Present (%)")
    ax1.set_xlim(0, 108)
    for i, v in enumerate(e_pcts):
        ax1.text(v + 1.2, i, f"{v:.1f}%", va="center", color=INKD, fontsize=8, fontweight="bold")

    # Nodes
    node_sorted = sorted(node_completeness, key=lambda x: x["pct_present"])
    n_names = [n["attribute"] for n in node_sorted]
    n_pcts = [n["pct_present"] for n in node_sorted]
    n_colors = [GREEN if p >= 90 else AMBER if p >= 20 else RED for p in n_pcts]

    ax2.barh(n_names, n_pcts, color=n_colors, height=0.5, edgecolor=SURFACE)
    style_ax(ax2, title="Node Attribute Completeness (%)", xlabel="Percentage Present (%)")
    ax2.set_xlim(0, 108)
    for i, v in enumerate(n_pcts):
        ax2.text(v + 1.2, i, f"{v:.2f}%", va="center", color=INKD, fontsize=8, fontweight="bold")

    # Legend for tier colors
    from matplotlib.patches import Patch
    tier_legend = [
        Patch(facecolor=GREEN, label="High (≥90%)"),
        Patch(facecolor=AMBER, label="Moderate (20-89%)"),
        Patch(facecolor=RED, label="Sparse (<20%)"),
    ]
    ax1.legend(handles=tier_legend, loc="lower right", facecolor=SURFACE, edgecolor=SPINE,
               labelcolor=INKD, fontsize=9)

    fig.suptitle(f"OpenStreetMap Data Completeness — {CITY} Road Network", color=INKD, fontsize=15, fontweight="bold")
    fig1_path = FIGURES_DIR / f"{TAG}_attribute_completeness_dark.png"
    fig.savefig(fig1_path, dpi=180, facecolor=SURFACE)
    plt.close(fig)
    log(f"Saved: {fig1_path.relative_to(PROJECT_ROOT)}")

    # Figure 2: Road Hierarchy (Counts & Length km)
    top_hw = hw_counts.most_common(10)
    hw_labels = [k for k, _ in top_hw]
    hw_edge_cnt = [cnt for _, cnt in top_hw]
    hw_edge_len = [hw_lengths[k] for k in hw_labels]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5), facecolor=SURFACE, constrained_layout=True)
    bars1 = ax1.bar(hw_labels, hw_edge_cnt, color=BLUE, edgecolor=SURFACE, width=0.65)
    style_ax(ax1, title="Street Segments by Classification", ylabel="Count of Directed Segments (log scale)")
    ax1.set_yscale("log")
    ax1.set_xticks(range(len(hw_labels)))
    ax1.set_xticklabels(hw_labels, rotation=40, ha="right", color=INKD, fontsize=9)
    for bar in bars1:
        yval = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, yval * 1.15, f"{yval:,}", ha="center", va="bottom",
                 color=INKD, fontsize=7.5, rotation=0)

    bars2 = ax2.bar(hw_labels, hw_edge_len, color=ORANGE, edgecolor=SURFACE, width=0.65)
    style_ax(ax2, title="Total Road Length by Classification", ylabel="Length (km, log scale)")
    ax2.set_yscale("log")
    ax2.set_xticks(range(len(hw_labels)))
    ax2.set_xticklabels(hw_labels, rotation=40, ha="right", color=INKD, fontsize=9)
    for bar in bars2:
        yval = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, yval * 1.15, f"{yval:.0f} km", ha="center", va="bottom",
                 color=INKD, fontsize=7.5, rotation=0)

    fig.suptitle(f"Road Network Hierarchy Breakdown — {CITY}", color=INKD, fontsize=15, fontweight="bold")
    fig2_path = FIGURES_DIR / f"{TAG}_road_hierarchy_distribution_dark.png"
    fig.savefig(fig2_path, dpi=180, facecolor=SURFACE)
    plt.close(fig)
    log(f"Saved: {fig2_path.relative_to(PROJECT_ROOT)}")

    # Figure 3: Junction Types & Traffic Signals
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), facecolor=SURFACE, constrained_layout=True)
    
    # Degrees
    deg_labels = ["3-way (T-junction)", "1-way (Dead end)", "4-way (Crossroad)", "2-way (Loop/Link)", "5-way+", "6-way"]
    deg_keys = [3, 1, 4, 2, 5, 6]
    deg_vals = [street_counts.get(k, 0) for k in deg_keys]
    deg_colors = [BLUE, CYAN, AMBER, PURPLE, RED, "#e91e63"]

    bars = ax1.bar(deg_labels, deg_vals, color=deg_colors, width=0.6, edgecolor=SURFACE)
    style_ax(ax1, title="Intersection Topologies (Physical Degree)", ylabel="Number of Intersections")
    ax1.set_xticks(range(len(deg_labels)))
    ax1.set_xticklabels(deg_labels, rotation=30, ha="right", color=INKD, fontsize=8.5)
    for bar in bars:
        h = bar.get_height()
        pct = (h / n_nodes) * 100
        ax1.text(bar.get_x() + bar.get_width()/2, h + 1500, f"{h:,}\n({pct:.1f}%)", ha="center", va="bottom",
                 color=INKD, fontsize=8)
    ax1.set_ylim(0, max(deg_vals) * 1.15)

    # Traffic signals & features
    sig_labels = [k for k, _ in node_hw_counts.most_common(7)]
    sig_vals = [node_hw_counts[k] for k in sig_labels]
    ax2.barh(sig_labels, sig_vals, color=GREEN, edgecolor=SURFACE, height=0.6)
    style_ax(ax2, title="Specialized Traffic Control Nodes", xlabel="Node Count")
    for i, v in enumerate(sig_vals):
        ax2.text(v + 3, i, f"{v}", va="center", color=INKD, fontsize=9, fontweight="bold")
    ax2.set_xlim(0, max(sig_vals) * 1.18 if sig_vals else 10)

    fig.suptitle(f"Intersection Topologies & Traffic Control Infrastructure — {CITY}",
                 color=INKD, fontsize=15, fontweight="bold")
    fig3_path = FIGURES_DIR / f"{TAG}_node_features_and_junctions_dark.png"
    fig.savefig(fig3_path, dpi=180, facecolor=SURFACE)
    plt.close(fig)
    log(f"Saved: {fig3_path.relative_to(PROJECT_ROOT)}")

    # Figure 4: Physical Road Properties (4-panel summary)
    fig, axs = plt.subplots(2, 2, figsize=(14, 11), facecolor=SURFACE, constrained_layout=True)
    
    # 4a: Oneway vs Bidirectional (Length share)
    ax = axs[0, 0]
    ax.set_facecolor(SURFACE)
    labels = ["Bidirectional", "One-Way"]
    sizes = [oneway_lengths[False], oneway_lengths[True]]
    colors = [BLUE, ORANGE]
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, autopct="%1.1f%%", colors=colors, startangle=140,
        textprops=dict(color=INKD, fontsize=10),
        wedgeprops=dict(width=0.45, edgecolor=SURFACE, linewidth=2)
    )
    for at in autotexts:
        at.set_color(INKD)
        at.set_fontweight("bold")
    ax.set_title(f"Street Directionality (% of {total_km:,.0f} km)", color=INKD, fontsize=12, fontweight="bold")

    # 4b: Named vs Unnamed
    ax = axs[0, 1]
    ax.set_facecolor(SURFACE)
    labels = ["Unnamed Street", "Named Street"]
    sizes = [total_km - named_len_km, named_len_km]
    colors = ["#55554f", GREEN]
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, autopct="%1.1f%%", colors=colors, startangle=140,
        textprops=dict(color=INKD, fontsize=10),
        wedgeprops=dict(width=0.45, edgecolor=SURFACE, linewidth=2)
    )
    for at in autotexts:
        at.set_color(INKD)
        at.set_fontweight("bold")
    ax.set_title("Street Naming Coverage (% of Length)", color=INKD, fontsize=12, fontweight="bold")

    # 4c: Lanes Distribution
    ax = axs[1, 0]
    top_lanes = lanes_counts.most_common(5)
    l_keys = [str(k) for k, _ in top_lanes]
    l_vals = [cnt for _, cnt in top_lanes]
    bars = ax.bar(l_keys, l_vals, color=CYAN, edgecolor=SURFACE, width=0.55)
    style_ax(ax, title=f"Reported Lane Counts (n = {sum(lanes_counts.values()):,} edges with data)",
             xlabel="Number of Lanes", ylabel="Segment Count")
    for b in bars:
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 300, f"{int(b.get_height()):,}",
                ha="center", va="bottom", color=INKD, fontsize=8)

    # 4d: Grade Separation: Bridges & Tunnels
    ax = axs[1, 1]
    struct_names = ["Elevated Flyovers\n& Bridges", "Tunnels &\nUnderpasses"]
    struct_lens = [bridge_len_km, tunnel_len_km]
    struct_counts = [bridge_count, tunnel_count]
    bars = ax.bar(struct_names, struct_lens, color=[AMBER, RED], width=0.45, edgecolor=SURFACE)
    style_ax(ax, title="Grade-Separated Infrastructure", ylabel="Total Length (km)")
    for i, b in enumerate(bars):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 1.5,
                f"{b.get_height():.1f} km\n({struct_counts[i]:,} segments)",
                ha="center", va="bottom", color=INKD, fontsize=9, fontweight="bold")
    ax.set_ylim(0, max(struct_lens) * 1.25 if struct_lens else 10)

    fig.suptitle(f"Physical Street Network Characteristics — {CITY}", color=INKD, fontsize=15, fontweight="bold")
    fig4_path = FIGURES_DIR / f"{TAG}_road_attributes_profile_dark.png"
    fig.savefig(fig4_path, dpi=180, facecolor=SURFACE)
    plt.close(fig)
    log(f"Saved: {fig4_path.relative_to(PROJECT_ROOT)}")

    log("Data analysis pipeline complete!")


if __name__ == "__main__":
    main()
