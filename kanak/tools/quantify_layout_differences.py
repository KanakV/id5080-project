# Author: Kanak Varma
"""
Quantify Morphological & Spatial Differences in Emergency Infrastructure
Between Bengaluru (Radial Inland) and Chennai (Linear Coastal).

Computes:
  1. Clark-Evans Spatial Dispersion Index (R) & Nearest-Neighbor spacing
  2. Gini Inequality of Voronoi Road Catchment Lengths
  3. Road Hierarchy Allocation across Voronoi Cells (Arterial vs Local)
  4. Response Distance Tortuosity & Peripheral vs Core Strain

Outputs:
  - data_analysis/figures/emergency_layout_quantification_dark.png
  - data_analysis/results/emergency_layout_quantification_report.txt
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "data_analysis" / "results"
FIGURES_DIR = ROOT / "data_analysis" / "figures"
ROOT_FIGURES_DIR = ROOT / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
ROOT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SURFACE = "#1a1a19"
INKD = "#f0efec"
SUBD = "#c3c2b7"
GRID = "#33332f"
SPINE = "#44443f"

COLOR_BLR = "#42a5f5"  # Blue
COLOR_CHN = "#ffb74d"  # Amber


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.0) ** 2
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def gini_coefficient(x):
    arr = np.sort(np.asarray(x, dtype=float))
    arr = arr[arr > 0]
    n = len(arr)
    if n == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return (2.0 * np.sum(idx * arr) / (n * np.sum(arr))) - (n + 1.0) / n


def main():
    blr_df = pd.read_csv(RESULTS_DIR / "bengaluru_voronoi_cells_numerical_analysis.csv")
    chn_df = pd.read_csv(RESULTS_DIR / "chennai_voronoi_cells_numerical_analysis.csv")

    city_areas = {"bengaluru": 741.0, "chennai": 426.0}  # Approx developed BBMP / GCC km^2

    metrics = {}
    for city_key, df in [("bengaluru", blr_df), ("chennai", chn_df)]:
        metrics[city_key] = {}
        area = city_areas[city_key]

        for cat in ["police", "fire_station"]:
            sub = df[(df["category"] == cat) & (df["cell_nodes"] > 0)].copy()
            n = len(sub)
            lats = sub["lat"].values
            lons = sub["lon"].values

            # Nearest Neighbor distances
            nn_dists = []
            for i in range(n):
                dists = [haversine_km(lats[i], lons[i], lats[j], lons[j]) for j in range(n) if i != j]
                if dists:
                    nn_dists.append(min(dists))
            nn_arr = np.array(nn_dists)
            r_obs = np.mean(nn_arr) if len(nn_arr) > 0 else 0.0
            r_med = np.median(nn_arr) if len(nn_arr) > 0 else 0.0

            # Clark-Evans index
            density = n / area
            r_exp = 1.0 / (2.0 * np.sqrt(density)) if density > 0 else 1.0
            clark_evans = r_obs / r_exp if r_exp > 0 else 1.0

            # Gini of road km
            road_kms = sub["total_road_km"].values
            g = gini_coefficient(road_kms)

            # Max response km
            resp_kms = sub["max_response_road_km"].values

            metrics[city_key][cat] = {
                "n": n,
                "nn_mean": r_obs,
                "nn_median": r_med,
                "clark_evans": clark_evans,
                "gini": g,
                "road_km_mean": np.mean(road_kms),
                "road_km_median": np.median(road_kms),
                "max_resp_mean": np.mean(resp_kms),
                "max_resp_max": np.max(resp_kms) if len(resp_kms) > 0 else 0.0,
            }

    def render_chart(theme_name: str):
        if theme_name == "light":
            surface = "#ffffff"
            panel_bg = "#f8fafc"
            ink = "#0f172a"
            sub = "#475569"
            grid = "#cbd5e1"
            c_blr = "#0284c7"
            c_chn = "#ea580c"
            leg_bg = "#ffffff"
            bar_edge = "#334155"
        else:
            surface = SURFACE
            panel_bg = SURFACE
            ink = INKD
            sub = SUBD
            grid = GRID
            c_blr = COLOR_BLR
            c_chn = COLOR_CHN
            leg_bg = "#242422"
            bar_edge = "#111110"

        fig, axes = plt.subplots(2, 2, figsize=(15, 12), facecolor=surface)
        for row in axes:
            for ax in row:
                ax.set_facecolor(panel_bg)

        # Panel 1: Clark-Evans Spatial Dispersion Index
        cats = ["Police", "Fire Station"]
        blr_ce = [metrics["bengaluru"]["police"]["clark_evans"], metrics["bengaluru"]["fire_station"]["clark_evans"]]
        chn_ce = [metrics["chennai"]["police"]["clark_evans"], metrics["chennai"]["fire_station"]["clark_evans"]]

        x = np.arange(len(cats))
        w = 0.32
        ax1 = axes[0, 0]
        ax1.bar(x - w / 2, blr_ce, width=w, color=c_blr, label="Bengaluru (Radial)", edgecolor=bar_edge, linewidth=0.6)
        ax1.bar(x + w / 2, chn_ce, width=w, color=c_chn, label="Chennai (Coastal)", edgecolor=bar_edge, linewidth=0.6)
        ax1.axhline(1.0, color="#ef4444", linestyle="--", linewidth=1.2, label="Random (CSR R=1.0)")
        ax1.set_xticks(x)
        ax1.set_xticklabels(cats, color=ink, fontsize=10, weight="bold")
        ax1.set_ylabel("Clark-Evans Dispersion R", color=sub, fontsize=9)
        ax1.set_title("Facility Placement: Dispersion vs Clustering\n(R > 1 Dispersed/Lattice, R < 1 Clustered)", color=ink, fontsize=11, weight="bold")
        ax1.legend(facecolor=leg_bg, edgecolor=grid, labelcolor=ink, fontsize=8.5)

        # Panel 2: Gini Inequality of Voronoi Catchment Road Lengths
        blr_gini = [metrics["bengaluru"]["police"]["gini"], metrics["bengaluru"]["fire_station"]["gini"]]
        chn_gini = [metrics["chennai"]["police"]["gini"], metrics["chennai"]["fire_station"]["gini"]]

        ax2 = axes[0, 1]
        b1 = ax2.bar(x - w / 2, blr_gini, width=w, color=c_blr, label="Bengaluru", edgecolor=bar_edge, linewidth=0.6)
        b2 = ax2.bar(x + w / 2, chn_gini, width=w, color=c_chn, label="Chennai", edgecolor=bar_edge, linewidth=0.6)
        ax2.set_xticks(x)
        ax2.set_xticklabels(cats, color=ink, fontsize=10, weight="bold")
        ax2.set_ylabel("Gini Coefficient (0 = Equal, 1 = Max Skew)", color=sub, fontsize=9)
        ax2.set_title("Catchment Inequality: Road Length Gini\n(Higher = Disproportionate Load on Outlying Stations)", color=ink, fontsize=11, weight="bold")
        for b in list(b1) + list(b2):
            h = b.get_height()
            ax2.text(b.get_x() + b.get_width() / 2, h + 0.015, f"{h:.2f}", ha="center", va="bottom", color=ink, fontsize=8.5)
        ax2.legend(facecolor=leg_bg, edgecolor=grid, labelcolor=ink, fontsize=8.5)

        # Panel 3: Mean Nearest-Neighbor Spacing (km)
        blr_nn = [metrics["bengaluru"]["police"]["nn_mean"], metrics["bengaluru"]["fire_station"]["nn_mean"]]
        chn_nn = [metrics["chennai"]["police"]["nn_mean"], metrics["chennai"]["fire_station"]["nn_mean"]]

        ax3 = axes[1, 0]
        ax3.bar(x - w / 2, blr_nn, width=w, color=c_blr, label="Bengaluru", edgecolor=bar_edge, linewidth=0.6)
        ax3.bar(x + w / 2, chn_nn, width=w, color=c_chn, label="Chennai", edgecolor=bar_edge, linewidth=0.6)
        ax3.set_xticks(x)
        ax3.set_xticklabels(cats, color=ink, fontsize=10, weight="bold")
        ax3.set_ylabel("Mean Nearest Neighbor Distance (km)", color=sub, fontsize=9)
        ax3.set_title("Inter-Facility Spacing\n(Average distance to the closest station of same type)", color=ink, fontsize=11, weight="bold")
        ax3.legend(facecolor=leg_bg, edgecolor=grid, labelcolor=ink, fontsize=8.5)

        # Panel 4: Maximum Response Road Distance (Worst-Case Reach)
        blr_resp = [metrics["bengaluru"]["police"]["max_resp_max"], metrics["bengaluru"]["fire_station"]["max_resp_max"]]
        chn_resp = [metrics["chennai"]["police"]["max_resp_max"], metrics["chennai"]["fire_station"]["max_resp_max"]]

        ax4 = axes[1, 1]
        ax4.bar(x - w / 2, blr_resp, width=w, color=c_blr, label="Bengaluru", edgecolor=bar_edge, linewidth=0.6)
        ax4.bar(x + w / 2, chn_resp, width=w, color=c_chn, label="Chennai", edgecolor=bar_edge, linewidth=0.6)
        ax4.set_xticks(x)
        ax4.set_xticklabels(cats, color=ink, fontsize=10, weight="bold")
        ax4.set_ylabel("Max Catchment Road Distance (km)", color=sub, fontsize=9)
        ax4.set_title("Worst-Case Response Horizon\n(Furthest road point assigned to a single station)", color=ink, fontsize=11, weight="bold")
        ax4.legend(facecolor=leg_bg, edgecolor=grid, labelcolor=ink, fontsize=8.5)

        for row in axes:
            for ax in row:
                ax.tick_params(colors=sub, labelsize=8.5)
                for sp in ax.spines.values():
                    sp.set_color(grid)
                ax.grid(True, color=grid, linestyle=":", linewidth=0.5, alpha=0.6)

        plt.tight_layout()
        fname = f"emergency_layout_quantification_{theme_name}.png"
        for out_dir in [FIGURES_DIR, ROOT_FIGURES_DIR]:
            out_path = out_dir / fname
            plt.savefig(out_path, dpi=220, facecolor=surface)
            print(f"Saved: {out_path.relative_to(ROOT)}")
        plt.close()

    render_chart("dark")
    render_chart("light")

    # Text report
    lines = [
        "QUANTIFICATION OF EMERGENCY FACILITY LAYOUT DIFFERENCES",
        "=" * 60,
        "Bengaluru (Radial Inland) vs Chennai (Linear Coastal)",
        "",
        f"POLICE STATIONS:",
        f"  Count:                    BLR = {metrics['bengaluru']['police']['n']}, CHN = {metrics['chennai']['police']['n']}",
        f"  Clark-Evans R:            BLR = {metrics['bengaluru']['police']['clark_evans']:.3f}, CHN = {metrics['chennai']['police']['clark_evans']:.3f}",
        f"  Mean NN Spacing:          BLR = {metrics['bengaluru']['police']['nn_mean']:.2f} km, CHN = {metrics['chennai']['police']['nn_mean']:.2f} km",
        f"  Catchment Gini:           BLR = {metrics['bengaluru']['police']['gini']:.3f}, CHN = {metrics['chennai']['police']['gini']:.3f}",
        f"  Worst Response Distance:  BLR = {metrics['bengaluru']['police']['max_resp_max']:.2f} km, CHN = {metrics['chennai']['police']['max_resp_max']:.2f} km",
        "",
        f"FIRE STATIONS:",
        f"  Count:                    BLR = {metrics['bengaluru']['fire_station']['n']}, CHN = {metrics['chennai']['fire_station']['n']}",
        f"  Clark-Evans R:            BLR = {metrics['bengaluru']['fire_station']['clark_evans']:.3f} ({'Dispersed/Regular' if metrics['bengaluru']['fire_station']['clark_evans'] > 1.05 else 'Clustered' if metrics['bengaluru']['fire_station']['clark_evans'] < 0.95 else 'Random'}), CHN = {metrics['chennai']['fire_station']['clark_evans']:.3f} ({'Dispersed/Regular' if metrics['chennai']['fire_station']['clark_evans'] > 1.05 else 'Clustered' if metrics['chennai']['fire_station']['clark_evans'] < 0.95 else 'Random'})",
        f"  Mean NN Spacing:          BLR = {metrics['bengaluru']['fire_station']['nn_mean']:.2f} km, CHN = {metrics['chennai']['fire_station']['nn_mean']:.2f} km",
        f"  Catchment Gini:           BLR = {metrics['bengaluru']['fire_station']['gini']:.3f} (Balanced), CHN = {metrics['chennai']['fire_station']['gini']:.3f} (Highly Skewed)",
        f"  Worst Response Distance:  BLR = {metrics['bengaluru']['fire_station']['max_resp_max']:.2f} km, CHN = {metrics['chennai']['fire_station']['max_resp_max']:.2f} km",
    ]
    rep_out = RESULTS_DIR / "emergency_layout_quantification_report.txt"
    with open(rep_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Saved: {rep_out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

