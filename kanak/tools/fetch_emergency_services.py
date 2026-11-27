# Author: Kanak Varma
"""
Fetch and export emergency services (Police, Hospitals, Fire Stations) for
Bengaluru and Chennai via OpenStreetMap Overpass API.

Amenities fetched:
  - Police: amenity=police
  - Hospitals: amenity=hospital, amenity=clinic (with emergency=yes or designated hospital)
  - Fire Stations: amenity=fire_station

Outputs:
  - data/bengaluru_emergency_services.csv
  - data/chennai_emergency_services.csv
  - data_analysis/results/emergency_services_report.txt
  - data_analysis/figures/bengaluru_emergency_services_dark.png
  - data_analysis/figures/chennai_emergency_services_dark.png (if graphml available)

Usage:
  python tools/fetch_emergency_services.py
  python tools/fetch_emergency_services.py bengaluru
  python tools/fetch_emergency_services.py chennai
"""

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import osmnx as ox
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

# ── Paths & Setup ─────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "data_analysis" / "results"
FIGURES_DIR = ROOT / "data_analysis" / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

try:
    from tools.curated_fire_stations import CHENNAI_FIRE_STATIONS, BENGALURU_FIRE_STATIONS
except ImportError:
    from curated_fire_stations import CHENNAI_FIRE_STATIONS, BENGALURU_FIRE_STATIONS

# Dark theme palette
SURFACE = "#1a1a19"
INKD = "#f0efec"
SUBD = "#c3c2b7"
GRID = "#33332f"
ROAD_COLOR = "#2c2c29"

# Facility Colors
COLOR_POLICE = "#29b6f6"       # Bright Cyan / Blue
COLOR_HOSPITAL = "#ef5350"     # Crimson / Red
COLOR_FIRE = "#ffa726"         # Vivid Orange / Amber

OVERPASS_URLS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
    "https://maps.mail.ru/osm/tools/overpass/api/interpreter",
]

# Bounding boxes (min_lat, min_lon, max_lat, max_lon)
CITY_BBOXES = {
    "bengaluru": {
        "bbox": (12.78, 77.42, 13.18, 77.78),
        "graphml": DATA_DIR / "bengaluru_city_bengaluru_drive.graphml",
        "title": "Bengaluru",
    },
    "chennai": {
        "bbox": (12.85, 80.05, 13.25, 80.35),
        "graphml": DATA_DIR / "chennai_city_chennai_drive.graphml",
        "title": "Chennai",
    },
}


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def query_overpass(query_str: str) -> dict:
    data = urllib.parse.urlencode({"data": query_str}).encode("utf-8")
    for endpoint in OVERPASS_URLS:
        try:
            log(f"Querying {endpoint} ...")
            req = urllib.request.Request(
                endpoint,
                data=data,
                headers={"User-Agent": "AntigravityEmergencyAnalysis/1.0"},
            )
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            log(f"Endpoint {endpoint} failed ({e}), trying next ...")
    raise RuntimeError("All Overpass endpoints failed.")


def build_overpass_query(bbox):
    s, w, n, e = bbox
    return f"""
    [out:json][timeout:90];
    (
      node["amenity"="police"]({s}, {w}, {n}, {e});
      way["amenity"="police"]({s}, {w}, {n}, {e});
      relation["amenity"="police"]({s}, {w}, {n}, {e});

      node["amenity"="hospital"]({s}, {w}, {n}, {e});
      way["amenity"="hospital"]({s}, {w}, {n}, {e});
      relation["amenity"="hospital"]({s}, {w}, {n}, {e});

      node["amenity"="fire_station"]({s}, {w}, {n}, {e});
      way["amenity"="fire_station"]({s}, {w}, {n}, {e});
      relation["amenity"="fire_station"]({s}, {w}, {n}, {e});
    );
    out center tags;
    """


def fetch_city_amenities(city_key: str):
    cfg = CITY_BBOXES[city_key]
    bbox = cfg["bbox"]
    cache_file = DATA_DIR / f"{city_key}_emergency_raw.json"

    if cache_file.exists():
        log(f"Loading cached {city_key} raw data from {cache_file.name} ...")
        with open(cache_file, "r", encoding="utf-8") as f:
            data = json.load(f)
    else:
        log(f"Fetching {city_key} emergency facilities from Overpass...")
        query = build_overpass_query(bbox)
        data = query_overpass(query)
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(data, f)
        log(f"Cached raw response: {cache_file.relative_to(ROOT)}")

    elements = data.get("elements", [])
    log(f"Parsed {len(elements):,} raw elements for {city_key}.")

    records = []
    seen = set()

    for el in elements:
        tags = el.get("tags", {})
        amenity = tags.get("amenity")
        if amenity not in ("police", "hospital", "fire_station"):
            continue

        # Get coordinates (node has lat/lon, way/relation with 'out center' has center dict)
        if "lat" in el and "lon" in el:
            lat, lon = el["lat"], el["lon"]
        elif "center" in el:
            lat, lon = el["center"]["lat"], el["center"]["lon"]
        else:
            continue

        name = tags.get("name") or tags.get("operator") or f"Unnamed {amenity.replace('_', ' ').title()}"
        emergency = tags.get("emergency", "yes" if amenity in ("fire_station", "police") else tags.get("emergency", "unknown"))
        phone = tags.get("phone") or tags.get("contact:phone") or ""
        addr = tags.get("addr:street") or tags.get("addr:suburb") or tags.get("addr:full") or ""

        # Dedup key
        key = (amenity, round(lat, 4), round(lon, 4))
        if key in seen:
            continue
        seen.add(key)

        records.append({
            "city": city_key,
            "osm_id": el.get("id"),
            "osm_type": el.get("type"),
            "category": amenity,
            "name": name,
            "lat": lat,
            "lon": lon,
            "emergency": emergency,
            "phone": phone,
            "address": addr,
        })

    # Integrate official curated fire stations
    curated_map = {"bengaluru": BENGALURU_FIRE_STATIONS, "chennai": CHENNAI_FIRE_STATIONS}
    if city_key in curated_map:
        records = [r for r in records if r["category"] != "fire_station"]
        for idx, fs in enumerate(curated_map[city_key]):
            records.append({
                "city": city_key,
                "osm_id": f"curated_fs_{idx+1:02d}",
                "osm_type": "node",
                "category": "fire_station",
                "name": fs["name"],
                "lat": fs["lat"],
                "lon": fs["lon"],
                "emergency": "yes",
                "phone": "101",
                "address": fs.get("address", ""),
            })
        log(f"Injected {len(curated_map[city_key])} official curated fire stations for {city_key}.")

    # Save CSV
    csv_path = DATA_DIR / f"{city_key}_emergency_services.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["city", "osm_id", "osm_type", "category", "name", "lat", "lon", "emergency", "phone", "address"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    log(f"Wrote {len(records)} clean records to {csv_path.relative_to(ROOT)}")
    return records


def render_map(city_key: str, records: list):
    cfg = CITY_BBOXES[city_key]
    graph_path = cfg["graphml"]
    fig_path = FIGURES_DIR / f"{city_key}_emergency_services_dark.png"

    log(f"Rendering emergency services map for {cfg['title']} ...")
    fig, ax = plt.subplots(figsize=(11, 11), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    has_graph = graph_path.exists()
    if has_graph:
        log(f"Loading road network from {graph_path.name} for background...")
        G = ox.load_graphml(str(graph_path))
        lines = []
        for u, v, data in G.edges(data=True):
            if "geometry" in data:
                xs, ys = data["geometry"].xy
                lines.append(list(zip(xs, ys)))
            else:
                lines.append([
                    (G.nodes[u]["x"], G.nodes[u]["y"]),
                    (G.nodes[v]["x"], G.nodes[v]["y"])
                ])
        lc = LineCollection(lines, colors=ROAD_COLOR, linewidths=0.5, alpha=0.5, zorder=1)
        ax.add_collection(lc)

    # Plot facility categories
    cats = {
        "police": {"color": COLOR_POLICE, "label": "Police Station", "size": 32, "zorder": 4, "marker": "o"},
        "hospital": {"color": COLOR_HOSPITAL, "label": "Hospital", "size": 36, "zorder": 5, "marker": "s"},
        "fire_station": {"color": COLOR_FIRE, "label": "Fire Station", "size": 55, "zorder": 6, "marker": "^"},
    }

    counts = Counter(r["category"] for r in records)

    for cat_name, style in cats.items():
        sub = [r for r in records if r["category"] == cat_name]
        if not sub:
            continue
        xs = [r["lon"] for r in sub]
        ys = [r["lat"] for r in sub]
        ax.scatter(
            xs, ys,
            s=style["size"],
            c=style["color"],
            marker=style["marker"],
            edgecolors="#111110",
            linewidth=0.6,
            alpha=0.9,
            zorder=style["zorder"],
            label=f"{style['label']} ({len(sub)})",
        )

    s, w, n, e = cfg["bbox"]
    if has_graph:
        xs = [data["x"] for _, data in G.nodes(data=True) if "x" in data]
        ys = [data["y"] for _, data in G.nodes(data=True) if "y" in data]
        if xs and ys:
            pad_x = (max(xs) - min(xs)) * 0.02
            pad_y = (max(ys) - min(ys)) * 0.02
            ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
            ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)
        else:
            ax.set_xlim(w, e)
            ax.set_ylim(s, n)
    else:
        ax.set_xlim(w, e)
        ax.set_ylim(s, n)

    ax.set_aspect("equal")
    ax.set_title(
        f"{cfg['title']} — Emergency Services Infrastructure (OSM)\n"
        f"Police: {counts['police']}  |  Hospitals: {counts['hospital']}  |  Fire Stations: {counts['fire_station']}",
        fontsize=13,
        color=INKD,
        weight="bold",
        pad=16,
    )

    leg = ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="#242422",
        edgecolor=GRID,
        labelcolor=INKD,
        fontsize=10,
    )
    for text in leg.get_texts():
        text.set_color(INKD)

    ax.tick_params(colors=SUBD, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.grid(True, color=GRID, linestyle="--", linewidth=0.4, alpha=0.5)
    ax.set_xlabel("Longitude", color=SUBD, fontsize=9)
    ax.set_ylabel("Latitude", color=SUBD, fontsize=9)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=250, facecolor=SURFACE)
    plt.close()
    log(f"Saved figure: {fig_path.relative_to(ROOT)}")


def main():
    target = sys.argv[1].lower() if len(sys.argv) > 1 else "all"
    cities = [target] if target in CITY_BBOXES else ["bengaluru", "chennai"]

    summary_lines = ["EMERGENCY INFRASTRUCTURE DATA SUMMARY (OpenStreetMap Overpass API)", "=" * 65, ""]

    for city in cities:
        records = fetch_city_amenities(city)
        counts = Counter(r["category"] for r in records)
        summary_lines.append(f"City: {city.upper()}")
        summary_lines.append(f"  Total Facilities Tagged: {len(records)}")
        summary_lines.append(f"  - Police Stations:     {counts.get('police', 0)}")
        summary_lines.append(f"  - Hospitals:           {counts.get('hospital', 0)}")
        summary_lines.append(f"  - Fire Stations:       {counts.get('fire_station', 0)}")
        summary_lines.append("")
        try:
            render_map(city, records)
        except Exception as e:
            log(f"Could not render map for {city}: {e}")

    report_path = RESULTS_DIR / "emergency_services_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))
    log(f"Saved report to {report_path.relative_to(ROOT)}")

    # Generate side-by-side comparison figure if both city maps exist
    blr_fig = FIGURES_DIR / "bengaluru_emergency_services_dark.png"
    chn_fig = FIGURES_DIR / "chennai_emergency_services_dark.png"
    if blr_fig.exists() and chn_fig.exists():
        try:
            import matplotlib.image as mpimg
            fig, axes = plt.subplots(1, 2, figsize=(20, 10), facecolor=SURFACE)
            img_blr = mpimg.imread(str(blr_fig))
            img_chn = mpimg.imread(str(chn_fig))

            axes[0].imshow(img_blr)
            axes[0].axis("off")
            axes[0].set_title("Bengaluru Emergency Infrastructure", color=INKD, fontsize=14, weight="bold", pad=12)

            axes[1].imshow(img_chn)
            axes[1].axis("off")
            axes[1].set_title("Chennai Emergency Infrastructure", color=INKD, fontsize=14, weight="bold", pad=12)

            plt.tight_layout()
            out = FIGURES_DIR / "compare_emergency_services_dark.png"
            plt.savefig(out, dpi=200, facecolor=SURFACE)
            plt.close()
            log(f"Saved side-by-side comparison: {out.relative_to(ROOT)}")
        except Exception as e:
            log(f"Could not build side-by-side comparison: {e}")


if __name__ == "__main__":
    main()

