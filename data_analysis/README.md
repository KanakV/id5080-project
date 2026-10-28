# OpenStreetMap GraphML Data Profile & Completeness Analysis

This directory contains the data quality, completeness, and distribution analysis of the municipal road network GraphML files (e.g. `bengaluru_city_bengaluru_drive.graphml`).

---

## Directory Organization

```
data_analysis/
├── analyze_data.py          # Data profiling pipeline script
├── figures/                 # High-resolution dark-theme charts
│   ├── bengaluru_city_bengaluru_drive_attribute_completeness_dark.png
│   ├── bengaluru_city_bengaluru_drive_road_hierarchy_distribution_dark.png
│   ├── bengaluru_city_bengaluru_drive_node_features_and_junctions_dark.png
│   └── bengaluru_city_bengaluru_drive_road_attributes_profile_dark.png
└── results/                 # Detailed text reports & raw tabular CSVs
    ├── bengaluru_city_bengaluru_drive_data_profile_report.txt
    ├── bengaluru_city_bengaluru_drive_attribute_completeness.csv
    ├── bengaluru_city_bengaluru_drive_road_hierarchy.csv
    ├── bengaluru_city_bengaluru_drive_node_features.csv
    └── bengaluru_city_bengaluru_drive_road_attributes_summary.csv
```

---

## 1. Known vs. Unknown Data (Completeness)

| Category | Attribute | Known Count | Unknown / Missing | % Complete | Notes |
|---|---|---|---|---|---|
| **Node** | `x`, `y` | 155,779 | 0 | **100.0%** | Exact WGS84 geographic coordinates |
| **Node** | `street_count` | 155,779 | 0 | **100.0%** | Physical incident street count before simplification |
| **Node** | `highway` | 328 | 155,451 | **0.21%** | Traffic signals, motorway junctions, crossings |
| **Node** | `junction` | 30 | 155,749 | **0.02%** | Explicit roundabout/circular tags |
| **Edge** | `osmid`, `highway`, `oneway`, `reversed`, `length` | 394,401 | 0 | **100.0%** | Complete topological backbone & street length |
| **Edge** | `name` | 117,983 | 276,418 | **29.9%** | **70.1% of street segments are unnamed** in OSM |
| **Edge** | `geometry` | 113,156 | 281,245 | **28.7%** | Curved road centerlines; 71.3% are straight chords |
| **Edge** | `lanes` | 24,028 | 370,373 | **6.1%** | 93.9% of streets lack lane-count data |
| **Edge** | `maxspeed` | 8,825 | 385,576 | **2.2%** | 97.8% of streets lack posted speed limits |
| **Edge** | `access` | 5,049 | 389,352 | **1.3%** | Private/permissive/gated access tags |
| **Edge** | `bridge` | 2,248 | 392,153 | **0.57%** | Elevated flyovers, viaducts, railway overpasses |
| **Edge** | `ref` | 1,694 | 392,707 | **0.43%** | National/State Highway codes (e.g., NH 44) |
| **Edge** | `width` / `est_width` | 1,545 | 392,856 | **0.39%** | Physical roadway width measurements |
| **Edge** | `tunnel` | 188 | 394,213 | **0.05%** | Underpasses, subterranean vehicular tunnels |

---

## 2. Road Network Hierarchy

| OSM Classification | Role | Segments | % of Edges | Total Length | % of Length |
|---|---|---|---|---|---|
| **residential** | Local neighborhood streets | 327,564 | **83.1%** | 19,798.9 km | **81.6%** |
| **tertiary** | Neighborhood connectors | 39,711 | **10.1%** | 2,322.7 km | **9.6%** |
| **secondary** | Intermediate arterials | 13,509 | **3.4%** | 875.8 km | **3.6%** |
| **living_street** | Pedestrian-shared residential | 5,439 | **1.4%** | 276.4 km | **1.1%** |
| **primary** | Main city trunk arterials | 4,215 | **1.1%** | 464.8 km | **1.9%** |
| **unclassified** | Minor interconnecting roads | 1,348 | **0.3%** | 129.9 km | **0.5%** |
| **trunk** | High-capacity express corridors | 1,070 | **0.3%** | 253.0 km | **1.0%** |
| **links / ramps** | Grade-separated ramps (all classes) | 1,458 | **0.4%** | 82.8 km | **0.3%** |
| **motorway** | Controlled-access expressways | 61 | **0.02%** | 56.0 km | **0.2%** |

> **Key takeaway**: Bengaluru's street network is over **84% residential / living streets**. High-capacity corridors (`primary`, `trunk`, `motorway`) make up less than **3% of the physical length**, explaining why congestion and betweenness loads concentrate so severely on a few links.

---

## 3. Intersection Topologies & Traffic Signals

### Intersection Types (by Incident Street Count)
- **T-junctions (`street_count = 3`)**: **96,348 (61.8%)** — the dominant Indian urban junction type.
- **Dead ends (`street_count = 1`)**: **38,262 (24.6%)** — 1 in 4 nodes is a cul-de-sac.
- **4-way Crossroads (`street_count = 4`)**: **20,894 (13.4%)** — only ~13% form regular grid intersections.
- **Complex Multi-way (`street_count = 5 or 6`)**: **126 (0.08%)** — max physical junction degree is 6.

### Specialized Infrastructure Tags (Nodes)
- **`traffic_signals`**: **232 nodes** explicitly mapped with traffic lights.
- **`motorway_junction`**: **19 nodes** tagged as highway interchanges.
- **`crossing` / `bus_stop`**: **60 nodes** mapped at intersections.
- **`mini_roundabout` / `turning_circle`**: **13 nodes**.

---

## 4. How to Run

```bash
# Analyze Bengaluru (default)
.venv/bin/python data_analysis/analyze_data.py

# Analyze any other city graphml
.venv/bin/python data_analysis/analyze_data.py data/chennai_city_chennai_drive.graphml
```
