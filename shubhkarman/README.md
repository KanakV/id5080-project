# Complex Networks: Road Network Analysis

This project compares two ways of describing the structure of road networks in Chennai and Bengaluru:

1. **Modularity and detours** groups intersections into communities and measures how much trips are lengthened when a road inside a community is removed.
2. **Thermal conduction** treats intersections as thermal nodes and road segments as conductors, then measures how strongly each node influences the rest of the network.

The precomputed comparison plots and analysis outputs are included. Notebook outputs are preserved. No analysis needs to be rerun just to view the figures.

## Which code should I use?

| Goal | Recommended file | Use it when |
|---|---|---|
| Understand the original modularity experiment | `modularity.py` or `modularity.ipynb` | You want the earlier, compact Leiden and detour prototype. |
| Reproduce the modularity-based detour analysis | `modularity_detour_analysis.py` or `modularity_detour_analysis.ipynb` | You want the complete, documented workflow for Chennai and Bengaluru, including communities, edge detours and comparison maps. |
| Reproduce the thermal network analysis | `road_network_thermal_conduction_lean.py` or `road_network_thermal_conduction_lean.ipynb` | You want the streamlined thermal influence workflow and its city comparison. |

The notebooks are easiest to follow: explanatory text is interleaved with code and saved results. The `.py` files are plain-source exports of the corresponding notebooks; their section comments match the notebook sections. Edit the notebook or its Python export, but keep the two versions in sync if both are maintained.

## How to read the comparison plots

Open `comparison_outputs/` for the side-by-side city maps:

- `bengaluru_chennai_detour_comparison.png` / `.svg` compares detour ratios. Brightness and road thickness indicate the mean increase in shortest-path travel time for trips that used a tested road. Both cities use the same colour scale.
- `bengaluru_chennai_thermal_comparison.png` / `.svg` compares size-scaled thermal influence. The common colour scale allows direct comparison between cities; roads are shown as lines and intersections as coloured points.

The `.png` files are convenient for viewing; the `.svg` files are vector graphics for zooming or publication layouts. The underlying methods differ, so compare the patterns and interpretations rather than treating detour and thermal values as the same measurement.

## Output folders

| Folder | Contents |
|---|---|
| `modularity_detour_output/` | Community maps, per-edge detour CSVs, per-community summaries, and the saved city detour maps. |
| `thermal_road_network_output/` | Thermal figures, node and edge data, configuration files and research summaries for each saved city run. |
| `comparison_outputs/` | Final Chennai–Bengaluru comparison figures for the two policies. |

The detour notebook writes to `modularity_detour_output/`. The thermal notebook writes to `thermal_road_network_output/<city>/`, using the selected `PLACE_QUERY`. The cross-city thermal plot reads the Chennai and Bengaluru run data from those folders.

## Running the analyses

These workflows use Python notebooks and geospatial/network packages, including `osmnx`, `igraph`, `networkx`, `numpy`, `pandas`, `scipy`, `geopandas`, `pyproj`, `shapely`, `matplotlib`, `joblib` and `tqdm`. Use an environment with compatible versions of these packages and Jupyter installed.

The network inputs are GraphML files named for each city (for example, `chennai_roads.graphml` and `bengaluru_roads.graphml`). The detour workflow can also use its saved `*_detour_cache.npz` files. Keep the input files beside the notebook or adjust the paths in the configuration section. The thermal notebook can load a GraphML network or download a place through OSMnx by changing `NETWORK_SOURCE`.

The exact detour computation can take a long time. Review and adjust the configuration cell before running it; the defaults test every eligible edge and origin. The thermal workflow also performs sparse solves over the network. Existing notebook results and output folders are supplied for inspection without rerunning either workflow.

## Suggested reading order

1. Read this file to choose a workflow.
2. Open the relevant notebook and read its introduction and configuration section.
3. Follow the numbered method sections, then inspect the output folder and comparison plot.
4. Use the `.py` export when you prefer browsing the analysis as a single source file.
