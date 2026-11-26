# Urban Road Network Centrality Comparison Tool

A Python pipeline that downloads real-world drivable road networks for any two cities using OpenStreetMap (OSM), computes key graph centrality metrics, compares structural network resilience and flow between the two cities, and identifies the exact real-world street junction with the highest traffic bottleneck potential (betweenness centrality) in each city.

---

## Project Structure

```text
├── main.py                    # Entry point: orchestrates inputs, caching checks, and pipeline execution
├── Graph_downloader.py        # Downloads OSM road graphs and exports node/edge CSVs
├── Centrality_finder.py       # Computes Betweenness, Closeness, Eigenvector, and PageRank centralities
├── Winner_determination.py    # Normalizes and compares Top-10 centrality metrics across two cities
└── Junction_finder.py         # Reverse-geocodes the top OSM Node ID into a physical street location
```

---

## Prerequisites & Installation

Ensure you have Python 3.9+ installed, then install the required dependencies:

```bash
pip install osmnx pandas networkx requests numpy scipy
```

---

## How to Run

Execute the main orchestrator script from your terminal:

```bash
python main.py
```

When prompted, enter two city names (including country or state for accurate geocoding):

```text
Enter the first city (e.g., Chennai, India): Chennai, India
Enter the second city (e.g., Bangalore, India): Bangalore, India
```

### Smart Caching / Skip Logic
`main.py` automatically checks your working directory before running expensive network downloads or graph calculations:
1. **If `{City}_centralities.csv` exists:** Skips both the OSM download and centrality calculations, proceeding straight to city comparison.
2. **If `{City}_nodes.csv` and `{City}_edges.csv` exist:** Skips the OSM download and jumps directly to centrality calculation.
3. **If neither exists:** Downloads the road network from scratch and executes the full pipeline.

---

## Module Breakdown

### 1. `main.py` (Pipeline Orchestrator)
Coordinates the end-to-end workflow across all helper modules:
* Prompts the user for two cities and extracts clean filename prefixes (e.g., `"Chennai, India"` $\rightarrow$ `"Chennai"`).
* Checks for existing CSV files to avoid redundant API calls or recomputations.
* Passes the generated `{City}_centralities.csv` files to `Winner_determination.py`.
* Extracts the winning `nodeid` for each city and calls `Junction_finder.py` to print its physical address and coordinates.

---

### 2. `Graph_downloader.py` (Network Extraction)
* **Function:** `download_city_graph(place_name)`
* Uses `osmnx.graph_from_place(..., network_type="drive")` to fetch the directed drivable road network for the target city.
* Converts the graph into GeoDataFrames and exports two minimal, Gephi-compatible CSV files:
  * **`{City}_nodes.csv`**: Contains `Id` (OSM Node ID), `x` (longitude), `y` (latitude), and `Label`.
  * **`{City}_edges.csv`**: Contains `Source`, `Target`, `Type` (`Directed`), and `Weight` (physical road segment length in meters).

---

### 3. `Centrality_finder.py` (Graph Analytics)
* **Function:** `calculate_centralities(nodes_csv, edges_csv, output_csv, k_pivots=300, seed=42)`
* Reconstructs a directed graph (`networkx.DiGraph`) from the exported CSVs and calculates four centrality metrics, saving the results to **`{City}_centralities.csv`**:
  * **Betweenness Centrality (`between`):** Approximated via Brandes source-sampling using $k = 300$ random pivot nodes and physical road distance (`weight`). Identifies critical bridge intersections that lie on the highest number of shortest paths.
  * **Closeness Centrality (`closeness`):** Approximated via Eppstein–Wang pivot sampling on the reversed graph with Wasserman–Faust normalization for disconnected components. Measures how quickly a node can be reached from all other nodes.
  * **Eigenvector Centrality (`eigenvector`):** Computed using inverted edge weights ($w_{\text{inv}} = \frac{1}{\text{length}}$) so shorter road segments represent stronger connectivity. Includes an automatic NumPy fallback if power iteration does not converge.
  * **PageRank (`pagerank`):** Computed using inverted edge weights ($w_{\text{inv}}$) to measure flow concentration toward major arterial junctions.

> **Performance Tip:** You can adjust `K_PIVOTS` at the top of `Centrality_finder.py`. Increasing it improves approximation accuracy on large cities at the cost of longer runtime.

---

### 4. `Winner_determination.py` (Cross-City Comparison)
* **Function:** `compare_top10_centralities(csv1_path, csv2_path)`
* Standardizes column headers and normalizes metrics so cities with different total node counts ($N_1$ vs. $N_2$) can be fairly compared:
  * **PageRank** is size-adjusted by multiplying raw scores by the network size ($N$).
  * **Eigenvector** is max-normalized ($\frac{x_i}{\max(x)}$).
  * **Betweenness** and **Closeness** are already normalized to $[0, 1]$.
* Prints a comparison table of the **Top 10 Mean** scores for each centrality metric and returns the single highest-betweenness `nodeid` for both cities.

---

### 5. `Junction_finder.py` (Reverse Geocoding)
* **Function:** `get_junction_name(node_id)`
* Queries the OpenStreetMap v0.6 REST API (`https://api.openstreetmap.org/api/0.6/node/{node_id}.json`) to retrieve the exact latitude and longitude of the top bottleneck node.
* Passes those coordinates to the Nominatim Reverse Geocoding API to display the human-readable street/junction name.
* *(Standalone Feature)*: Running `python Junction_finder.py` directly also provides an interactive CLI to look up road names between two adjacent node IDs (`get_road_name(u_id, v_id)`) via the Overpass API.

---

## Generated Output Files

For each analyzed city (e.g., `Chennai`), the pipeline produces:
1. `Chennai_nodes.csv` — Road intersections and coordinates.
2. `Chennai_edges.csv` — Directed road segments and lengths.
3. `Chennai_centralities.csv` — Computed centrality scores for every node in the network.
