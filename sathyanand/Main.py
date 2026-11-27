import os
from Graph_downloader import download_city_graph
from Centrality_finder import calculate_centralities
from Winner_determination import compare_top10_centralities
from Junction_finder import get_junction_name

DATA_DIR = "data"

def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    city1_input = input("Enter the first city (e.g., Chennai, India): ").strip()
    city2_input = input("Enter the second city (e.g., Bengaluru, India): ").strip()

    cities = [city1_input, city2_input]
    centrality_csvs = []

    for place_name in cities:
        city_prefix = place_name.split(',')[0].strip()
        nodes_csv = os.path.join(DATA_DIR, f"{city_prefix}_nodes.csv")
        edges_csv = os.path.join(DATA_DIR, f"{city_prefix}_edges.csv")
        centrality_csv = os.path.join(DATA_DIR, f"{city_prefix}_centralities.csv")

        # If centralities CSV already exists in data/, skip download and calculation
        if os.path.exists(centrality_csv):
            print(f"[{city_prefix}] Found existing '{centrality_csv}'. Skipping download and calculation.")
        else:
            # 1. Check if nodes and edges CSVs already exist in data/ before downloading
            if os.path.exists(nodes_csv) and os.path.exists(edges_csv):
                print(f"[{city_prefix}] Found existing '{nodes_csv}' and '{edges_csv}'. Skipping download.")
            else:
                _, nodes_csv, edges_csv = download_city_graph(place_name, output_dir=DATA_DIR)

            # 2. Compute centralities and export data/{city}_centralities.csv
            calculate_centralities(nodes_csv, edges_csv, centrality_csv)

        centrality_csvs.append((city_prefix, centrality_csv))

    (city1_name, csv1_path), (city2_name, csv2_path) = centrality_csvs

    # 3. Compare the two cities' centralities
    print(f"\nComparing {city1_name} (CSV 1) vs. {city2_name} (CSV 2)...\n")
    _, highest_nodes = compare_top10_centralities(csv1_path, csv2_path)

    # 4. Resolve and print the physical location of the highest-betweenness node for each city
    print("\n=== Highest Betweenness Junction Locations ===")
    for label, city_name in [("CSV 1", city1_name), ("CSV 2", city2_name)]:
        if label in highest_nodes:
            top_node_id = int(highest_nodes[label]["nodeid"])
            top_score = highest_nodes[label]["score"]
            print(f"\n{city_name} — Node ID: {top_node_id} (Betweenness: {top_score:.6f})")
            print(get_junction_name(top_node_id))


if __name__ == "__main__":
    main()
