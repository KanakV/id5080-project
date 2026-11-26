import os
import osmnx as ox
import pandas as pd


def download_city_graph(place_name, output_dir="data"):
    """Downloads the drive network for a city and exports node and edge CSVs to output_dir."""
    os.makedirs(output_dir, exist_ok=True)

    file_name = place_name.split(',')[0].strip()
    print(f"Searching for {file_name}...")
    G = ox.graph_from_place(place_name, network_type="drive")
    print("Obtained graph successfully")

    # Extract nodes and edges as GeoDataFrames
    nodes_gdf, edges_gdf = ox.graph_to_gdfs(G)

    # Prepare minimal nodes CSV
    nodes_df = nodes_gdf.reset_index()[['osmid', 'x', 'y']].copy()
    nodes_df.rename(columns={'osmid': 'Id'}, inplace=True)
    nodes_df['Label'] = nodes_df['Id'].astype(str)

    # Prepare minimal edges CSV
    edges_df = edges_gdf.reset_index()[['u', 'v', 'length']].copy()
    edges_df.rename(columns={
        'u': 'Source',
        'v': 'Target',
        'length': 'Weight'
    }, inplace=True)
    edges_df['Type'] = 'Directed'
    edges_df = edges_df[['Source', 'Target', 'Type', 'Weight']]

    # Export minimal CSVs to output_dir
    nodes_csv = os.path.join(output_dir, f"{file_name}_nodes.csv")
    edges_csv = os.path.join(output_dir, f"{file_name}_edges.csv")
    nodes_df.to_csv(nodes_csv, index=False)
    edges_df.to_csv(edges_csv, index=False)

    print(f"Export complete: {nodes_csv} and {edges_csv} generated.")
    return file_name, nodes_csv, edges_csv