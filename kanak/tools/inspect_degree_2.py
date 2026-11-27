# Author: Kanak Varma
"""
Inspect Degree-2 Nodes in Bengaluru and Chennai Street Networks.
Explores why OSMnx simplification retained these nodes and produces detailed
zoomed-in visual plots of the street geometry around them.
"""

import osmnx as ox
import networkx as nx
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

for city, fpath in [("Bengaluru", "data/bengaluru_city_bengaluru_drive.graphml"),
                    ("Chennai", "data/chennai_city_chennai_drive.graphml")]:
    G = ox.load_graphml(fpath)
    S = nx.Graph()
    S.add_nodes_from(G.nodes(data=True))
    S.add_edges_from((u, v) for u, v, _ in G.edges(keys=True) if u != v)
    
    deg2_nodes = [u for u, d in S.degree() if d == 2]
    sc2_nodes = [u for u, d in G.nodes(data=True) if str(d.get("street_count")) == "2"]
    
    print(f"\n{'='*60}")
    print(f"=== {city} Degree-2 Nodes Analysis ===")
    print(f"{'='*60}")
    print(f"Total graph nodes: {G.number_of_nodes():,}")
    print(f"Nodes with undirected graph degree == 2: {len(deg2_nodes)} ({len(deg2_nodes)/G.number_of_nodes():.3%})")
    print(f"Nodes with street_count == 2 attribute: {len(sc2_nodes)}")
    
    reasons = []
    sample_details = []
    
    for u in deg2_nodes:
        ndata = G.nodes[u]
        in_edges = list(G.in_edges(u, data=True))
        out_edges = list(G.out_edges(u, data=True))
        nbrs = list(S.neighbors(u))
        
        # Check reasons:
        # 1. Self loop at node?
        has_self_loop = u in G.neighbors(u)
        
        # 2. Reversal / oneway transition (e.g. 2 in, 0 out or 0 in, 2 out)
        is_source_or_sink = (len(in_edges) == 0 or len(out_edges) == 0)
        
        # 3. Change in highway class across the node
        hw_in = {d.get("highway") for _, _, d in in_edges}
        hw_out = {d.get("highway") for _, _, d in out_edges}
        hw_diff = (hw_in != hw_out)
        
        # 4. Change in oneway attribute
        ow_in = {d.get("oneway") for _, _, d in in_edges}
        ow_out = {d.get("oneway") for _, _, d in out_edges}
        ow_diff = (ow_in != ow_out)
        
        # 5. Node has specialized OSM tag (e.g. traffic signal, barrier, toll, crossing)
        node_hw = ndata.get("highway")
        node_barrier = ndata.get("barrier")
        has_node_tag = bool(node_hw or node_barrier or ndata.get("junction"))
        
        # Classify primary reason
        reason = "Other"
        if is_source_or_sink:
            reason = "Oneway Reversal / Convergence (Cul-de-sac or Direction Change)"
        elif has_node_tag:
            reason = f"Special Node Tag: {node_hw or node_barrier}"
        elif hw_diff:
            reason = "Road Classification Boundary (e.g. primary -> secondary)"
        elif ow_diff:
            reason = "Oneway-to-Twoway Transition"
        else:
            reason = "Complex Geometry / Topology Boundary"
            
        reasons.append(reason)
        
        names = []
        for _, _, d in in_edges + out_edges:
            nm = d.get("name")
            if nm and nm not in names:
                names.append(str(nm))
        st_name = "; ".join(names) if names else "Unnamed street"
        
        sample_details.append({
            "node_id": u,
            "lat": ndata.get("y"),
            "lon": ndata.get("x"),
            "reason": reason,
            "in_deg": len(in_edges),
            "out_deg": len(out_edges),
            "street_count": ndata.get("street_count"),
            "node_highway": node_hw,
            "streets": st_name
        })
        
    df_reasons = pd.Series(reasons).value_counts()
    print("\nBreakdown of Why Degree-2 Nodes Were Retained:")
    for r, cnt in df_reasons.items():
        print(f"  - {r}: {cnt} ({cnt/len(deg2_nodes):.1%})")
        
    print("\nSample Top 5 Degree-2 Nodes:")
    df_sample = pd.DataFrame(sample_details)
    for idx, row in df_sample.head(5).iterrows():
        print(f"  #{idx+1} Node {row['node_id']} at ({row['lat']:.5f}, {row['lon']:.5f}):")
        print(f"     Reason: {row['reason']}")
        print(f"     Directed edges: {row['in_deg']} in, {row['out_deg']} out | street_count tag = {row['street_count']}")
        print(f"     Streets: {row['streets']}")

