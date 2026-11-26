import random
import heapq
from collections import deque
import pandas as pd
import networkx as nx
import os

# Configuration for pivot sampling
K_PIVOTS = 300  # Increase for higher precision; decrease for faster runtime
SEED = 42


def pivot_betweenness_centrality(G, k=500, weight='weight', seed=42):
    """Approximate normalized betweenness centrality via Brandes source sampling."""
    rng = random.Random(seed)
    nodes = list(G.nodes())
    n = len(nodes)
    if n <= 2:
        return {v: 0.0 for v in nodes}

    k = min(k, n)
    pivots = rng.sample(nodes, k)
    betweenness = {v: 0.0 for v in nodes}

    for s in pivots:
        stack = []
        pred = {v: [] for v in nodes}
        sigma = {v: 0.0 for v in nodes}
        dist = {v: float('inf') for v in nodes}

        sigma[s] = 1.0
        dist[s] = 0.0

        # 1. Forward pass: shortest-path DAG
        if weight is None:
            queue = deque([s])
            while queue:
                v = queue.popleft()
                stack.append(v)
                for w in G.neighbors(v):
                    if dist[w] == float('inf'):
                        dist[w] = dist[v] + 1
                        queue.append(w)
                    if dist[w] == dist[v] + 1:
                        sigma[w] += sigma[v]
                        pred[w].append(v)
        else:
            pq = [(0.0, s)]
            while pq:
                d, v = heapq.heappop(pq)
                if d > dist[v]:
                    continue
                stack.append(v)
                for w, edge_data in G[v].items():
                    vw_dist = d + edge_data.get(weight, 1.0)
                    if vw_dist < dist[w]:
                        dist[w] = vw_dist
                        heapq.heappush(pq, (vw_dist, w))
                        sigma[w] = 0.0
                        pred[w] = []
                    if vw_dist == dist[w]:
                        sigma[w] += sigma[v]
                        pred[w].append(v)

        # 2. Backward pass: dependency accumulation
        delta = {v: 0.0 for v in nodes}
        while stack:
            w = stack.pop()
            for v in pred[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                betweenness[w] += delta[w]

    # 3. Scale and normalize for directed graph: N / (k * (N - 1) * (N - 2))
    scale = n / (k * (n - 1) * (n - 2))
    for v in betweenness:
        betweenness[v] *= scale

    return betweenness


def pivot_closeness_centrality(G, k=500, weight='weight', seed=42):
    """Approximate normalized closeness centrality via Eppstein-Wang pivot sampling."""
    rng = random.Random(seed)
    nodes = list(G.nodes())
    n = len(nodes)
    if n <= 1:
        return {v: 0.0 for v in nodes}

    k = min(k, n)
    pivots = rng.sample(nodes, k)

    dist_sums = {v: 0.0 for v in nodes}
    reach_counts = {v: 0 for v in nodes}

    # Reverse directed graph so SSSP from pivot p computes inward distances d(v, p)
    search_graph = G.reverse(copy=False) if G.is_directed() else G

    for p in pivots:
        lengths = nx.single_source_dijkstra_path_length(search_graph, p, weight=weight)
        for v, d in lengths.items():
            dist_sums[v] += d
            reach_counts[v] += 1

    approx_closeness = {}
    for v in nodes:
        r_v = reach_counts[v]
        total_pivot_dist = dist_sums[v]

        # Handle isolated or unreachable nodes (Wasserman-Faust normalization)
        est_reachable = (r_v / k) * n
        if r_v == 0 or est_reachable <= 1.0 or total_pivot_dist <= 0.0:
            approx_closeness[v] = 0.0
        else:
            est_total_dist = (n / k) * total_pivot_dist
            approx_closeness[v] = ((est_reachable - 1.0) ** 2) / ((n - 1) * est_total_dist)

    return approx_closeness


def calculate_centralities(nodes_csv, edges_csv, output_csv, k_pivots=K_PIVOTS, seed=SEED):
    """Reads node and edge CSVs, computes four centrality metrics, and writes the results to output_csv."""
    # 1. Read input CSV files
    nodes_df = pd.read_csv(nodes_csv)
    edges_df = pd.read_csv(edges_csv)

    # 2. Initialize directed graph
    G = nx.DiGraph()
    G.add_nodes_from(nodes_df['Id'])

    # Add directed edges with both standard weight (distance) and inverted weight (strength)
    for _, row in edges_df.iterrows():
        w = float(row['Weight'])
        inv_w = 1.0 / w if w != 0 else 0.0
        G.add_edge(row['Source'], row['Target'], weight=w, inv_weight=inv_w)

    # 3. Compute centralities
    betweenness = pivot_betweenness_centrality(G, k=k_pivots, weight='weight', seed=seed)
    closeness = pivot_closeness_centrality(G, k=k_pivots, weight='weight', seed=seed)

    try:
        eigenvector = nx.eigenvector_centrality(G, weight='inv_weight', max_iter=1000)
    except nx.PowerIterationFailedConvergence:
        eigenvector = nx.eigenvector_centrality_numpy(G, weight='inv_weight')

    pagerank = nx.pagerank(G, weight='inv_weight')

    # 4. Assemble output DataFrame using the specified column headers
    output_df = pd.DataFrame({
        'nodeid': list(G.nodes()),
        'between': [betweenness[n] for n in G.nodes()],
        'closeness': [closeness[n] for n in G.nodes()],
        'eigenvector': [eigenvector[n] for n in G.nodes()],
        'pagerank': [pagerank[n] for n in G.nodes()]
    })

    # 5. Save output CSV
    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    output_df.to_csv(output_csv, index=False)
    print(f"Centralities successfully calculated and saved to {output_csv}")
    return output_df