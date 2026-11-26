import pandas as pd
import os

def _standardize_columns(df):
    """Strip whitespace, lowercase headers, and map common centrality aliases."""
    df.columns = df.columns.str.strip().str.lower()
    aliases = {
        'id': 'nodeid',
        'node_id': 'nodeid',
        'betweenness': 'between',
        'betweenness_centrality': 'between',
        'close': 'closeness',
        'closeness_centrality': 'closeness',
        'eigenvector_centrality': 'eigenvector',
        'pagerank_centrality': 'pagerank',
    }
    return df.rename(columns=aliases)


def compare_top10_centralities(csv1_path, csv2_path):
    """
    Displays only the Top 10 Mean for each centrality across two networks,
    followed by the node with the highest betweenness centrality for both networks.
    """
    df1 = _standardize_columns(pd.read_csv(csv1_path))
    df2 = _standardize_columns(pd.read_csv(csv2_path))

    N1, N2 = len(df1), len(df2)

    candidate_metrics = [
        ('between_norm', 'between', 'Betweenness'),
        ('closeness_norm', 'closeness', 'Closeness'),
        ('eigenvector_norm', 'eigenvector', 'Eigenvector (Max-Normalized)'),
        ('pagerank_norm', 'pagerank', 'PageRank (Size-Adjusted)')
    ]

    # Include only metrics present in both CSVs
    metrics = []
    for norm_col, raw_col, name in candidate_metrics:
        if raw_col in df1.columns and raw_col in df2.columns:
            metrics.append((norm_col, raw_col, name))
        else:
            print(f"Warning: Skipping '{raw_col}' (not found in both CSVs)")

    # Apply normalizations
    if 'pagerank' in df1.columns and 'pagerank' in df2.columns:
        df1['pagerank_norm'] = df1['pagerank'] * N1
        df2['pagerank_norm'] = df2['pagerank'] * N2

    if 'eigenvector' in df1.columns and 'eigenvector' in df2.columns:
        df1['eigenvector_norm'] = df1['eigenvector'] / df1['eigenvector'].max()
        df2['eigenvector_norm'] = df2['eigenvector'] / df2['eigenvector'].max()

    if 'between' in df1.columns and 'between' in df2.columns:
        df1['between_norm'] = df1['between']
        df2['between_norm'] = df2['between']

    if 'closeness' in df1.columns and 'closeness' in df2.columns:
        df1['closeness_norm'] = df1['closeness']
        df2['closeness_norm'] = df2['closeness']

    # 1. Top 10 Mean comparison table
    summary_rows = []
    for norm_col, _, name in metrics:
        mean10_1 = df1.nlargest(10, norm_col)[norm_col].mean()
        mean10_2 = df2.nlargest(10, norm_col)[norm_col].mean()

        summary_rows.append({
            'Centrality': name,
            'CSV 1 Top 10 Mean': round(mean10_1, 6),
            'CSV 2 Top 10 Mean': round(mean10_2, 6),
        })

    summary_df = pd.DataFrame(summary_rows)
    print("=== Top 10 Mean Centrality ===")
    print(summary_df.to_string(index=False))

    # 2. Highest betweenness centrality node for both networks
    highest_betweenness = {}
    if 'between_norm' in df1.columns and 'between_norm' in df2.columns:
        top1 = df1.nlargest(1, 'between_norm').iloc[0]
        top2 = df2.nlargest(1, 'between_norm').iloc[0]

        highest_betweenness = {
            'CSV 1': {'nodeid': top1['nodeid'], 'score': round(top1['between_norm'], 6)},
            'CSV 2': {'nodeid': top2['nodeid'], 'score': round(top2['between_norm'], 6)},
        }

        print("\n=== Highest Betweenness Centrality Node ===")
        print(f"CSV 1 ({csv1_path}): NodeID = {top1['nodeid']} | Betweenness = {top1['between_norm']:.6f}")
        print(f"CSV 2 ({csv2_path}): NodeID = {top2['nodeid']} | Betweenness = {top2['between_norm']:.6f}")

    return summary_df, highest_betweenness