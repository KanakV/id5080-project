# Namma Graph: Bengaluru vs Chennai road networks

ID5080 Complex Networks course project. Each member's work is in its own folder, with its own
README that says what it does, how to run it and what it found. Every folder runs on its own:
`cd` into it and follow its README.

| Folder | Member | What it does | Report section |
|---|---|---|---|
| [`omkar/`](omkar/README.md) | Omkar Atkuri | Robustness under random failure and targeted attack; path structure, spatial layout and growth over time | *(Omkar to fill)* |
| [`sathyanand/`](sathyanand/README.md) | S Sathyanand | Betweenness, closeness, PageRank and eigenvector centrality, normalised for comparison across the two cities | Centrality normalisation and comparison |
| [`aakash/`](aakash/README.md) | Aakash Rajan | Degree assortativity (global and local) and hub-to-hub bottlenecks | Degree mixing and hub-to-hub bottlenecks |
| `kanak/` | Kanak Varma | *Pending:* data pipeline, data profiles, emergency-service coverage | Data, network construction, emergency-service coverage |
| `jayvardhan/` | Jayvardhan More | *Pending:* Scotland Yard pursuit–evasion simulation | Scotland Yard simulation |
| `shubhkarman/` | Shubhkarman | *Pending* | |

## Setup

```bash
pip install -r requirements.txt
```

The road graphs are downloaded from OpenStreetMap by each folder's scripts (`data/` and
`*.graphml` are git-ignored because they are large and regenerable). Map data © OpenStreetMap
contributors (ODbL).
