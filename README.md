# Namma Graph: Bengaluru vs Chennai road networks

ID5080 Complex Networks course project. Each member's work is in its own folder, with its own
README that says what it does, how to run it and what it found. Every folder runs on its own:
`cd` into it and follow its README.

| Folder | Member | What it does |
|---|---|---|
| [`kanak/`](kanak/README.md) | Kanak Varma | OSM download pipeline, data cleaning and profiling, degree-2 root causes, arterial backbone, emergency-service (police, fire) network-Voronoi coverage |
| [`omkar/`](omkar/README.md) | Omkar Atkuri | Robustness under random failure and targeted attack; path structure, spatial layout and growth over time |
| [`sathyanand/`](sathyanand/README.md) | S Sathyanand | Betweenness, closeness, PageRank and eigenvector centrality, normalised for comparison across the two cities |
| [`aakash/`](aakash/README.md) | Aakash Rajan | Degree assortativity (global and local) and hub-to-hub bottlenecks |
| [`jayvardhan/`](jayvardhan/README.md) | Jayvardhan More | Scotland Yard pursuit–evasion simulation on the game board, model networks and both cities |
| [`shubhkarman/`](shubhkarman/README.md) | Shubhkarman Singh Sandhu | Modularity communities and detours; thermal-conduction influence of junctions |

## Setup

```bash
pip install -r requirements.txt
```

`jayvardhan/` also needs `flask` (see its `requirements.txt`). The road graphs are downloaded
from OpenStreetMap by `kanak/tools/fetch_namma.py` (`data/` and `*.graphml` are git-ignored
because they are large and regenerable). Map data © OpenStreetMap contributors (ODbL).
