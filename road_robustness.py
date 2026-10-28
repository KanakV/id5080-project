#Author: Omkar
"""
Robustness of Chennai & Bengaluru road networks vs ER and scale-free networks
-----------------------------------------------------------------------------
Two removal strategies (nodes):
  1. Random failure
  2. Targeted attack on betweenness centrality

Metric (same as the "Robustness Lab" slides/artifact):
  P_inf(f) / P_inf(0)  = size of giant component after removing a fraction f
                         of nodes, relative to its initial size
  f_c = first f where the giant component drops below 1% of its initial size
  R   = Schneider robustness measure = mean over f of (giant size / N)
        (area under the curve; larger = more robust)

Comparison networks (same N and same number of links as each city):
  ER           : G(N, M) random graph
  Scale-free   : configuration model with p_k ~ k^-gamma

Install once:   pip install osmnx networkx numpy matplotlib
Run:            python road_robustness.py
"""

import os
import time
import numpy as np
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------
CITIES = {
    "Chennai":   "Chennai, Tamil Nadu, India",
    "Bengaluru": "Bengaluru, Karnataka, India",
}
SEED = 42
GAMMA = 2.5            # degree exponent of the scale-free comparison network
N_RANDOM_RUNS = 10     # random failure is averaged over this many runs
BC_SAMPLES = 500       # betweenness is approximated with k sampled sources
                       # (exact betweenness on 20k+ nodes is very slow)
ADAPTIVE = False       # False: rank nodes once by initial betweenness (static)
                       # True : recompute betweenness after every batch removal
ADAPTIVE_STEP = 0.05   # batch size (fraction of N) when ADAPTIVE = True
GRID = np.linspace(0, 1, 401)   # common f grid for plotting/averaging
OUT_DIR = "robustness_plots"

# colours taken from the Robustness Lab artifact
C_SF = "#D62F8C"       # scale-free (pink)
C_ER = "#2F6BFF"       # ER (blue)
C_CITY = {"Chennai": "#FF7A1A", "Bengaluru": "#0B8F62"}

os.makedirs(OUT_DIR, exist_ok=True)
ox.settings.use_cache = True
ox.settings.log_console = False


# ----------------------------------------------------------------------------
# BUILD NETWORKS
# ----------------------------------------------------------------------------
def load_city(place):
    """Drivable road network -> simple undirected graph, largest component,
    nodes relabelled 0..N-1."""
    print(f"  downloading / loading '{place}' ...")
    Gm = ox.graph_from_place(place, network_type="drive")  # simplified by default
    U = nx.Graph()
    for u, v in Gm.edges():
        if u != v:
            U.add_edge(u, v)
    lcc = max(nx.connected_components(U), key=len)
    U = U.subgraph(lcc).copy()
    return nx.convert_node_labels_to_integers(U)


def make_er(n, m, seed):
    return nx.gnm_random_graph(n, m, seed=seed)


def make_scale_free(n, m, gamma, seed):
    """Configuration model with p_k ~ k^-gamma, mean degree tuned to 2m/n."""
    rng = np.random.default_rng(seed)
    u = rng.random(n)
    target = 2.0 * m / n
    kmax = int(n ** (1.0 / (gamma - 1.0)))      # natural cutoff
    kmin = 1.0
    for _ in range(60):                          # tune kmin to hit mean degree
        k = np.clip(np.rint(kmin * (1 - u) ** (-1.0 / (gamma - 1.0))), 1, kmax)
        kmin *= target / k.mean()
    deg = k.astype(int)
    if deg.sum() % 2:
        deg[np.argmax(deg)] += 1
    G = nx.Graph(nx.configuration_model(deg.tolist(), seed=seed))
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


# ----------------------------------------------------------------------------
# REMOVAL ORDERS
# ----------------------------------------------------------------------------
def betweenness_order(G, k, seed, adaptive=False, step=0.05):
    n = G.number_of_nodes()
    if not adaptive:
        bc = nx.betweenness_centrality(G, k=min(k, n), seed=seed)
        return np.array(sorted(bc, key=bc.get, reverse=True))
    H = G.copy()
    batch = max(1, int(step * n))
    order = []
    while H.number_of_nodes() > 0:
        bc = nx.betweenness_centrality(H, k=min(k, H.number_of_nodes()), seed=seed)
        top = sorted(bc, key=bc.get, reverse=True)[:batch]
        order.extend(top)
        H.remove_nodes_from(top)
    return np.array(order)


# ----------------------------------------------------------------------------
# GIANT COMPONENT vs NODES REMOVED (exact, one pass, reverse union-find)
# ----------------------------------------------------------------------------
def giant_curve(G, order):
    """giant[r] = size of the largest component after removing order[:r]."""
    n = G.number_of_nodes()
    adj = [list(G.neighbors(i)) for i in range(n)]
    parent = np.arange(n)
    size = np.ones(n, dtype=np.int64)
    present = np.zeros(n, dtype=bool)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    giant = np.zeros(n + 1, dtype=np.int64)
    best = 0
    for r in range(n, 0, -1):          # add nodes back, last removed first
        giant[r] = best
        v = order[r - 1]
        present[v] = True
        for w in adj[v]:
            if present[w]:
                a, b = find(v), find(w)
                if a != b:
                    if size[a] < size[b]:
                        a, b = b, a
                    parent[b] = a
                    size[a] += size[b]
        best = max(best, size[find(v)])
    giant[0] = best
    return giant


def curve_on_grid(giant):
    """Return (P_inf(f)/P_inf(0) on GRID, R, f_c)."""
    n = len(giant) - 1
    f = np.arange(n + 1) / n
    rel = giant / giant[0]
    R = giant.mean() / n
    below = np.where(rel < 0.01)[0]
    fc = f[below[0]] if len(below) else 1.0
    return np.interp(GRID, f, rel), R, fc


def random_failure(G, runs, rng):
    n = G.number_of_nodes()
    curves, Rs, fcs = [], [], []
    for _ in range(runs):
        c, R, fc = curve_on_grid(giant_curve(G, rng.permutation(n)))
        curves.append(c); Rs.append(R); fcs.append(fc)
    return np.mean(curves, axis=0), np.mean(Rs), np.mean(fcs)


def targeted_attack(G, seed):
    t = time.time()
    order = betweenness_order(G, BC_SAMPLES, seed, ADAPTIVE, ADAPTIVE_STEP)
    print(f"    betweenness ranking done in {time.time() - t:.0f}s")
    return curve_on_grid(giant_curve(G, order))


# ----------------------------------------------------------------------------
# PLOTTING
# ----------------------------------------------------------------------------
def style(ax, title):
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("fraction of nodes removed  f")
    ax.set_ylabel(r"$P_\infty(f)\,/\,P_\infty(0)$")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_city(city, res):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, mode, name in zip(axes, ("random", "attack"),
                              ("Random failure", "Targeted attack (betweenness)")):
        ax.plot(GRID, res["Scale-free"][mode][0], color=C_SF, lw=2.5, label="Scale-free")
        ax.plot(GRID, res["ER"][mode][0], color=C_ER, lw=2.5, label="Erdős–Rényi")
        ax.plot(GRID, res[city][mode][0], color=C_CITY[city], lw=3.2, label=f"{city} roads")
        style(ax, f"{city}: {name}")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"{city}_robustness.png"), dpi=200)


def plot_cities(all_res):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, mode, name in zip(axes, ("random", "attack"),
                              ("Random failure", "Targeted attack (betweenness)")):
        for city in CITIES:
            ax.plot(GRID, all_res[city][city][mode][0], color=C_CITY[city], lw=3, label=city)
        style(ax, f"Chennai vs Bengaluru: {name}")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "chennai_vs_bengaluru.png"), dpi=200)


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
def main():
    rng = np.random.default_rng(SEED)
    all_res, summary = {}, []

    for city, place in CITIES.items():
        print(f"\n=== {city} ===")
        G = load_city(place)
        n, m = G.number_of_nodes(), G.number_of_edges()
        print(f"  N = {n:,}   M = {m:,}   <k> = {2*m/n:.2f}")

        nets = {
            city: G,
            "ER": make_er(n, m, SEED),
            "Scale-free": make_scale_free(n, m, GAMMA, SEED),
        }
        for name, g in nets.items():
            print(f"  {name:<11} N={g.number_of_nodes():,} M={g.number_of_edges():,} "
                  f"<k>={2*g.number_of_edges()/g.number_of_nodes():.2f}")

        res = {}
        for name, g in nets.items():
            print(f"  -> {name}")
            res[name] = {"random": random_failure(g, N_RANDOM_RUNS, rng),
                         "attack": targeted_attack(g, SEED)}
            for mode in ("random", "attack"):
                _, R, fc = res[name][mode]
                summary.append((city, name, mode, R, fc))
        all_res[city] = res
        plot_city(city, res)

    plot_cities(all_res)

    # ---------------- summary table ----------------
    print("\n" + "=" * 72)
    print(f"{'Network':<22}{'Removal':<10}{'R (area)':>12}{'f_c (1% rule)':>16}")
    print("-" * 72)
    for city, name, mode, R, fc in summary:
        label = f"{city} roads" if name == city else f"{name} ({city} match)"
        print(f"{label:<22}{mode:<10}{R:>12.3f}{fc:>16.3f}")
    print("=" * 72)

    # ---------------- verdict ----------------
    print("\nWhich city is more robust? (higher R and higher f_c = more robust)")
    cities = list(CITIES)
    for mode in ("random", "attack"):
        R = {c: all_res[c][c][mode][1] for c in cities}
        fc = {c: all_res[c][c][mode][2] for c in cities}
        winner = max(R, key=R.get)
        print(f"  {mode:<7}: R  " + "  ".join(f"{c}={R[c]:.3f}" for c in cities) +
              "  |  f_c  " + "  ".join(f"{c}={fc[c]:.3f}" for c in cities) +
              f"   ->  {winner} more robust")

    print(f"\nPlots saved in ./{OUT_DIR}/")
    plt.show()


if __name__ == "__main__":
    main()