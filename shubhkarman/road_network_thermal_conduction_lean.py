"""
Python source export of `road_network_thermal_conduction_lean.ipynb`.
The notebook is the primary, sectioned version of this analysis.
"""

# %% [markdown]
# # Road Network as a Thermal Conduction System#
# Every intersection is a thermal node, every road segment a conductor with relative conductance $g_{ij}=1/L_{ij}$, and every node leaks heat to the environment through $h=\beta\,g_{\mathrm{ref}}$ ($g_{\mathrm{ref}}$ = median edge conductance).# With $A=L_g+H$ and $P=A^{-1}$, holding node $s$ at the source temperature gives the normalised equilibrium response $K_{is}=P_{is}/P_{ss}$.# Repeating this for **every** node as the source gives the thermal influence centrality#
# $$C_{\mathrm{in}}(i)=\frac{1}{N-1}\sum_{s\neq i}K_{is}.$$#
# *Theoretical analogy only: relative conductances from road length, not real traffic or real temperatures.*#
# **Outputs:** thermal influence map, influence vs. heat loss $\beta$ (4 panels), and a Bengaluru vs. Chennai comparison (run the notebook once per city first).
# %% [markdown]
# ## Setup
# %%
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

import geopandas as gpd
import osmnx as ox
import pyproj
from shapely.geometry import LineString

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from tqdm.auto import tqdm

# %% [markdown]
# ## Configuration
# %%
# ── Network ───────────────────────────────────────────────────────────────────────────────
NETWORK_SOURCE = "graphml"       # "graphml": load GRAPHML_PATH
                                 # "place"  : download PLACE_QUERY with OSMnx (needs a geocodable name; cached to GRAPHML_PATH)
c = ["chennai_tambaram", "chennai", "bengaluru"]
PLACE_QUERY = c[2]               # e.g. "Piedmont, California, USA" (small, fast test)
GRAPHML_PATH = f"{PLACE_QUERY}_roads.graphml"
NETWORK_TYPE = "drive"           # OSMnx network_type (only used when downloading)
CUSTOM_FILTER = None
RETAIN_ALL = False
TRUNCATE_BY_EDGE = True
MIN_EDGE_LENGTH_M = 1.0          # floor on segment length (avoids infinite conductance)

# ── Thermal model ─────────────────────────────────────────────────────────────────────────
BETA_MAIN = 1e-2                 # heat-loss parameter for the main influence map (> 0)
BETA_PANELS = [1e-4, 3e-3, 3e-2, 1.0]   # the four beta values shown in the influence-vs-heat-loss figure
MEMORY_BUDGET_MB = 1500          # working memory for one block of sparse solves

# ── Output / style ────────────────────────────────────────────────────────────────────────
OUTPUT_DIR = Path(f"thermal_road_network_output/{PLACE_QUERY}")
THEME = "dark"                   # "dark" | "light"
FIG_DPI = 300
SHOW_FIGS = True

# %% [markdown]
# ### Output folders and plotting helpers
# %%
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR, DATA_DIR = OUTPUT_DIR / "figures", OUTPUT_DIR / "data"
for directory in (FIGURE_DIR, DATA_DIR):
    directory.mkdir(parents=True, exist_ok=True)
CITY = PLACE_QUERY.replace("_", " ").title()


def _trunc(name, lo, hi):
    return LinearSegmentedColormap.from_list(name, mpl.colormaps[name](np.linspace(lo, hi, 256)))

THEMES = {
    "dark":  dict(bg="#0B1020", fg="#E9EDF5", muted="#8E99B3", road="#5F6C8F", road_alpha=0.65, cmap=_trunc("inferno", 0.12, 1.0)),
    "light": dict(bg="#FFFFFF", fg="#1B1F2A", muted="#6B7280", road="#9AA3B5", road_alpha=0.70, cmap=mpl.colormaps["viridis"]),
}
TH = THEMES[THEME]
mpl.rcParams.update({
    "figure.facecolor": TH["bg"], "axes.facecolor": TH["bg"], "savefig.facecolor": TH["bg"],
    "text.color": TH["fg"], "axes.labelcolor": TH["fg"], "axes.titlecolor": TH["fg"],
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42, "figure.dpi": 110,
})


def save_fig(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"{stem}.{ext}", dpi=FIG_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.show() if SHOW_FIGS else plt.close(fig)


def auto_node_size(n):
    return float(np.clip(25.0 * (500.0 / n) ** 0.75, 0.3, 40.0))


def draw_roads(ax, segments, lw=0.4, alpha=None):
    ax.add_collection(LineCollection(segments, colors=TH["road"], linewidths=lw,
                                     alpha=TH["road_alpha"] if alpha is None else alpha, rasterized=True, zorder=1))


def finish_map(ax, xy, pad=0.02):
    dx, dy = np.ptp(xy[:, 0]), np.ptp(xy[:, 1])
    ax.set_xlim(xy[:, 0].min() - pad * dx, xy[:, 0].max() + pad * dx)
    ax.set_ylim(xy[:, 1].min() - pad * dy, xy[:, 1].max() + pad * dy)
    ax.set_aspect("equal")
    ax.set_axis_off()


def style_cbar(cb, label, fontsize=10):
    cb.set_label(label, fontsize=fontsize, color=TH["fg"])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=TH["muted"], labelsize=8, length=2)


def inset_cbar(ax, mappable, label, loc=(0.60, 0.05, 0.36, 0.022), **kw):
    """Slim horizontal colour bar inside the map."""
    cax = ax.inset_axes(list(loc))
    cb = ax.figure.colorbar(mappable, cax=cax, orientation="horizontal", **kw)
    style_cbar(cb, label, fontsize=9)
    cax.xaxis.set_label_position("top")

# %% [markdown]
# ## 1. Road network → undirected thermal network#
# OSMnx returns a directed graph in which a two-way street appears as two reciprocal edges. Each *physical* road is kept once: two directed edges are the same road when they join the same nodes and carry the same OSM way id, length and geometry (read in opposite directions). Genuinely distinct parallel roads are kept, self-loops are dropped, and topology is taken only from OSM nodes, so bridges and tunnels that merely cross other roads do not create intersections. Only the largest connected component is analysed.
# %%
def load_road_network():
    if NETWORK_SOURCE == "graphml":
        if not Path(GRAPHML_PATH).exists():
            raise FileNotFoundError(f"GRAPHML_PATH not found: {GRAPHML_PATH}")
        G = ox.load_graphml(GRAPHML_PATH)
    elif NETWORK_SOURCE == "place":
        G = ox.graph_from_place(PLACE_QUERY, network_type=NETWORK_TYPE, simplify=True, retain_all=RETAIN_ALL,
                                truncate_by_edge=TRUNCATE_BY_EDGE, custom_filter=CUSTOM_FILTER)
    else:
        raise ValueError("NETWORK_SOURCE must be 'graphml' or 'place'.")
    if pyproj.CRS.from_user_input(G.graph["crs"]).is_geographic:
        G = ox.project_graph(G)                                  # metres, needed for lengths and plotting
    if NETWORK_SOURCE == "place" and not Path(GRAPHML_PATH).exists():
        ox.save_graphml(G, GRAPHML_PATH)
    return G


def build_thermal_network(G):
    """Undirected thermal network: returns node table, edge GeoDataFrame (u, v, length, geometry) and counts."""
    ids = sorted(G.nodes)
    idx = {n: i for i, n in enumerate(ids)}
    xs = np.array([G.nodes[n]["x"] for n in ids], float)
    ys = np.array([G.nodes[n]["y"] for n in ids], float)
    seen, recs = set(), []
    n_loops = n_dups = 0
    for u, v, d in G.edges(data=True):
        iu, iv = idx[u], idx[v]
        if iu == iv:                                             # a self-loop carries no heat
            n_loops += 1
            continue
        a, b = min(iu, iv), max(iu, iv)
        geom = d.get("geometry")
        coords = None
        if geom is not None and not geom.is_empty:
            coords = np.asarray(geom.coords)[:, :2]
            coords = coords[::-1] if iu > iv else coords         # canonical orientation a -> b
        length = float(d.get("length") or 0.0)
        if length <= 0:
            length = float(geom.length) if geom is not None else float(np.hypot(xs[a] - xs[b], ys[a] - ys[b]))
        osm = d.get("osmid")
        osm = "|".join(sorted(map(str, osm))) if isinstance(osm, (list, tuple, set)) else str(osm)
        sig = (a, b, osm, round(length, 1), None if coords is None else hash((np.round(coords, 1) + 0.0).tobytes()))
        if sig in seen:                                          # reciprocal copy of an already-kept road
            n_dups += 1
            continue
        seen.add(sig)
        if coords is None:
            coords = np.array([[xs[a], ys[a]], [xs[b], ys[b]]])
        recs.append((a, b, max(length, MIN_EDGE_LENGTH_M), coords))

    u0 = np.array([r[0] for r in recs]); v0 = np.array([r[1] for r in recs])
    n0 = len(ids)
    _, lab = connected_components(sp.coo_matrix((np.ones(len(recs)), (u0, v0)), shape=(n0, n0)), directed=False)
    keep = lab == np.bincount(lab).argmax()                      # largest connected component
    new = np.cumsum(keep) - 1
    recs = [r for r in recs if keep[r[0]]]
    edges = gpd.GeoDataFrame(
        {"u": [new[r[0]] for r in recs], "v": [new[r[1]] for r in recs], "length": [r[2] for r in recs],
         "geometry": [LineString(r[3]) for r in recs]}, geometry="geometry", crs=G.graph["crs"])
    nodes = pd.DataFrame({"osmid": np.array(ids, dtype=object)[keep], "x": xs[keep], "y": ys[keep]})
    return nodes, edges, dict(directed_edges=G.number_of_edges(), reciprocal_duplicates=n_dups, self_loops=n_loops)


nodes, edges, counts = build_thermal_network(load_road_network())
N, E = len(nodes), len(edges)
XY = nodes[["x", "y"]].to_numpy()
SEGMENTS = [np.asarray(g.coords)[:, :2] for g in edges.geometry]
print(f"{CITY}: {N:,} intersections, {E:,} undirected road segments "
      f"({counts['directed_edges']:,} directed edges; {counts['reciprocal_duplicates']:,} reciprocal duplicates merged, "
      f"{counts['self_loops']:,} self-loops dropped)")

# %% [markdown]
# ## 2. Thermal influence centrality#
# $A=L_g+H$ is factorised once per $\beta$ with a sparse LU; columns of $P=A^{-1}$ are then obtained in blocks (no dense inverse is ever formed). For each block, $K_{is}=P_{is}/P_{ss}$ is accumulated into the row sums that give $C_{\mathrm{in}}$.
# %%
def conductance_laplacian(nodes, edges):
    n = len(nodes)
    u, v, g = edges["u"].to_numpy(), edges["v"].to_numpy(), 1.0 / edges["length"].to_numpy()
    W = sp.coo_matrix((np.r_[g, g], (np.r_[u, v], np.r_[v, u])), shape=(n, n)).tocsr()   # parallel roads add up
    return sp.diags(np.asarray(W.sum(axis=1)).ravel()) - W, float(np.median(g))


def influence_centrality(Lg, g_ref, beta):
    """C_in(i) = mean over sources s != i of K_is = P_is / P_ss, with P = (Lg + beta*g_ref*I)^-1."""
    n = Lg.shape[0]
    lu = splu((Lg + sp.identity(n) * (beta * g_ref)).tocsc(), permc_spec="MMD_AT_PLUS_A",
              diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
    block = int(np.clip(MEMORY_BUDGET_MB * 1e6 / (8 * n * 4), 4, 1024))
    row_sum = np.zeros(n)
    for start in tqdm(range(0, n, block), desc=f"beta = {beta:g}", leave=False):
        src = np.arange(start, min(start + block, n))
        rhs = np.zeros((n, len(src)))
        rhs[src, np.arange(len(src))] = 1.0
        K = lu.solve(rhs)                                        # columns of P
        K /= K[src, np.arange(len(src))]                         # divide column s by P_ss
        row_sum += K.sum(axis=1)
    return (row_sum - 1.0) / (n - 1)                             # remove the unit diagonal


L_G, G_REF = conductance_laplacian(nodes, edges)
BETAS = sorted(set(BETA_PANELS) | {BETA_MAIN})
C_IN = {b: influence_centrality(L_G, G_REF, b) for b in tqdm(BETAS, desc="heat-loss values")}
print(f"g_ref = {G_REF:.3e} 1/m;  mean C_in at beta = {BETA_MAIN:g}: {C_IN[BETA_MAIN].mean():.4f}")

# %% [markdown]
# ## Figure 1 · Thermal influence map# Nodes coloured by $C_{\mathrm{in}}$, the mean equilibrium response of each node to all other possible sources.
# %%
def plot_influence_map(c_in, beta):
    fig, ax = plt.subplots(figsize=(10.5, 9.3))
    draw_roads(ax, SEGMENTS)
    o = np.argsort(c_in)
    sc = ax.scatter(XY[o, 0], XY[o, 1], c=c_in[o], s=auto_node_size(N), cmap=TH["cmap"],
                    norm=Normalize(*np.percentile(c_in, [1, 99.5])), linewidths=0, zorder=3, rasterized=True)
    finish_map(ax, XY)
    ax.set_title("Thermal influence map", loc="left", fontsize=14, fontweight="semibold", pad=22)
    ax.text(0, 1.012, f"{CITY}  ·  mean equilibrium response of each node to all other possible sources  ·  β = {beta:g}",
            transform=ax.transAxes, fontsize=8.5, color=TH["muted"], va="bottom")
    inset_cbar(ax, sc, r"Thermal influence centrality  $C_{\mathrm{in}}$", extend="both")
    save_fig(fig, "fig01_thermal_influence_map")

plot_influence_map(C_IN[BETA_MAIN], BETA_MAIN)

# %% [markdown]
# ## Figure 2 · Thermal influence centrality vs. heat loss $\beta$# Each panel has its own colour range so the *pattern* is comparable. Weak loss: heat spreads city-wide and centrality varies smoothly with position. Strong loss: influence stays local.
# %%
def plot_influence_vs_beta(betas):
    fig, axes = plt.subplots(1, len(betas), figsize=(4.9 * len(betas), 5.6))
    for ax, b in zip(np.atleast_1d(axes), betas):
        c = C_IN[b]
        o = np.argsort(c)
        draw_roads(ax, SEGMENTS, lw=0.25, alpha=0.4)
        sc = ax.scatter(XY[o, 0], XY[o, 1], c=c[o], s=auto_node_size(N), cmap=TH["cmap"],
                        norm=Normalize(*np.percentile(c, [1, 99.5])), linewidths=0, zorder=3, rasterized=True)
        finish_map(ax, XY)
        ax.set_title(rf"$\beta$ = {b:g}", loc="left", fontsize=13, pad=18)
        ax.text(0, 1.012, f"mean $C_{{in}}$ = {c.mean():.3g}  ·  ≈{1 / np.sqrt(b):.0f} segments", transform=ax.transAxes,
                fontsize=8.5, color=TH["muted"])
        style_cbar(fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.02, shrink=0.8), r"$C_{\mathrm{in}}$", 9)
    fig.suptitle(f"{CITY}: thermal influence centrality vs. heat loss", x=0.01, y=1.04, ha="left", fontsize=15, fontweight="semibold")
    save_fig(fig, "fig06b_beta_spatial_maps")

plot_influence_vs_beta(BETA_PANELS)

# %% [markdown]
# ## Export# The three files below feed the city comparison at the end.
# %%
pd.DataFrame({"node": np.arange(N), "osmid": nodes["osmid"], "x": XY[:, 0], "y": XY[:, 1],
              "C_in": C_IN[BETA_MAIN]}).to_csv(DATA_DIR / "node_thermal_metrics.csv", index=False)
edges.to_file(DATA_DIR / "thermal_network.gpkg", layer="edges", driver="GPKG")
(OUTPUT_DIR / "run_config.json").write_text(json.dumps(dict(place=PLACE_QUERY, beta_main=BETA_MAIN, beta_panels=BETA_PANELS,
                                                     N=N, E=E), indent=2))
print(f"Saved to {OUTPUT_DIR.resolve()}")

# %% [markdown]
# ## Bengaluru vs. Chennai# Both cities side by side on one colour scale (run the notebook for each city first so their output folders exist).
# %%
ROOT = Path("thermal_road_network_output")
CITIES = {"Bengaluru": "bengaluru", "Chennai": "chennai"}      # label -> output sub-folder
METRIC = "C_in"
SHARED_SCALE = True                                            # True: one colour scale for both cities
OUT_STEM = ROOT / "bengaluru_vs_chennai_thermal_influence"


def plot_city_comparison():
    missing = [k for k, f in CITIES.items() if not (ROOT / f / "data" / "thermal_network.gpkg").exists()]
    if missing:
        print(f"Run the notebook first for: {missing}")
        return
    data = {}
    for label, folder in CITIES.items():
        d = ROOT / folder
        gdf = gpd.read_file(d / "data" / "thermal_network.gpkg", layer="edges")
        data[label] = dict(nodes=pd.read_csv(d / "data" / "node_thermal_metrics.csv"),
                           segs=[np.asarray(g.coords)[:, :2] for g in gdf.geometry],
                           beta=json.loads((d / "run_config.json").read_text())["beta_main"])
    betas = {k: v["beta"] for k, v in data.items()}
    if len(set(betas.values())) > 1:
        print(f"Warning: different beta values {betas}; the cities are not strictly comparable.")

    # Scale each city's centrality by its graph size, then use one pooled
    # maximum so both cities share a [0, 1] colour scale.
    scaled = {k: len(v["nodes"]) * v["nodes"][METRIC].to_numpy() for k, v in data.items()}
    common_max = max(float(np.nanmax(x)) for x in scaled.values())
    if not np.isfinite(common_max) or common_max <= 0:
        raise ValueError("Cannot normalize thermal centrality: pooled maximum must be positive.")
    normalized = {k: x / common_max for k, x in scaled.items()}
    shared = Normalize(0.0, 1.0)
    aspects = [np.ptp(v["nodes"]["x"]) / np.ptp(v["nodes"]["y"]) for v in data.values()]
    H = 8.5
    fig = plt.figure(figsize=(H * sum(aspects) + 1.0, H + 1.6))
    gs = GridSpec(1, len(data), width_ratios=aspects, wspace=0.03, left=0.02, right=0.98, top=0.88, bottom=0.14)
    for i, (label, v) in enumerate(data.items()):
        ax = fig.add_subplot(gs[i])
        n = v["nodes"]
        xy, c = n[["x", "y"]].to_numpy(), normalized[label]
        draw_roads(ax, v["segs"])
        o = np.argsort(c)
        sc = ax.scatter(xy[o, 0], xy[o, 1], c=c[o], s=auto_node_size(len(n)), cmap=TH["cmap"],
                        norm=shared, linewidths=0, zorder=3, rasterized=True)
        finish_map(ax, xy)
        ax.set_title(label, loc="left", fontsize=18, fontweight="semibold", pad=26)
        ax.text(0, 1.012, f"N = {len(n):,} intersections  ·  β = {v['beta']:g}", transform=ax.transAxes,
                fontsize=10, color=TH["muted"], va="bottom")

    cb = fig.colorbar(sc, cax=fig.add_axes([0.30, 0.075, 0.40, 0.022]), orientation="horizontal")
    style_cbar(cb, r"Normalized size-scaled thermal centrality  $N C_{\mathrm{in}} / \max(N C_{\mathrm{in}})$", 11)
    fig.legend(handles=[Line2D([0], [0], color=TH["road"], lw=1.6, label="Road segment"),
                        Line2D([0], [0], marker="o", ls="", ms=6, mfc=TH["cmap"](0.75), mec="none", label="Intersection (coloured by C_in)")],
               loc="lower center", ncol=2, frameon=False, labelcolor=TH["fg"], fontsize=10, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Road networks as thermal conduction systems: Bengaluru vs. Chennai", x=0.02, y=0.975, ha="left",
                 fontsize=20, fontweight="semibold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{OUT_STEM}.{ext}", dpi=FIG_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.show()
    print(f"Saved: {OUT_STEM}.png / .pdf / .svg")

plot_city_comparison()

