"""
Hub-to-hub connectivity compared across cities, normalised by ROAD distance.

Hub pairs sit at different distances in each city (Chennai's hubs span ~33 km of road,
Bengaluru's ~30 km), so every measure is plotted against the shortest road distance
between the two hubs, and cities are compared only where their distances overlap.

  A  Resistance per km of road  R_len / d_road
       R_len = effective resistance with each segment's resistance = its length (km).
       A single road of length L gives R_len = L, so the ratio is 1; k equal parallel
       roads give L/k, so 1/k. Lower = more parallel capacity between the hubs.
  B  Min cut between the hubs' 1 km areas (road segments) = number of edge-disjoint
     routes (Menger). Hollow = the choke point sits in one hub's own catchment.
  C  Corridor Fiedler conductance phi (worst balanced cut of the corridor between the
     hubs, relative to its size). Lower = sharper bottleneck.
  D  The three measures at MATCHED road distance (the range all cities share),
     as the city's median divided by the first city's (Bengaluru's) median.
Lines in A-C are medians in road-distance bins.

Road distances: Dijkstra on segment length between the snapped hub nodes, computed
once per city and cached to results/<tag>_hub_road_km.csv (needs data/<tag>.graphml).

Usage:  python tools/hub_compare.py [graphml ...]   (default: Bengaluru, Chennai)
Writes figures/compare_hubs_dark.png, results/compare_hubs_report.txt
"""
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
import matplotlib.patheffects
from matplotlib.collections import LineCollection

from graph_utils import classify_cuts, hub_road_km, load_simple, log, snap_hubs, to_km

sys.stdout.reconfigure(encoding="utf-8")
SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
CITY_STYLE = {"Bengaluru": ("#3987e5", "o"), "Chennai": ("#c98500", "s")}

paths = sys.argv[1:] or ["data/bengaluru_city_bengaluru_drive.graphml", "data/chennai_city_chennai_drive.graphml"]
tags = [os.path.splitext(os.path.basename(p))[0] for p in paths]
cities = [t.split("_")[0].capitalize() for t in tags]
REF = cities[0]   # panel D baseline
D = {}
for p, t in zip(paths, tags):
    f = f"results/{t}_hub_bottlenecks.csv"
    if not os.path.exists(f):
        sys.exit(f"{f} missing — run hub_bottlenecks.py on that graph first")
    with open(f, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    D[t] = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    for k in ("dist_km", "R_unit", "R_len_km", "lam2", "n_lam2", "phi", "min_cut"):
        D[t][k] = D[t][k].astype(float)
    where = classify_cuts(list(zip(D[t]["hub_a"], D[t]["hub_b"], D[t]["min_cut"].astype(int))))
    D[t]["cut_limited_by"] = np.array([where[(a, b)] for a, b in zip(D[t]["hub_a"], D[t]["hub_b"])])

    # shortest road distance between each hub pair, and the base map for the hub-pair panel (both cached)
    rf, mf = f"results/{t}_hub_road_km.csv", f"results/{t}_hub_map.npz"
    if not (os.path.exists(rf) and os.path.exists(mf)):
        log(f"{t}: computing road distances between hubs / map geometry ...")
        S = load_simple(p)
        hubs = snap_hubs(S, t.split("_")[0])
        road = hub_road_km(S, hubs)
        with open(rf, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["hub_a", "hub_b", "road_km"])
            w.writerows([(a, b, f"{v:.4f}") for (a, b), v in road.items()])
        lon = np.array([S.nodes[u]["x"] for u in S.nodes])
        lat = np.array([S.nodes[u]["y"] for u in S.nodes])
        P = dict(zip(S.nodes, to_km(lon, lat, lon.mean(), lat.mean())))
        np.savez_compressed(mf, segs=np.array([[P[u], P[v]] for u, v in S.edges], dtype=np.float32),
                            hub_names=np.array([h for h, _, _ in hubs]), hub_xy=np.array([P[n] for _, n, _ in hubs]))
        del S
    zm = np.load(mf)
    D[t]["map"] = (zm["segs"], dict(zip(zm["hub_names"].tolist(), zm["hub_xy"])))
    with open(rf, newline="", encoding="utf-8") as fh:
        road = {(r["hub_a"], r["hub_b"]): float(r["road_km"]) for r in csv.DictReader(fh)}
    D[t]["road_km"] = np.array([road[(a, b)] for a, b in zip(D[t]["hub_a"], D[t]["hub_b"])])
    D[t]["r_per_km"] = D[t]["R_len_km"] / D[t]["road_km"]
    D[t]["choked"] = np.array([s.startswith("near ") for s in D[t]["cut_limited_by"]])

# matched range: road distances every city has
lo = max(D[t]["road_km"].min() for t in tags)
hi = min(D[t]["road_km"].max() for t in tags)
inr = {t: (D[t]["road_km"] >= lo) & (D[t]["road_km"] <= hi) for t in tags}
METRICS = [("r_per_km", "resistance per road km\n(lower = better)"), ("min_cut", "min cut\n(higher = better)"),
           ("phi", "corridor conductance φ\n(higher = better)")]

# ---------------------------------------------------------------- report
rep = ("================= Hub-to-hub connectivity, normalised by road distance =================\n"
       f"{'':44s}" + "".join(f"{c:>13s}" for c in cities) + "\n")


def row(name, fn, fmt):
    global rep
    rep += f"{name:44s}" + "".join(f"{format(fn(t), fmt):>13s}" for t in tags) + "\n"


row("hub pairs", lambda t: len(D[t]["road_km"]), "d")
row("road distance range (km)", lambda t: f"{D[t]['road_km'].min():.1f}-{D[t]['road_km'].max():.1f}", "s")
row("median road / straight-line distance", lambda t: np.median(D[t]["road_km"] / D[t]["dist_km"]), ".2f")
row("median resistance per road km, all pairs", lambda t: np.median(D[t]["r_per_km"]), ".3f")
row("median min cut, all pairs", lambda t: np.median(D[t]["min_cut"]), ".0f")
row("median φ x1e-3, all pairs", lambda t: 1e3 * np.median(D[t]["phi"]), ".2f")
rep += f"\nMatched road distance {lo:.1f}-{hi:.1f} km (the range all cities share):\n"
row("  pairs in range", lambda t: int(inr[t].sum()), "d")
for k, lab in METRICS:
    sc = 1e3 if k == "phi" else 1
    lab = lab.split("\n")[0]
    row(f"  median {lab}" + (" x1e-3" if k == "phi" else ""), lambda t, k=k, sc=sc: sc * np.median(D[t][k][inr[t]]),
        ".3f" if k == "r_per_km" else ".2f")
if REF in cities:
    tr = tags[cities.index(REF)]
    for k, lab in METRICS:
        row(f"  {lab.split(chr(10))[0]} / {REF}",
            lambda t, k=k: np.median(D[t][k][inr[t]]) / np.median(D[tr][k][inr[tr]]), ".2f")
rep += "\nChoked hubs (same min cut to nearly every other hub -> choke point in their own catchment):\n"
for c, t in zip(cities, tags):
    ch = {s[5:]: int(mc) for s, mc in zip(D[t]["cut_limited_by"], D[t]["min_cut"]) if s.startswith("near ")}
    rep += f"  {c:10s} " + (", ".join(f"{h} (cut {k})" for h, k in sorted(ch.items(), key=lambda x: x[1])) or "none") + "\n"
rep += ("\nResistance per road km = length-weighted effective resistance / shortest road distance:\n"
        "1 = as good as a single road, 1/k = as good as k parallel roads of that length.\n")
print(rep)
with open("results/compare_hubs_report.txt", "w", encoding="utf-8") as f:
    f.write(rep)

# ---------------------------------------------------------------- figure
bins = np.geomspace(min(D[t]["road_km"].min() for t in tags) * 0.99, max(D[t]["road_km"].max() for t in tags) * 1.01, 8)


def style_ax(ax, title, xlabel="shortest road distance between hubs (km)"):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=10)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, color=INKD, fontsize=13)
    ax.set_xlabel(xlabel, color=SUBD, fontsize=10)


def road_axis(ax):
    ax.set_xscale("log")
    ax.set_xticks([2, 3, 5, 10, 20, 30])
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.axvspan(lo, hi, color="#2a2a28", zorder=0)


def scatter_with_medians(ax, key, scale=1.0, hollow_choked=False):
    for c, t in zip(cities, tags):
        col, mk = CITY_STYLE[c]
        x, y = D[t]["road_km"], scale * D[t][key]
        hol = D[t]["choked"] if hollow_choked else np.zeros(len(x), bool)
        ax.scatter(x[~hol], y[~hol], color=col, marker=mk, s=42, edgecolors=SURFACE, lw=0.8, zorder=3, alpha=0.85,
                   label=c)
        ax.scatter(x[hol], y[hol], facecolors="none", edgecolors=col, marker=mk, s=42, lw=1.4, zorder=3)
        idx = np.digitize(x, bins)
        bx, by = [], []
        for b in np.unique(idx):
            m = idx == b
            if m.sum() >= 2:
                bx.append(np.sqrt(bins[b - 1] * bins[b]))
                by.append(np.median(y[m]))
        ax.plot(bx, by, "-", color=col, lw=2.5, zorder=4)


fig = plt.figure(figsize=(20, 19), dpi=160)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(3, 6, height_ratios=[1.5 if len(tags) == 2 else 0.95, 1, 1], hspace=0.3, wspace=0.6)

# ---- maps: hub pairs coloured by resistance per road km, ONE shared scale across cities
HIGH_BAD = mpl.colors.LinearSegmentedColormap.from_list("hb", ["#1e88e5", "#ffb300", "#e53935"])
allr = np.concatenate([D[t]["r_per_km"] for t in tags])
mnorm = mpl.colors.Normalize(np.percentile(allr, 2), np.percentile(allr, 98))
halo = [mpl.patheffects.withStroke(linewidth=3, foreground=SURFACE)]
for i, (c, t) in enumerate(zip(cities, tags)):
    span = 6 // len(tags)
    ax = fig.add_subplot(gs[0, span * i:span * i + span])
    ax.set_facecolor(SURFACE)
    segs, hxy = D[t]["map"]
    ax.add_collection(LineCollection(segs, colors="#2f2f2c", linewidths=0.2, zorder=1))
    for k in np.argsort(D[t]["r_per_km"]):                  # worst pairs drawn on top
        a, b, v = D[t]["hub_a"][k], D[t]["hub_b"][k], D[t]["r_per_km"][k]
        ax.plot(*zip(hxy[a], hxy[b]), color=HIGH_BAD(mnorm(v)), lw=1.2 + 3.5 * mnorm(v), zorder=2,
                solid_capstyle="round")
    for j, (n_, xy) in enumerate(sorted(hxy.items(), key=lambda kv: kv[1][0])):
        ax.scatter(*xy, s=70, color="#ff8a3d", zorder=4, edgecolors=SURFACE, linewidths=1.2)
        ax.annotate(n_, xy, xytext=(6, 6 if j % 2 == 0 else -13), textcoords="offset points", color=INKD,
                    fontsize=9, fontweight="bold", zorder=5, path_effects=halo)   # stagger: neighbours alternate
    ax.set_aspect("equal")
    ax.autoscale()
    ax.margins(0.06)
    ax.tick_params(colors=SUBD, labelsize=8)
    ax.set_xlabel("km", color=SUBD, fontsize=9)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.set_title(f"{c}: hub pairs", color=INKD, fontsize=13)
sm = mpl.cm.ScalarMappable(norm=mnorm, cmap=HIGH_BAD)
cb = fig.colorbar(sm, ax=fig.axes[:len(tags)], orientation="horizontal", shrink=0.45, pad=0.06, aspect=40)
cb.set_label("resistance per road km between the pair (same scale for every city; red / thick = more constricted)",
             color=SUBD)
cb.ax.tick_params(colors=SUBD)

axs = np.array([[fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])],
                [fig.add_subplot(gs[2, 0:3]), fig.add_subplot(gs[2, 3:6])]])
ax = axs[0, 0]
style_ax(ax, "A  Resistance between hubs per km of road  (lower = better)")
scatter_with_medians(ax, "r_per_km")
road_axis(ax)
ax.set_ylabel("R_len / road distance   (1 = a single road, 1/k = k parallel roads)", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=10, loc="upper right")

ax = axs[0, 1]
style_ax(ax, "B  Min cut between hub areas  (higher = better)")
scatter_with_medians(ax, "min_cut", hollow_choked=True)
road_axis(ax)
ax.set_ylabel("road segments that must fail to cut the hubs apart", color=SUBD)
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=SUBD, label="choke point in the corridor"),
                   plt.Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=SUBD,
                              label="choke point at one hub (same cut to every hub)")],
          frameon=False, labelcolor=INKD, fontsize=9, loc="upper right")

ax = axs[1, 0]
style_ax(ax, "C  Corridor bottleneck: Fiedler conductance φ  (higher = better)")
scatter_with_medians(ax, "phi", scale=1e3)
road_axis(ax)
ax.set_yscale("log")
ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%g"))
ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(subs=(2, 3, 5)))
ax.yaxis.set_minor_formatter(plt.FormatStrFormatter("%g"))
ax.set_ylabel("φ × 10³  (log)", color=SUBD)

ax = axs[1, 1]
style_ax(ax, f"D  At matched road distance ({lo:.1f}–{hi:.1f} km): relative to {REF}", xlabel="")
others = [(c, t) for c, t in zip(cities, tags) if c != REF]    # the baseline is 1 by definition
w = 0.8 / len(others)
tr = tags[cities.index(REF)]
for i, (c, t) in enumerate(others):
    col, _ = CITY_STYLE[c]
    rel = [np.median(D[t][k][inr[t]]) / np.median(D[tr][k][inr[tr]]) for k, _ in METRICS]
    vals = [1 / v if k == "r_per_km" else v for v, (k, _) in zip(rel, METRICS)]   # all bars: higher = better
    xs = np.arange(len(METRICS)) + (i - (len(others) - 1) / 2) * w
    ax.bar(xs, vals, width=w * 0.9, color=col, edgecolor=SURFACE, linewidth=2, label=f"{c} (n={int(inr[t].sum())})")
    for x_, v in zip(xs, vals):
        ax.text(x_, v, f"{v:.2f}", ha="center", va="bottom", color=INKD, fontsize=9)
ax.axhline(1, color=SUBD, ls=":", lw=1)
ax.set_xticks(range(len(METRICS)))
ax.set_xticklabels([f"resistance per road km\n({REF} ÷ city)", f"min cut\n(city ÷ {REF})",
                    f"corridor conductance φ\n(city ÷ {REF})"], color=INKD, fontsize=10)
ax.set_ylabel(f"performance relative to {REF}   (1 = equal, higher = better)", color=SUBD)
ax.grid(False, axis="x")
ax.legend(frameon=False, labelcolor=INKD, fontsize=10, loc="upper left")

fig.suptitle("Hub-to-hub connectivity, normalised by road distance — " + " vs ".join(cities)
             + "\n(shaded band = road distances all cities share; lines = medians per distance bin)",
             color=INKD, fontsize=16, y=0.995)
fig.savefig("figures/compare_hubs_dark.png", facecolor=SURFACE, bbox_inches="tight")
print("wrote figures/compare_hubs_dark.png, results/compare_hubs_report.txt")
