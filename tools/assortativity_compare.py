"""
Side-by-side degree-correlation comparison of several cities, from the caches
degree_correlations.py writes (results/<tag>_assortativity.npz) — no graph reload.

Panels: global r, degree distribution P(k), k_nn(k) / (<k^2>/<k>) (1 = uncorrelated),
e_jk difference (2nd city minus 1st), radial profile of local assortativity,
mean local assortativity by road category.

Usage:  python assortativity_compare.py data/bengaluru_city_bengaluru_drive.graphml data/chennai_city_chennai_drive.graphml [...]
Writes figures/compare_assortativity_dark.png, results/compare_assortativity_report.txt
"""
import os
import sys

import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np

SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
# fixed entity -> colour/marker (dark categorical slots; never re-ordered by rank)
CITY_STYLE = {"Bengaluru": ("#3987e5", "o"), "Chennai": ("#c98500", "s"), "Barcelona": ("#d55181", "D")}
FALLBACK = [("#199e70", "^"), ("#9085e9", "v")]
KMAX_SHOW = 6

paths = sys.argv[1:] or ["data/bengaluru_city_bengaluru_drive.graphml", "data/chennai_city_chennai_drive.graphml"]
tags = [os.path.splitext(os.path.basename(p))[0] for p in paths]
Z = {}
for t in tags:
    f = f"results/{t}_assortativity.npz"
    if not os.path.exists(f):
        sys.exit(f"{f} missing — run degree_correlations.py on that graph first")
    Z[t] = np.load(f, allow_pickle=True)
city = {t: str(Z[t]["city"]) for t in tags}
style = {}
for i, t in enumerate(tags):
    style[t] = CITY_STYLE.get(city[t], FALLBACK[i % len(FALLBACK)])

# ---------------------------------------------------------------- report
rep = "================= Degree correlations — city comparison =================\n"
rep += f"{'':28s}" + "".join(f"{city[t]:>14s}" for t in tags) + "\n"
rows = [("nodes (giant, simple)", lambda z: f"{int(z['n']):,}"),
        ("edges", lambda z: f"{int(z['m']):,}"),
        ("<k>", lambda z: f"{float(z['k_mean']):.3f}"),
        ("<k^2>/<k>", lambda z: f"{float(z['baseline']):.3f}"),
        ("assortativity r", lambda z: f"{float(z['r']):+.4f}"),
        ("slope k_nn vs k", lambda z: f"{float(z['slope']):+.4f}"),
        ("dead ends (k=1) share", lambda z: f"{z['deg_hist'][1] / z['deg_hist'].sum():.1%}"),
        ("k=3 share", lambda z: f"{z['deg_hist'][3] / z['deg_hist'].sum():.1%}"),
        ("k=4 share", lambda z: f"{z['deg_hist'][4] / z['deg_hist'].sum():.1%}" if len(z['deg_hist']) > 4 else "0"),
        ("boundary clusters (>10 km)", lambda z: f"{len(z['clusters'])}")]
for name, fn in rows:
    rep += f"{name:28s}" + "".join(f"{fn(Z[t]):>14s}" for t in tags) + "\n"

rep += "\nk_nn(k) / (<k^2>/<k>)   (1 = no degree correlation; <1 at high k = disassortative hubs)\n"
for k in range(1, KMAX_SHOW + 1):
    rep += f"  k={k}:" + "".join(
        f"{(Z[t]['knn'][k - 1] / float(Z[t]['baseline'])) if k <= len(Z[t]['knn']) else np.nan:>14.3f}"
        for t in tags) + "\n"

rep += "\nMean local assortativity rho*N by road category (95% bootstrap CI)\n"
cats_all = [str(c) for c in Z[tags[0]]["cat_names"]]
for c in cats_all:
    rep += f"  {c:34s}"
    for t in tags:
        names = [str(x) for x in Z[t]["cat_names"]]
        if c in names:
            n_, mu, lo, hi, _ = Z[t]["cat_stats"][names.index(c)]
            rep += f"  {mu:+.3f} [{lo:+.3f},{hi:+.3f}] n={int(n_):,}" if n_ > 1 else "  (none)"
    rep += "\n"

rep += "\nRadial profile: mean rho*N by distance from CBD (km bin centre)\n"
for t in tags:
    rep += f"  {city[t]:10s} " + "  ".join(f"{c:.0f}km:{v:+.2f}" for c, v, _ in Z[t]["radial"]) + "\n"
print(rep)
with open("results/compare_assortativity_report.txt", "w", encoding="utf-8") as f:
    f.write(rep)

# ---------------------------------------------------------------- figure
def style_ax(ax, title):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=9)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, color=INKD, fontsize=12)


fig = plt.figure(figsize=(19, 11), dpi=180)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.42)
nT = len(tags)
w = 0.8 / nT

# 1: global r
ax = fig.add_subplot(gs[0, 0])
style_ax(ax, "Global assortativity r")
for i, t in enumerate(tags):
    c, _ = style[t]
    rv = float(Z[t]["r"])
    ax.bar(i, rv, width=0.6, color=c, edgecolor=SURFACE, linewidth=2)
    ax.annotate(f"{rv:+.4f}", (i, rv), xytext=(0, 6 if rv >= 0 else -6), textcoords="offset points",
                ha="center", va="bottom" if rv >= 0 else "top", color=INKD, fontsize=10)
rs = [float(Z[t]["r"]) for t in tags]
span = max(rs + [0]) - min(rs + [0])
ax.set_ylim(min(rs + [0]) - 0.25 * span, max(rs + [0]) + 0.25 * span)   # room for the labels
ax.axhline(0, color=SUBD, lw=1)
ax.set_xticks(range(nT))
ax.set_xticklabels([city[t] for t in tags], color=INKD)
ax.set_ylabel("Newman r (eq. 7.82)", color=SUBD)

# 2: degree distribution
ax = fig.add_subplot(gs[0, 1])
style_ax(ax, "Degree distribution P(k)")
ks = np.arange(1, KMAX_SHOW + 1)
for i, t in enumerate(tags):
    c, _ = style[t]
    h = Z[t]["deg_hist"].astype(float)
    p = np.array([h[k] / h.sum() if k < len(h) else 0 for k in ks])
    ax.bar(ks + (i - (nT - 1) / 2) * w, p, width=w * 0.92, color=c, label=city[t], edgecolor=SURFACE, linewidth=1)
ax.set_xticks(ks)
ax.set_xlabel("degree k", color=SUBD)
ax.set_ylabel("share of nodes", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

# 3: normalised k_nn
ax = fig.add_subplot(gs[0, 2])
style_ax(ax, "Neighbour degree  $k_{nn}(k)$ / (<k²>/<k>)")
for t in tags:
    c, mk = style[t]
    z = Z[t]
    ok = z["knn_cnt"] >= 5
    kk = z["ks"][ok]
    v = z["knn"][ok] / float(z["baseline"])
    ax.plot(kk, v, "-", color=c, lw=2, marker=mk, ms=8, label=f"{city[t]}  (r = {float(z['r']):+.3f})")
    ax.annotate(city[t], (kk[-1], v[-1]), xytext=(6, 0), textcoords="offset points", color=INKD,
                fontsize=9, va="center")
ax.axhline(1, color=SUBD, ls="--", lw=1.2, label="no correlation")
ax.set_xlabel("node degree k", color=SUBD)
ax.set_ylabel("relative mean neighbour degree", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9, loc="lower left")

# 4: e_jk difference (2nd - 1st)
ax = fig.add_subplot(gs[1, 0])
ax.set_facecolor(SURFACE)
for s_ in ax.spines.values():
    s_.set_color(SPINE)
K = KMAX_SHOW
e = [np.pad(Z[t]["ejk"], ((0, max(0, K + 1 - Z[t]["ejk"].shape[0])),) * 2)[1:K + 1, 1:K + 1] for t in tags[:2]]
if nT >= 2:
    d = e[1] - e[0]
    lim = np.abs(d).max()
    DIV = mpl.colors.LinearSegmentedColormap.from_list("div", ["#3987e5", "#2b2b29", "#e66767"])
    im = ax.imshow(d, origin="lower", extent=(0.5, K + 0.5, 0.5, K + 0.5), cmap=DIV,
                   norm=mpl.colors.TwoSlopeNorm(0, -lim, lim))
    for j in range(K):
        for k in range(K):
            if abs(d[k, j]) > 0.15 * lim:
                ax.text(j + 1, k + 1, f"{d[k, j] * 100:+.1f}", ha="center", va="center", color=INKD, fontsize=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(f"share of edge ends: {city[tags[1]]} − {city[tags[0]]}", color=SUBD)
    cb.ax.tick_params(colors=SUBD)
ax.set_xticks(range(1, K + 1))
ax.set_yticks(range(1, K + 1))
ax.tick_params(colors=SUBD)
ax.set_xlabel("degree j at one end", color=SUBD)
ax.set_ylabel("degree k at the other end", color=SUBD)
ax.set_title(f"$e_{{jk}}$ difference (labels in %-points)", color=INKD, fontsize=12)

# 5: radial profile
ax = fig.add_subplot(gs[1, 1])
style_ax(ax, "Local assortativity vs distance from CBD")
for t in tags:
    c, mk = style[t]
    rx, ry, re = Z[t]["radial"].T
    ax.errorbar(rx, ry, yerr=re, fmt="-", marker=mk, color=c, ecolor=SUBD, capsize=3, ms=7, lw=2,
                label=city[t])
ax.axhline(0, color=SUBD, lw=1)
ax.set_xlabel("distance from CBD (km)", color=SUBD)
ax.set_ylabel("mean ρ·N", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

# 6: by road category
ax = fig.add_subplot(gs[1, 2])
ax.set_facecolor(SURFACE)
for s_ in ax.spines.values():
    s_.set_color(SPINE)
show = [c for c in cats_all if "Important" not in c]
y0 = np.arange(len(show))[::-1]
for i, t in enumerate(tags):
    col, mk = style[t]
    names = [str(x) for x in Z[t]["cat_names"]]
    for yy, cname in zip(y0, show):
        if cname not in names:
            continue
        n_, mu, lo, hi, _ = Z[t]["cat_stats"][names.index(cname)]
        if n_ < 2:
            continue
        yv = yy + (i - (nT - 1) / 2) * 0.25
        ax.errorbar(mu, yv, xerr=[[mu - lo], [hi - mu]], fmt=mk, color=col, ecolor=SUBD, capsize=3, ms=7,
                    label=city[t] if yy == y0[0] else None)
ax.axvline(0, color=SUBD, lw=1)
ax.set_yticks(y0)
ax.set_yticklabels(show, color=INKD, fontsize=9)
ax.tick_params(colors=SUBD)
ax.grid(True, axis="x", color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.set_xlabel("mean ρ·N (95% CI)", color=SUBD)
ax.set_title("By road category", color=INKD, fontsize=12)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=nT)

fig.suptitle("Degree correlations — " + " vs ".join(city[t] for t in tags), color=INKD, fontsize=15, y=0.98)
fig.savefig("figures/compare_assortativity_dark.png", bbox_inches="tight", facecolor=SURFACE)
print("wrote figures/compare_assortativity_dark.png, results/compare_assortativity_report.txt")
