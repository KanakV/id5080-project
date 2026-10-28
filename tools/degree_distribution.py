# Author: Kanak Varma
"""
Degree distribution of the city vs ER and WS — grouped bars, two panels:
raw counts and normalized P(d). Reads the histograms cached by
model_comparison.py in results/<tag>_model_degrees.npz (render-only).

Usage:  python degree_distribution.py [path/to/graph.graphml]
Writes figures/<tag>_degree_distribution_dark.png
"""
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
tag = os.path.splitext(os.path.basename(GRAPH))[0]
CITY = tag.split("_")[0].capitalize()
os.makedirs("figures", exist_ok=True)

SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
COLORS = {CITY: "#3987e5", "ER": "#d95926", "WS": "#199e70"}
ORDER = [CITY, "ER", "WS"]

z = np.load(f"results/{tag}_model_degrees.npz")
dmax = max(len(z[g]) for g in ORDER) - 1
hists = {g: np.pad(z[g].astype(float), (0, dmax + 1 - len(z[g]))) for g in ORDER}
degrees = np.arange(1, dmax + 1)          # GCCs have no degree-0 nodes

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=200)
fig.patch.set_facecolor(SURFACE)
W = 0.26
for ax, normed in zip(axes, (False, True)):
    ax.set_facecolor(SURFACE)
    for i, g in enumerate(ORDER):
        h = hists[g][1:]
        y = h / h.sum() if normed else h
        ax.bar(degrees + (i - 1) * W, y, width=W, color=COLORS[g], label=g)
    ax.set_xlabel("degree  d", color=SUBD)
    ax.set_ylabel("fraction of nodes  P(d)" if normed else "number of nodes", color=SUBD)
    ax.set_title("Normalized" if normed else "Raw counts", color=INKD, fontsize=12)
    ax.set_xticks(degrees)
    ax.grid(True, axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=SUBD)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
axes[0].legend(frameon=False, labelcolor=INKD, fontsize=10)
axes[0].yaxis.set_major_formatter(lambda v, _: f"{v/1000:,.0f}k" if v >= 1000 else f"{v:,.0f}")
fig.suptitle(f"Degree distribution: {CITY} vs ER and WS   (same n; WS ~52% denser: k=4)",
             color=INKD, fontsize=13)
fig.tight_layout()
fig.savefig(f"figures/{tag}_degree_distribution_dark.png", bbox_inches="tight", facecolor=SURFACE)
print(f"[{time.strftime('%H:%M:%S')}] done: figures/{tag}_degree_distribution_dark.png")
