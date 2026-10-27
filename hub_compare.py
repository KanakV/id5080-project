"""
Direct comparison of hub-to-hub connectivity across cities, from the outputs
hub_bottlenecks.py writes (results/<tag>_hub_bottlenecks.csv + report) — no graph reload.

Hubs sit at different distances in each city, so every measure is plotted against
the straight-line hub distance: effective resistance (unit and length-weighted) with
its R = a + b·ln(d) fit (the 2-D grid law), min cut between hub areas (filled = a
bottleneck between the hubs, hollow = capped by a hub's own access roads), corridor
Fiedler conductance and size-normalised corridor n·λ₂, and each hub's access roads.

Usage:  python hub_compare.py data/bengaluru_city_bengaluru_drive.graphml data/chennai_city_chennai_drive.graphml [...]
Writes figures/compare_hubs_dark.png, results/compare_hubs_report.txt
"""
import csv
import os
import re
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from graph_utils import classify_cuts

sys.stdout.reconfigure(encoding="utf-8")
SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
CITY_STYLE = {"Bengaluru": ("#3987e5", "o"), "Chennai": ("#c98500", "s"), "Barcelona": ("#d55181", "D")}

paths = sys.argv[1:] or ["data/bengaluru_city_bengaluru_drive.graphml", "data/chennai_city_chennai_drive.graphml",
                         "data/barcelona_city_barcelona_drive.graphml"]
tags = [os.path.splitext(os.path.basename(p))[0] for p in paths]
cities = [t.split("_")[0].capitalize() for t in tags]
D, ACCESS = {}, {}
for t in tags:
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
    rep_txt = open(f"results/{t}_hub_bottlenecks_report.txt", encoding="utf-8").read()
    m = re.search(r"area \(an upper limit on every min cut involving it\):\n\s+(.+)\n", rep_txt)
    ACCESS[t] = [(nm.strip(), int(k)) for nm, k in re.findall(r"([^,]+?) (\d+)(?:,|$)", m.group(1))] if m else []


def fit(t):
    """R_unit = a + b ln d."""
    b, a = np.polyfit(np.log(D[t]["dist_km"]), D[t]["R_unit"], 1)
    return a, b


# ---------------------------------------------------------------- report
rep = "================= Hub-to-hub connectivity — city comparison =================\n"
rep += f"{'':40s}" + "".join(f"{c:>14s}" for c in cities) + "\n"


def row(name, fn, fmt):
    global rep
    rep += f"{name:40s}" + "".join(f"{format(fn(t), fmt):>14s}" for t in tags) + "\n"


row("hubs", lambda t: len({*D[t]["hub_a"], *D[t]["hub_b"]}), "d")
row("hub pairs", lambda t: len(D[t]["dist_km"]), "d")
row("median hub distance (km)", lambda t: np.median(D[t]["dist_km"]), ".1f")
row("R_unit fit: intercept a", lambda t: fit(t)[0], ".2f")
row("R_unit fit: slope b (per ln km)", lambda t: fit(t)[1], ".2f")
row("R_unit predicted at 5 km", lambda t: fit(t)[0] + fit(t)[1] * np.log(5), ".2f")
row("R_unit predicted at 10 km", lambda t: fit(t)[0] + fit(t)[1] * np.log(10), ".2f")
row("median R_len per km of distance", lambda t: np.median(D[t]["R_len_km"] / D[t]["dist_km"]), ".3f")
row("median min cut (segments)", lambda t: np.median(D[t]["min_cut"]), ".0f")
row("min cut set in one hub's catchment", lambda t: np.mean([s.startswith("near ") for s in D[t]["cut_limited_by"]]), ".0%")
row("min cut in the corridor between hubs", lambda t: np.mean(D[t]["cut_limited_by"] == "between hubs"), ".0%")
row("median corridor φ (x1e-3)", lambda t: 1e3 * np.median(D[t]["phi"]), ".2f")
row("median corridor n·λ₂", lambda t: np.median(D[t]["n_lam2"]), ".2f")
row("corridor cut separates the hubs", lambda t: np.mean(D[t]["splits"] == "True"), ".0%")
row("median access roads per hub (1 km)", lambda t: np.median([k for _, k in ACCESS[t]]) if ACCESS[t] else np.nan, ".0f")
rep += "\nChoked hubs (same min cut to nearly every other hub -> choke point in their own catchment):\n"
for c, t in zip(cities, tags):
    ch = {}
    for s, mc in zip(D[t]["cut_limited_by"], D[t]["min_cut"]):
        if s.startswith("near "):
            ch[s[5:]] = int(mc)
    rep += f"  {c:10s} " + (", ".join(f"{h} (cut {k})" for h, k in sorted(ch.items(), key=lambda x: x[1])) or "none") + "\n"
rep += "\nAccess roads leaving each hub's 1 km area (an upper limit on every min cut involving that hub):\n"
for c, t in zip(cities, tags):
    rep += f"  {c:10s} " + ", ".join(f"{n} {k}" for n, k in ACCESS[t]) + "\n"
rep += ("\nReading: R grows like ln(distance) on a 2-D grid; a larger slope b or intercept a means fewer parallel\n"
        "routes at the same distance. A larger min cut / φ / n·λ₂ means a wider corridor between the hubs.\n")
print(rep)
with open("results/compare_hubs_report.txt", "w", encoding="utf-8") as f:
    f.write(rep)

# ---------------------------------------------------------------- figure
def style_ax(ax, title, xlabel="straight-line distance between hubs (km)"):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=9)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, color=INKD, fontsize=11)
    ax.set_xlabel(xlabel, color=SUBD, fontsize=9)


def plain_log_x(ax):
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 5, 10, 20, 30])
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())


fig, axs = plt.subplots(2, 3, figsize=(19, 11), dpi=170)
fig.patch.set_facecolor(SURFACE)
xs = np.geomspace(min(D[t]["dist_km"].min() for t in tags), max(D[t]["dist_km"].max() for t in tags), 50)

ax = axs[0, 0]
style_ax(ax, "Effective resistance R (unit conductance) — line: a + b·ln d")
for c, t in zip(cities, tags):
    col, mk = CITY_STYLE[c]
    a, b = fit(t)
    ax.scatter(D[t]["dist_km"], D[t]["R_unit"], color=col, marker=mk, s=40, edgecolors=SURFACE, lw=0.8, zorder=3,
               label=f"{c}  (b = {b:.2f})")
    xx = xs[(xs >= D[t]["dist_km"].min()) & (xs <= D[t]["dist_km"].max())]
    ax.plot(xx, a + b * np.log(xx), color=col, lw=2)
plain_log_x(ax)
ax.set_ylabel("R between hubs", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

ax = axs[0, 1]
style_ax(ax, "Length-weighted resistance per km of distance")
for c, t in zip(cities, tags):
    col, mk = CITY_STYLE[c]
    ax.scatter(D[t]["dist_km"], D[t]["R_len_km"] / D[t]["dist_km"], color=col, marker=mk, s=40,
               edgecolors=SURFACE, lw=0.8, label=c)
plain_log_x(ax)
ax.set_ylabel("R_len / d   (lower = more parallel capacity)", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

ax = axs[0, 2]
style_ax(ax, "Min cut between 1 km hub areas (filled = choke point in the corridor,\nhollow = choke point in one hub's own catchment)")
for c, t in zip(cities, tags):
    col, mk = CITY_STYLE[c]
    btw = np.array([not s.startswith("near ") for s in D[t]["cut_limited_by"]])
    ax.scatter(D[t]["dist_km"][btw], D[t]["min_cut"][btw], color=col, marker=mk, s=45, edgecolors=SURFACE, lw=0.8,
               label=c, zorder=3)
    ax.scatter(D[t]["dist_km"][~btw], D[t]["min_cut"][~btw], facecolors="none", edgecolors=col, marker=mk, s=45,
               lw=1.4, zorder=3)
plain_log_x(ax)
ax.set_ylabel("road segments", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

ax = axs[1, 0]
style_ax(ax, "Corridor Fiedler conductance φ (lower = sharper bottleneck)")
for c, t in zip(cities, tags):
    col, mk = CITY_STYLE[c]
    ax.scatter(D[t]["dist_km"], D[t]["phi"], color=col, marker=mk, s=40, edgecolors=SURFACE, lw=0.8, label=c)
plain_log_x(ax)
ax.set_yscale("log")
ax.set_ylabel("φ", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

ax = axs[1, 1]
style_ax(ax, "Corridor n·λ₂ (size-normalised connectivity)")
for c, t in zip(cities, tags):
    col, mk = CITY_STYLE[c]
    ax.scatter(D[t]["dist_km"], D[t]["n_lam2"], color=col, marker=mk, s=40, edgecolors=SURFACE, lw=0.8, label=c)
plain_log_x(ax)
ax.set_ylabel("n·λ₂ of the corridor", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9)

ax = axs[1, 2]
style_ax(ax, "Access: road segments leaving each hub's 1 km area", xlabel="")
ypos, labels = 0, []
for c, t in zip(cities, tags):
    col, _ = CITY_STYLE[c]
    for n_, k_ in sorted(ACCESS[t], key=lambda x: x[1]):
        ax.barh(ypos, k_, color=col, edgecolor=SURFACE, linewidth=1.5, height=0.8)
        ax.text(k_, ypos, f" {k_}", va="center", color=INKD, fontsize=8)
        labels.append(f"{n_} ({c[:3]})")
        ypos += 1
    ypos += 0.6
yt, k = [], 0
for c, t in zip(cities, tags):
    yt += [k + i for i in range(len(ACCESS[t]))]
    k += len(ACCESS[t]) + 0.6
ax.set_yticks(yt)
ax.set_yticklabels(labels, color=INKD, fontsize=8)
ax.grid(True, axis="x", color=GRID, lw=0.6)
ax.grid(False, axis="y")
ax.set_xlabel("road segments", color=SUBD, fontsize=9)

fig.suptitle("Hub-to-hub connectivity — " + " vs ".join(cities), color=INKD, fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig("figures/compare_hubs_dark.png", facecolor=SURFACE)
print("wrote figures/compare_hubs_dark.png, results/compare_hubs_report.txt")
