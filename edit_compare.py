"""
Direct comparison of the edit-distance results across cities, from the caches
edit_distance.py writes (results/<tag>_edit_distance.npz + its report) — no graph reload.

  figures/compare_edit_metrics_dark.png  : each robustness / flow metric for every city window
                                           next to Barcelona, ER, scale-free and ring-spoke
  figures/compare_edit_distance_dark.png : edits needed per 1000 edges (lower + upper bound),
                                           and the gap D vs edits for each target, cities overlaid
  results/compare_edit_distance_report.txt

Usage:  python edit_compare.py data/bengaluru_city_bengaluru_drive.graphml data/chennai_city_chennai_drive.graphml
"""
import os
import re
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
SURFACE, INKD, SUBD, GRID, SPINE = "#1a1a19", "#f0efec", "#c3c2b7", "#33332f", "#44443f"
# fixed entity -> colour / marker (same slots as the other comparison figures)
STYLE = {"Bengaluru": ("#3987e5", "o"), "Chennai": ("#c98500", "s"), "barcelona": ("#d55181", "D"),
         "er": ("#d95926", "s"), "sf": ("#9085e9", "v"), "ringspoke": ("#008300", "o")}
LABELS = {"barcelona": "Barcelona", "er": "Erdős–Rényi", "sf": "scale-free", "ringspoke": "ring–spoke"}
MNAME = {"rob_rand": "robustness R, random failure ↑", "rob_attack": "robustness R, degree attack ↑",
         "bridges": "bridge share ↓", "n_lam2": "n·λ₂ spectral connectivity ↑ (log)",
         "mincut": "min cut / √n, flow capacity ↑ (log)", "efficiency": "efficiency d_eucl / d_net ↑",
         "detour": "median detour d_net / d_eucl ↓ (log)", "r": "assortativity r", "k_mean": "mean degree ⟨k⟩"}
LOG = {"n_lam2", "mincut", "detour"}

paths = sys.argv[1:] or ["data/bengaluru_city_bengaluru_drive.graphml", "data/chennai_city_chennai_drive.graphml"]
tags = [os.path.splitext(os.path.basename(p))[0] for p in paths]
cities = [t.split("_")[0].capitalize() for t in tags]
Z, REACHED, LB_RATE = {}, {}, {}
for t in tags:
    f = f"results/{t}_edit_distance.npz"
    if not os.path.exists(f):
        sys.exit(f"{f} missing — run edit_distance.py on that graph first")
    Z[t] = np.load(f, allow_pickle=True)
    # 'reached' and the full-city lower bound per 1000 edges are only in the report tables
    txt = open(f"results/{t}_edit_distance_report.txt", encoding="utf-8").read()
    lb_txt, ub_txt = txt.split("--- LOWER")[1].split("--- UPPER")
    REACHED[t], LB_RATE[t] = {}, {}
    for T, lab in LABELS.items():
        m = re.search(rf"^{re.escape(lab)}\s+(yes|no)\s", ub_txt, re.M)
        if m:
            REACHED[t][T] = m.group(1) == "yes"
        m = re.search(rf"^{re.escape(lab)}\s+[\d,]+\s+[\d,]+\s+([\d.]+)", lb_txt, re.M)
        if m:
            LB_RATE[t][T] = float(m.group(1))
M = [str(x) for x in Z[tags[0]]["metrics"]]
TG = [str(x) for x in Z[tags[0]]["targets"]]


def hist(t, T):
    h = Z[t][f"hist_{T}"]
    return h[:, 0], h[:, 1], h[:, 2], h[:, 3], h[:, 4:]          # edits, adds, removes, D, metrics


def start(t):
    return hist(t, TG[0])[4][0]


def tvec(t, T):
    return Z[t]["target_vectors"][TG.index(T)]


def style_ax(ax, title, grid_axis="y"):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=SUBD, labelsize=9)
    for s_ in ax.spines.values():
        s_.set_color(SPINE)
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, color=INKD, fontsize=11)


# ---------------------------------------------------------------- report
rep = "================= Edit distance — city comparison =================\n"
rep += "Windows: " + "; ".join(f"{c} n={int(Z[t]['n0']):,} m={int(Z[t]['m0']):,} (radius {float(Z[t]['R_win']):.2f} km)"
                               for c, t in zip(cities, tags)) + "\n\n"
cols = cities + [LABELS[T] for T in TG]
rep += f"{'metric':40s}" + "".join(f"{c:>13s}" for c in cols) + "\n"
for i, k in enumerate(M):
    vals = [start(t)[i] for t in tags] + [np.mean([tvec(t, T)[i] for t in tags]) for T in TG]
    rep += f"{MNAME[k]:40s}" + "".join(f"{v:>13.4g}" for v in vals) + "\n"
rep += "(model columns: mean of the per-city builds, each matched to that city's window n and m)\n"

rep += ("\nEdits per 1000 edges to reach each target's robustness & flow\n"
        "UB = greedy constructive edits on the window ((*) = not reached; value = where the move set stopped)\n"
        "LB = exact degree-sequence lower bound on the full city\n")
rep += f"{'target':14s}" + "".join(f"{c + ' UB':>18s}{c + ' LB':>14s}" for c in cities) + "\n"
for T in TG:
    rep += f"{LABELS[T]:14s}"
    for t in tags:
        ex, ad, rm, D, _ = hist(t, T)
        ub = 1000 * ex[-1] / float(Z[t]["m0"])
        mark = "" if REACHED[t].get(T) else " (*)"
        rep += f"{f'{ub:.0f}{mark}':>18s}{LB_RATE[t].get(T, float('nan')):>14.0f}"
    rep += "\n"
rep += "\nAdds / removes at the last checkpoint:\n"
for T in TG:
    rep += f"  {LABELS[T]:14s}" + "   ".join(
        f"{c}: +{int(hist(t, T)[1][-1]):,} / -{int(hist(t, T)[2][-1]):,} (D = {hist(t, T)[3][-1]:.2f})"
        for c, t in zip(cities, tags)) + "\n"
print(rep)
with open("results/compare_edit_distance_report.txt", "w", encoding="utf-8") as f:
    f.write(rep)

# ---------------------------------------------------------------- figure 1: metric vectors
ents = [(c, STYLE[c]) for c in cities] + [(LABELS[T], STYLE[T]) for T in TG]
fig, axs = plt.subplots(3, 3, figsize=(17, 12), dpi=170)
fig.patch.set_facecolor(SURFACE)
for ax, (i, k) in zip(axs.ravel(), enumerate(M)):
    style_ax(ax, MNAME[k])
    for x, t in enumerate(tags):
        col, mk = STYLE[cities[x]]
        ax.bar(x, start(t)[i], width=0.7, color=col, edgecolor=SURFACE, linewidth=2)
    for j, T in enumerate(TG):
        x = len(tags) + j
        v = np.array([tvec(t, T)[i] for t in tags])
        col, mk = STYLE[T]
        ax.bar(x, v.mean(), width=0.7, color=col, edgecolor=SURFACE, linewidth=2)
        if len(v) > 1 and np.ptp(v) > 0:                      # spread across the per-city builds
            ax.errorbar(x, v.mean(), yerr=[[v.mean() - v.min()], [v.max() - v.mean()]], color=SUBD,
                        capsize=3, lw=1)
    if k in LOG:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%g"))
        ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    if k == "r":
        ax.axhline(0, color=SUBD, lw=1)
    ax.set_xticks(range(len(ents)))
    ax.set_xticklabels([e[0] for e in ents], rotation=30, ha="right", color=INKD, fontsize=9)
fig.suptitle("Central city windows vs target networks — robustness and flow metrics  "
             "(↑ higher is better, ↓ lower is better; models matched to each city's n, m — whisker = spread)",
             color=INKD, fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig("figures/compare_edit_metrics_dark.png", facecolor=SURFACE)
plt.close(fig)

# ---------------------------------------------------------------- figure 2: how many edits
fig = plt.figure(figsize=(18, 10), dpi=170)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.38, wspace=0.28)

ax = fig.add_subplot(gs[0, :2])
style_ax(ax, "Edits per 1000 edges: bar = greedy upper bound (hatched = not reached), ◆ = exact degree lower bound")
w = 0.8 / len(tags)
for x, t in enumerate(tags):
    col, mk = STYLE[cities[x]]
    for j, T in enumerate(TG):
        ex = hist(t, T)[0]
        ub = 1000 * ex[-1] / float(Z[t]["m0"])
        xpos = j + (x - (len(tags) - 1) / 2) * w
        ax.bar(xpos, ub, width=w * 0.9, color=col, edgecolor=SURFACE if REACHED[t].get(T) else col,
               linewidth=1.5, hatch=None if REACHED[t].get(T) else "///", alpha=1 if REACHED[t].get(T) else 0.45)
        ax.scatter(xpos, LB_RATE[t].get(T, np.nan), marker="D", s=45, color=INKD, zorder=4, edgecolors=SURFACE)
ax.set_xticks(range(len(TG)))
ax.set_xticklabels([LABELS[T] for T in TG], color=INKD)
ax.set_ylabel("edits per 1000 existing edges", color=SUBD)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=STYLE[c][0], label=c) for c in cities]
          + [plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=SUBD, hatch="///", label="not reached"),
             plt.Line2D([], [], marker="D", ls="", color=INKD, label="exact degree lower bound")],
          frameon=False, labelcolor=INKD, fontsize=9)

ax = fig.add_subplot(gs[0, 2:])
style_ax(ax, "Adds vs removes at the last checkpoint (per 1000 edges)")
for x, t in enumerate(tags):
    col, _ = STYLE[cities[x]]
    for j, T in enumerate(TG):
        _, ad, rm, _, _ = hist(t, T)
        m0 = float(Z[t]["m0"])
        xpos = j + (x - (len(tags) - 1) / 2) * w
        ax.bar(xpos, 1000 * ad[-1] / m0, width=w * 0.9, color=col, edgecolor=SURFACE, linewidth=1.5,
               label=f"{cities[x]} added" if j == 0 else None)
        ax.bar(xpos, -1000 * rm[-1] / m0, width=w * 0.9, color=col, alpha=0.45, edgecolor=SURFACE,
               linewidth=1.5, label=f"{cities[x]} removed" if j == 0 else None)
ax.axhline(0, color=SUBD, lw=1)
ax.set_xticks(range(len(TG)))
ax.set_xticklabels([LABELS[T] for T in TG], color=INKD)
ax.set_ylabel("removed (below 0)  |  added (above 0)", color=SUBD)
ax.legend(frameon=False, labelcolor=INKD, fontsize=9, ncol=2)

for j, T in enumerate(TG):
    ax = fig.add_subplot(gs[1, j])
    style_ax(ax, f"gap D to {LABELS[T]}", grid_axis="both")
    for x, t in enumerate(tags):
        col, mk = STYLE[cities[x]]
        ex, _, _, D, _ = hist(t, T)
        ax.plot(1000 * ex / float(Z[t]["m0"]), D, "-", color=col, marker=mk, ms=6, lw=2, label=cities[x])
    ax.axhline(0.10, color=SUBD, ls=":", lw=1.2)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("edits per 1000 edges", color=SUBD, fontsize=9)
    if j == 0:
        ax.set_ylabel("normalised gap D (1 = start, 0 = target)", color=SUBD)
        ax.legend(frameon=False, labelcolor=INKD, fontsize=9)
fig.suptitle("How far " + " and ".join(cities) + " are from each target network", color=INKD, fontsize=14)
fig.savefig("figures/compare_edit_distance_dark.png", bbox_inches="tight", facecolor=SURFACE)
print("wrote figures/compare_edit_metrics_dark.png, figures/compare_edit_distance_dark.png, "
      "results/compare_edit_distance_report.txt")
