"""
Coordination friendliness: conflict points per intersection, Bengaluru vs Chennai.
Method and references: CONFLICT_POINTS.md

Each intersection (node with >=3 neighbouring "legs") is a circle with an entry
point and an exit point on every leg. A movement is a chord entry(a) -> exit(b), a != b,
allowed only if the directed graph has both edges (one-ways respected). Conflict points are
  DIVERGE : a stream entering a leg splits towards m exits -> m-1 diverge points
  MERGE   : m streams joining one exit -> m-1 merge points
  CROSS   : their chords intersect inside the circle
Indian traffic is left-hand, so on each leg the exit point comes before the entry
point in clockwise order (mirror-image gives identical counts).
"""
import itertools, math
import numpy as np, networkx as nx, osmnx as ox
import matplotlib.pyplot as plt

CITIES = {"Bengaluru": "data/bengaluru_city_bengaluru_drive.graphml",
          "Chennai": "data/chennai_city_chennai_drive.graphml"}


def conflicts(legs):
    """legs: list (clockwise order) of (has_entry, has_exit). Returns (diverge, merge, cross)."""
    pt = {}                                   # leg -> positions of its exit/entry on the circle
    for i, (has_in, has_out) in enumerate(legs):
        pt[i] = (2 * i, 2 * i + 1)            # exit first, entry second (left-hand traffic)
    mv = [(2 * a + 1, 2 * b) for a, (ai, _) in enumerate(legs) if ai
          for b, (_, bo) in enumerate(legs) if bo and a != b]
    # a stream that splits into m exits has m-1 diverge points; likewise m entries -> m-1 merge points
    entries, exits = {p for p, _ in mv}, {q for _, q in mv}
    div = sum(sum(p == e for p, _ in mv) - 1 for e in entries)
    mer = sum(sum(q == x for _, q in mv) - 1 for x in exits)
    cro = 0
    for (p, q), (r, s) in itertools.combinations(mv, 2):
        if p != r and q != s:
            lo, hi = sorted((p, q))
            cro += (lo < r < hi) != (lo < s < hi)
    return div, mer, cro


# Sanity anchors from the traffic-engineering literature (two-way roads, no U-turns):
assert sum(conflicts([(1, 1)] * 3)) == 9    # T-junction: 3 diverge + 3 merge + 3 cross
assert sum(conflicts([(1, 1)] * 4)) == 32   # 4-leg: 8 diverge + 8 merge + 16 cross
assert conflicts([(1, 1)] * 4) == (8, 8, 16)


def bearing(G, n, v):
    x0, y0, x1, y1 = G.nodes[n]["x"], G.nodes[n]["y"], G.nodes[v]["x"], G.nodes[v]["y"]
    return math.atan2((x1 - x0) * math.cos(math.radians(y0)), y1 - y0)


res = {}
for city, path in CITIES.items():
    G = nx.DiGraph(ox.load_graphml(path))
    rows = []
    for n in G.nodes:
        nb = (set(G.predecessors(n)) | set(G.successors(n))) - {n}
        if len(nb) < 3:
            continue
        nb = sorted(nb, key=lambda v: bearing(G, n, v))
        rows.append(conflicts([(G.has_edge(v, n), G.has_edge(n, v)) for v in nb]))
    a = np.array(rows, dtype=float)
    res[city] = dict(n=len(a), div=a[:, 0].mean(), mer=a[:, 1].mean(), cro=a[:, 2].mean(),
                     total=a.sum(1).mean(), median=np.median(a.sum(1)))
    print(city, {k: round(float(v), 2) for k, v in res[city].items()})

with open("results/conflict_points_report.txt", "w") as f:
    for c, r in res.items():
        f.write(f"{c}: intersections={r['n']:,}  mean conflicts={r['total']:.2f} "
                f"(diverge {r['div']:.2f}, merge {r['mer']:.2f}, cross {r['cro']:.2f})  "
                f"median={r['median']:.0f}\n")

fig, ax = plt.subplots(figsize=(6, 4.5))
cols = {"div": "#7fb3d5", "mer": "#f5b041", "cro": "#c0392b"}
labs = {"div": "Diverging", "mer": "Merging", "cro": "Crossing (incl. right turns)"}
for i, (city, r) in enumerate(res.items()):
    bottom = 0
    for k in ("div", "mer", "cro"):
        ax.bar(i, r[k], bottom=bottom, color=cols[k], width=0.5, label=labs[k] if i == 0 else None)
        ax.text(i, bottom + r[k] / 2, f"{r[k]:.1f}", ha="center", va="center", fontsize=9)
        bottom += r[k]
    ax.text(i, bottom + 0.3, f"{bottom:.1f}\n({r['n']:,} intersections)", ha="center", fontsize=10)
ax.set_xticks(range(len(res))); ax.set_xticklabels(res.keys())
ax.set_ylabel("Mean conflict points per intersection")
ax.set_title("Coordination load at intersections (lower = easier)")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8)
ax.set_ylim(0, max(r["total"] for r in res.values()) * 1.2)
fig.tight_layout()
fig.savefig("figures/conflict_points.png", dpi=200)
