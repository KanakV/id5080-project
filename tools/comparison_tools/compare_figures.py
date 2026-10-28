# Author: Kanak Varma
"""
Side-by-side city comparison figures.

Places matched figures for Bengaluru and Chennai next to each other, one
output PNG per topic, saved to figures/comparisons/.

Usage:
    python tools/comparison_tools/compare_figures.py [--cities blr chn] [--topics all]

Topics available (must exist for all requested cities):
    betweenness, cycles, degree_correlations, assortativity_by_road_type,
    assortativity_extremes, local_assortativity_map, local_assortativity_cbd,
    hub_bottlenecks, edit_distance, edit_map

Examples:
    python tools/comparison_tools/compare_figures.py
    python tools/comparison_tools/compare_figures.py --topics betweenness cycles
    python tools/comparison_tools/compare_figures.py --cities blr chn bcn
"""

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# ── paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]   # project root
FIGURES_IN  = ROOT / "figures"
FIGURES_OUT = ROOT / "figures" / "comparisons"

# ── city registry ──────────────────────────────────────────────────────────────
CITIES = {
    "blr": {
        "tag":   "bengaluru_city_bengaluru_drive",
        "label": "Bengaluru",
    },
    "chn": {
        "tag":   "chennai_city_chennai_drive",
        "label": "Chennai",
    },
    "bcn": {
        "tag":   "barcelona_city_barcelona_drive",
        "label": "Barcelona",
    },
}

# ── topic registry ─────────────────────────────────────────────────────────────
# Each entry: (suffix_in_filename, output_stem, human_title)
TOPICS = {
    "betweenness":             ("betweenness.png",                    "compare_betweenness",             "Betweenness Centrality"),
    "cycles":                  ("cycles.png",                         "compare_cycles",                  "Cyclic Structure"),
    "cycles_dark":             ("cycles_dark.png",                    "compare_cycles_dark",             "Cyclic Structure (dark)"),
    "degree_correlations":     ("degree_correlations_dark.png",       "compare_degree_correlations",     "Degree Correlations"),
    "assortativity_by_road":   ("assortativity_by_road_type_dark.png","compare_assortativity_by_road",   "Assortativity by Road Type"),
    "assortativity_extremes":  ("assortativity_extremes_dark.png",    "compare_assortativity_extremes",  "Assortativity Extremes"),
    "local_assortativity_map": ("local_assortativity_map_dark.png",   "compare_local_assortativity_map", "Local Assortativity Map"),
    "local_assortativity_cbd": ("local_assortativity_cbd_dark.png",   "compare_local_assortativity_cbd", "Local Assortativity (CBD profile)"),
    "hub_bottlenecks":         ("hub_bottlenecks_dark.png",           "compare_hub_bottlenecks",         "Hub Bottlenecks"),
    "edit_distance":           ("edit_distance_dark.png",             "compare_edit_distance",           "Edit Distance"),
    "edit_map":                ("edit_map_dark.png",                  "compare_edit_map",                "Edit Map"),
}

BG_DARK = "#1a1a19"


def load_image(path: Path):
    """Load a PNG and return as numpy array, or None if missing."""
    if not path.exists():
        return None
    return mpimg.imread(str(path))


def make_comparison(cities: list[dict], topic_key: str, suffix: str, out_stem: str, title: str, out_dir: Path):
    """Render one side-by-side figure and save it."""
    images = []
    labels = []
    for city in cities:
        p = FIGURES_IN / f"{city['tag']}_{suffix}"
        img = load_image(p)
        if img is None:
            print(f"  [skip] {p.name} not found")
            return False
        images.append(img)
        labels.append(city["label"])

    n = len(images)
    # Use the tallest image to set figure height; keep each panel at ~7 in wide
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7),
                             facecolor=BG_DARK,
                             constrained_layout=True,
                             gridspec_kw={"wspace": 0.04})
    if n == 1:
        axes = [axes]

    for ax, img, label in zip(axes, images, labels):
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(label, color="white", fontsize=14, pad=4, fontweight="bold")

    fig.suptitle(title, color="white", fontsize=16, fontweight="bold")

    out_path = out_dir / f"{out_stem}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=BG_DARK, edgecolor="none")
    plt.close(fig)
    print(f"  saved → {out_path.relative_to(ROOT)}")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cities",  nargs="+", default=["blr", "chn"],
                    metavar="CITY",
                    help="city keys to compare (default: blr chn); choices: " + ", ".join(CITIES))
    ap.add_argument("--topics",  nargs="+", default=["all"],
                    metavar="TOPIC",
                    help="topics to render (default: all); choices: " + ", ".join(TOPICS))
    args = ap.parse_args()

    # resolve cities
    unknown_cities = [c for c in args.cities if c not in CITIES]
    if unknown_cities:
        sys.exit(f"Unknown city keys: {unknown_cities}. Choose from: {list(CITIES)}")
    selected_cities = [CITIES[c] for c in args.cities]

    # resolve topics
    if args.topics == ["all"]:
        selected_topics = list(TOPICS.items())
    else:
        unknown_topics = [t for t in args.topics if t not in TOPICS]
        if unknown_topics:
            sys.exit(f"Unknown topics: {unknown_topics}. Choose from: {list(TOPICS)}")
        selected_topics = [(k, TOPICS[k]) for k in args.topics]

    FIGURES_OUT.mkdir(parents=True, exist_ok=True)

    city_str = " vs ".join(c["label"] for c in selected_cities)
    print(f"\nComparing: {city_str}")
    print(f"Output:    figures/comparisons/\n")

    ok = skipped = 0
    for topic_key, (suffix, out_stem, title) in selected_topics:
        print(f"[{topic_key}] {title}")
        if make_comparison(selected_cities, topic_key, suffix, out_stem, title, FIGURES_OUT):
            ok += 1
        else:
            skipped += 1

    print(f"\nDone — {ok} saved, {skipped} skipped (missing source figures).")


if __name__ == "__main__":
    main()

