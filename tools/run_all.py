# Author: Kanak Varma
"""
Plug-and-play pipeline runner: fetch a city's street network (once) and run
every analysis step on it, in dependency order.

    python tools/run_all.py bengaluru                 # whole pipeline
    python tools/run_all.py chennai                   # same thing for Chennai
    python tools/run_all.py chennai --steps core,maps # just those steps
    python tools/run_all.py data/foo.graphml          # any graphml, skip fetching
    python tools/run_all.py --list                    # show the steps

Each step is an ordinary script you can also run by hand (see README.md);
this runner just chains them with the right arguments. A step that fails
stops the run with the failing command printed.
"""
import argparse
import os
import subprocess
import sys
import time

# All analysis scripts live alongside this runner in tools/
TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))

STEPS = [  # (name, script, extra args, what it produces)
    ("core",        "analyze_network.py",     ["{k}"], "betweenness + cyclic structure, 2 maps"),
    ("expander",    "expander_bottleneck.py", [],      "lambda_2 / Cheeger / Ramanujan, Fiedler + min cuts, 3 maps"),
    ("models",      "model_comparison.py",    [],      "ER + WS null models, metrics table + panels"),
    ("depth",       "depth_diameter.py",      [],      "eccentricity, EXACT diameter (slow last step), gamma check"),
    ("percolation", "percolation.py",         [],      "robustness curves (random / degree / betweenness)"),
    ("perc-models", "percolation_models.py",  [],      "percolation vs ER + WS chart (needs models+percolation)"),
    ("degrees",     "degree_distribution.py", [],      "degree histogram chart (needs models)"),
    ("maps",        "render_dark_maps.py",    [],      "dark-theme betweenness / cycles / eccentricity maps"),
    ("edges-map",   "betweenness_edges.py",   [],      "light-theme edge-painted betweenness map"),
    ("assort",      "degree_correlations.py", [],      "assortativity r, e_jk, k_nn, local assortativity maps (needs core)"),
    ("hubs",        "hub_bottlenecks.py",     [],      "hub-to-hub effective resistance, corridor Fiedler cut, min cut"),
    ("edits",       "edit_distance.py",       [],      "edit distance to Barcelona / ER / scale-free / ring-spoke (needs assort)"),
    ("export",      "export_share_data.py",   [],      "plain-CSV hand-off bundle (share_<city>/ + zip)"),
]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("city", nargs="?", help="city key (e.g. bengaluru, chennai) or a .graphml path")
ap.add_argument("--steps", default="all", help="comma-separated subset of steps (default: all)")
ap.add_argument("--k", default="300", help="sampled sources for betweenness (default 300)")
ap.add_argument("--list", action="store_true", help="list steps and exit")
args = ap.parse_args()

if args.list or not args.city:
    print("steps, in run order:")
    for name, script, _, desc in STEPS:
        print(f"  {name:12s} tools/{script:26s} {desc}")
    sys.exit(0)

if args.city.endswith(".graphml"):
    graph = args.city
else:
    graph = f"data/{args.city}_city_{args.city}_drive.graphml"
    if not os.path.exists(graph):
        print(f"== fetch: {graph} missing, downloading via fetch_namma.py (needs internet once)")
        subprocess.run([sys.executable, "-u", os.path.join(TOOLS_DIR, "fetch_namma.py"), args.city, "city"], check=True)
if not os.path.exists(graph):
    sys.exit(f"graph not found: {graph}")

wanted = [s.strip() for s in args.steps.split(",")] if args.steps != "all" else [s[0] for s in STEPS]
unknown = set(wanted) - {s[0] for s in STEPS}
if unknown:
    sys.exit(f"unknown steps: {sorted(unknown)} (use --list)")

for name, script, extra, _ in STEPS:
    if name not in wanted:
        continue
    cmd = [sys.executable, "-u", os.path.join(TOOLS_DIR, script), graph] + [a.format(k=args.k) for a in extra]
    print(f"\n== {name}: {' '.join(cmd)}", flush=True)
    t = time.time()
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(f"step '{name}' failed (exit {r.returncode}) — fix and rerun with --steps {name},...")
    print(f"== {name} done [{time.time()-t:.0f}s]", flush=True)

print("\nall requested steps done — see figures/, results/, and ANALYSIS.md")
