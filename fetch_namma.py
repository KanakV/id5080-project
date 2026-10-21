"""
Namma Graph — fetch street networks for Chennai & Bengaluru at three scales.

    pip install osmnx          # needs osmnx >= 2.0 (brings networkx, geopandas, shapely)
    python fetch_namma.py                  # fetch everything (metro scale takes a while)
    python fetch_namma.py bengaluru        # only Bengaluru
    python fetch_namma.py bengaluru nbhd   # only Bengaluru neighbourhoods

Writes <city>_<scale>_<name>_drive.graphml into ./data/ and prints basic stats.
Nodes = intersections / dead ends, edges = street segments with a 'length' attribute (metres).
"""
import os, sys, time
import osmnx as ox

ox.settings.use_cache = True          # re-runs are free
ox.settings.log_console = False
ox.settings.requests_timeout = 600    # metro queries are big

NETWORK_TYPE = "drive"                # drivable public streets; use "all" if you also want lanes/footpaths

# --- the three scales ------------------------------------------------------
# nbhd : Nominatim place polygons (planned grid vs old organic core)
# city : municipal corporation boundary
# metro: a fixed 30 km radius around the city centre -> comparable across the two cities
CITIES = {
    "chennai": {
        "nbhd":  {"tnagar":   "T. Nagar, Chennai, Tamil Nadu, India",
                  "mylapore": "Mylapore, Chennai, Tamil Nadu, India"},
        "city":  {"chennai":  "Chennai, Tamil Nadu, India"},
        "metro": {"chennai30km": (13.0827, 80.2707)},
    },
    "bengaluru": {
        "nbhd":  {"jayanagar": "Jayanagar, Bengaluru, Karnataka, India",
                  "chickpet":  "Chickpet, Bengaluru, Karnataka, India"},
        "city":  {"bengaluru": "Bengaluru, Karnataka, India"},
        "metro": {"bengaluru30km": (12.9716, 77.5946)},
    },
}
METRO_RADIUS_M = 30_000


def fetch(city, scale, name, spec):
    out = f"data/{city}_{scale}_{name}_{NETWORK_TYPE}.graphml"
    if os.path.exists(out):
        print(f"  {out} exists, skipping"); return
    t = time.time()
    if scale == "metro":
        G = ox.graph_from_point(spec, dist=METRO_RADIUS_M, dist_type="bbox",
                                network_type=NETWORK_TYPE, simplify=True)
    else:
        G = ox.graph_from_place(spec, network_type=NETWORK_TYPE, simplify=True)
    ox.save_graphml(G, out)
    s = ox.basic_stats(G)
    print(f"  {out}: n={G.number_of_nodes():,}  m={G.number_of_edges():,}  "
          f"streets={s['street_length_total']/1000:,.0f} km  "
          f"mean segment={s['street_length_avg']:.0f} m  "
          f"circuity={s['circuity_avg']:.3f}  [{time.time()-t:.0f}s]")


if __name__ == "__main__":
    os.makedirs("data", exist_ok=True)
    want_city  = sys.argv[1] if len(sys.argv) > 1 else None
    want_scale = sys.argv[2] if len(sys.argv) > 2 else None
    for city, scales in CITIES.items():
        if want_city and city != want_city: continue
        for scale, items in scales.items():
            if want_scale and scale != want_scale: continue
            print(f"{city} / {scale}")
            for name, spec in items.items():
                try:
                    fetch(city, scale, name, spec)
                except Exception as e:
                    print(f"  {name} FAILED: {type(e).__name__}: {str(e)[:200]}")
    print("done — upload the data/*.graphml files to continue the analysis.")