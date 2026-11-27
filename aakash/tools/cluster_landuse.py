"""
What is on the ground at the high-assortativity boundary clusters? Queries OSM
(Overpass) for land use, building types, offices and named industrial/commercial
features within RADIUS_M of each cluster centroid found by degree_correlations.py.

Usage:  python cluster_landuse.py [path/to/graph.graphml]
Reads results/<tag>_assortativity.npz, writes results/<tag>_cluster_landuse.txt
(needs internet; OSM land use is patchy at city edges — treat labels as indicative).
"""
import collections
import os
import sys
import time

import numpy as np
import requests

GRAPH = sys.argv[1] if len(sys.argv) > 1 else "data/bengaluru_city_bengaluru_drive.graphml"
tag = os.path.splitext(os.path.basename(GRAPH))[0]
RADIUS_M = 700
URLS = ["https://overpass-api.de/api/interpreter",          # main server rate-limits bursts,
        "https://overpass.kumi.systems/api/interpreter"]    # so alternate with a public mirror
HEADERS = {"User-Agent": "id5080-namma-graph/1.0 (student network analysis)"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def overpass(la, lo):
    a = f"(around:{RADIUS_M},{la},{lo})"
    q = (f'[out:json][timeout:60];(way{a}[landuse];relation{a}[landuse];'
         f'way{a}[building~"industrial|commercial|office|residential|apartments|house|detached|retail"];'
         f'nwr{a}[office];nwr{a}[industrial];nwr{a}[shop];);out tags center;')
    for attempt in range(6):
        url = URLS[attempt % len(URLS)]
        try:
            r = requests.post(url, data={"data": q}, headers=HEADERS, timeout=120)
            if r.status_code == 200:
                return r.json()["elements"]
            log(f"  {url.split('/')[2]} HTTP {r.status_code}, retry {attempt + 1}")
        except Exception as e:
            log(f"  {url.split('/')[2]} {type(e).__name__}, retry {attempt + 1}")
        time.sleep(15 * (attempt // len(URLS) + 1))
    return None


def verdict(lu, bl, ind_named):
    """Coarse label from tag counts; 'unclear' when OSM has too little."""
    res = lu.get("residential", 0) + bl.get("apartments", 0) / 5 + bl.get("residential", 0) / 10
    ind = lu.get("industrial", 0) * 2 + bl.get("industrial", 0) / 2 + ind_named
    if res + ind < 3:
        return "unclear (sparse OSM)"
    if ind >= 2 * res:
        return "industrial"
    if res >= 2 * ind:
        return "residential layout"
    return "mixed industrial / residential"


z = np.load(f"results/{tag}_assortativity.npz", allow_pickle=True)
clusters = z["clusters"]            # id, lat, lon, n, km, <k>, k4, allj
out = (f"================ Land use at high-assortativity boundary clusters: {tag} ================\n"
       f"(OSM features within {RADIUS_M} m of each cluster centroid; counts are tagged features, not areas)\n")
for cid, la, lo, n_, km, kmean, k4, allj in clusters:
    log(f"cluster #{int(cid)} ({la:.4f}, {lo:.4f})")
    els = overpass(la, lo)
    if els is None:
        out += f"  #{int(cid):2d}  query failed\n"
        continue
    lu, bl, off = collections.Counter(), collections.Counter(), collections.Counter()
    named = []
    for e in els:
        t = e.get("tags", {})
        lu.update([t["landuse"]] if "landuse" in t else [])
        bl.update([t["building"]] if "building" in t else [])
        off.update([t["office"]] if "office" in t else [])
        if "name" in t and ("landuse" in t or "office" in t or "industrial" in t
                            or t.get("building") in ("industrial", "commercial", "office")):
            named.append(t["name"])
    ind_named = sum(any(w in nm.lower() for w in ("industr", "estate", "factory", "works", "plant", "sez"))
                    for nm in named)
    out += (f"  #{int(cid):2d} {int(n_):4d} nodes {km:4.1f} km  ({la:.4f}, {lo:.4f})  -> {verdict(lu, bl, ind_named)}\n"
            f"       landuse {dict(lu.most_common(6))}\n"
            f"       buildings {dict(bl.most_common(5))}  offices {dict(off)}\n"
            f"       named: {', '.join(list(dict.fromkeys(named))[:8]) or '-'}\n")
    time.sleep(2)                    # be polite to the public Overpass server
print(out)
with open(f"results/{tag}_cluster_landuse.txt", "w", encoding="utf-8") as f:
    f.write(out)
log(f"wrote results/{tag}_cluster_landuse.txt")
