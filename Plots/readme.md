# Plot guide: Chennai and Bengaluru road networks

How to read every plot the scripts produce, and the idea behind it.
Sections marked **Status: run** have results from your runs. Sections marked
**Status: not run yet** describe what to look for, not results.

## Conventions used in all plots

| Item | Meaning |
|---|---|
| Orange | Chennai |
| Green | Bengaluru |
| Blue | Erdős–Rényi (ER) random network with the same number of junctions and roads |
| Pink (solid) / purple (dashed) | Scale-free networks, γ = 2.5 and γ = 3.5 |
| Node | A junction or the end of a road |
| Edge | A road segment between two junctions |
| Edge-only maps | Every line is a road segment. No dots are drawn for junctions. |

Words used in the captions:
- **Giant component:** the largest group of junctions that can still reach each other. Think of it as the part of the city you can still drive through.
- **Hop:** one road segment. "Distance in hops" = number of segments on the route.

---

## 1. Robustness curves (`road_robustness_updated.py` → `robustness_plots/`)

**Status: run**

### `Chennai_robustness.png`, `Bengaluru_robustness.png`

**What it shows.** Three panels, each removing junctions in a different way:
random failure, targeted attack by betweenness, targeted attack by degree.
Each panel compares the city (thick line) with ER and with scale-free networks.

**How to read it.**
- x-axis: fraction of junctions removed, *f* (0 = nothing removed, 1 = everything).
- y-axis: size of the giant component compared with its starting size, P∞(f)/P∞(0). 1 means the city is still fully connected, and 0 means it has fallen apart.
- A curve that stays high for longer is more robust. The point where it reaches zero is the critical fraction f_c. The area under the curve is R, and a bigger area is more robust.

**Intuition.**
- *Random failure* models something like a flood or accident that hits anywhere.
- *Betweenness attack* removes the junctions that lie on the most shortest routes first (bridges, arterials, key crossings).
- *Degree attack* removes the junctions with the most roads first.
- Comparing the two attacks tells you *why* a network is fragile: if the betweenness attack is much worse than the degree attack, importance comes from position, not from the number of roads meeting at a junction.

**What we saw** (read from the plots, approximate).
- Random: roads collapse near f ≈ 0.3, about half of ER's (≈ 0.6).
- Betweenness attack: roads collapse near f ≈ 0.07, about the same as scale-free γ = 2.5, and far earlier than ER (≈ 0.35).
- Degree attack: roads collapse near f ≈ 0.2, close to ER. So road networks are vulnerable through their corridors, not through hubs.
- Scale-free with a lower γ survives random failure better and breaks earlier under attack.

**Why the road curves are jagged.** The steps are real. Removing one corridor junction can cut off a whole region at once. ER and scale-free curves are smooth because removing one node removes only a tiny piece. The random curve is an average over several runs, so it is smoother.

### `chennai_vs_bengaluru.png`

**What it shows.** The same three panels with both cities on one set of axes. Thick lines are the roads. Dashed lines are the scale-free networks.

**How to read it.** Look only at the thick orange and green lines and ask whether one stays above the other.

**What we saw.** Random failure: the curves overlap (tie). Betweenness attack: Bengaluru holds a little longer (around f ≈ 0.02–0.06). Degree attack: they nearly overlap.

---

## 2. Path structure (`road_path_structure.py` → `path_structure_plots/`)

**Status: run**

### `path_length_distributions.png`

**What it shows.** For each city, a histogram of the shortest-path length (in hops) between pairs of junctions. Orange is the road network, blue is the ER network of the same size.

**How to read it.**
- x-axis: number of hops. y-axis: fraction of junction pairs at that distance. Both shapes have the same total area (1), so the low orange bars only mean it is spread over many values.
- A tall narrow spike means almost every pair is at about the same distance. A wide flat shape means distances vary a lot.

**Intuition.** In a random network, links can join any two places, so everything is a few hops from everything. On roads, a link can only join nearby places, so the hop count between two points follows the physical distance between them.

**What we saw.** Chennai: roads average 107 hops, ER 11.2. Bengaluru: roads 101.7, ER 12.1. Roads are about 9× longer. Chennai's distribution is flatter with a longer tail (matches its long, narrow shape). Bengaluru's is a more symmetric bell near 100.

**Caution.** The averages come from 200 random start points, each measured to every junction. The diameter is a lower bound. Hop counts depend on mapping detail, so compare the two cities in km, not only in hops.

### `Chennai_hop_maps.png`, `Bengaluru_hop_maps.png`

**What it shows.** Two iterations. Each picks a random source junction (star) and colors every road by its hop distance from it. The red line is the longest shortest-path, ending at the farthest junction (X).

**How to read it.** Dark = close to the source, yellow = far. Look for smooth rings of color spreading out from the star.

**Intuition.** The rings are the signature of a spatial network: distance grows steadily as you move away. In a random network there would be almost no rings, because after a few hops you could be anywhere.

**What we saw.** The farthest points sit at the edges of the city (for example Chennai's long north–south extent reaches 236 hops from one source). The colors change smoothly with geography.

---

## 3. Spatial analysis (`road_spatial_analysis.py` → `spatial_plots/`)

**Status: run**

### `Chennai_spatial.png`, `Bengaluru_spatial.png`

Four panels per city.

**Top left: street-orientation map.** Each road is colored by its compass direction (0–180°).
- Reading: a planned grid shows large patches of only one or two colors. A mixed, organic layout shows a blend of all colors.
- Intuition: the map shows where streets align and where they don't.
- Seen: both cities are a thick mix, so both are organic. The long straight lines in other colors are main roads crossing the older pattern.

**Top right: polar histogram.** A bar chart wrapped around a compass. Each bar's length is the share of total road length pointing in that direction. Every street is counted in both directions, so the plot is symmetric.
- Reading: long bars at 0°/180° and 90°/270° mean most roads run north-south and east-west. The title gives the orientation entropy H.
- Intuition: H measures how spread out the directions are. A perfect grid gives about 1.4, and fully mixed directions give 3.58.
- Seen: H = 3.16 (Chennai) and 3.28 (Bengaluru). Both are close to fully mixed. Chennai is slightly more aligned along the main axes.

**Bottom: two circuity iterations.** Each picks a random origin (star) and destination (X) about 5 km apart in a straight line. Dashed = straight line, red = the real shortest road route.
- Reading: the title gives the circuity = route length ÷ straight-line length.
- Intuition: it shows how much longer a real trip is than the straight line. 1.0 is a straight route, and 1.3 is a 30% detour.
- Seen: Bengaluru 1.18 and 1.19. Chennai 1.55 and 1.29. Two examples only illustrate the idea.

### `circuity_comparison.png`

**What it shows.** Histogram of circuity for about 300 random origin-destination pairs per city.

**How to read it.** x-axis: circuity. y-axis: density (each city's bars scaled to the same area). A peak close to 1.0 means direct routes. A long tail to the right means some trips need big detours.

**What we saw.** Both peak between 1.15 and 1.30. Chennai mean 1.28 and median 1.24. Bengaluru mean 1.22 and median 1.20. Chennai has a few very circuitous pairs (up to about 3). The bulk of the two distributions overlap, so the gap is small.

---



## 4. Growth over time (`road_growth_over_time.py` → `growth_plots/`)

**Status: not run yet** (the historical downloads need to finish first). Remember that part of any "growth" may be better mapping, not construction.



### `{city}_1_robustness_over_time.png`

- Left: random-failure curves, one per year. Middle: betweenness-attack curves. Right: R against year.
- Reading: curves that move up over the years mean the network becomes more robust. Curves that stay on top of each other mean no change.

### `{city}_2_new_roads.png`

- Left: map with new roads in red over the old network in grey. Right: bars of new road length ÷ old road length by distance from the city centre.
- Reading: red in the outskirts means expansion. Red inside the grey means infill. A taller bar at a larger distance means growth is concentrated at the edge.



### `cities_over_time.png`

- Robustness R over time for both cities and the attachment exponent α for each.

---


