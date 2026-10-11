# NG911 Spatial QA Pipeline

Five spatial quality checks for NG911 address data that catch the errors
schema validation misses. Every record in a typical county E-911 layer
already passes field-level (NENA-style) validation; these checks ask the
spatial questions — *is the point on the street it claims, in the number
range it claims, near any road at all, unique, and consistent across the
county line?* — that decide whether a call routes and geocodes correctly.

Presented at the TNGIC East Tennessee Fall Forum 2026 (Nicholas J. Sokol,
NovelGeo Technologies): *Beyond Schema Checks: Catching the NG911 Errors
That Pass Validation but Misroute Calls.*

## The checklist

The pipeline runs five checks over a county's site-address points against
its road centerlines. Classes are mutually exclusive by construction: the
outlier test runs first, so a point too far from any centerline is never
also counted as a name or range failure.

1. **Street-name agreement** — the address point's street name is compared
   against the *candidate set* of named centerlines within a **45 m**
   window, not just the single nearest segment. Names on both sides are
   normalized first with a suffix alias table (STREET→ST, AVENUE→AVE,
   ROAD→RD, …) and pre/post-directionals, so the comparison is between
   like forms. A point with no same-name centerline in its window is a
   name mismatch.
2. **Address-range parity** — for points that pass check 1, the house
   number must fall inside a candidate segment's published address range
   (left/right from/to). Out-of-range numbers flag even when the street
   name is right.
3. **Spatial outliers** — a point more than **60 m** from the nearest
   centerline of any name is an outlier: geocoded to the wrong place,
   snapped to a driveway or parcel centroid, or a genuinely missing road.
4. **Duplicates** — records sharing the same normalized street + house
   number and sitting within 5 m of each other are duplicate extras.
5. **County-seam queue** — address points within **250 m** of the shared
   county boundary are queued for cross-county review (the stats file
   also reports the innermost 50 m band). Seam cases are a review queue,
   not automatic errors: two counties' independently maintained layers
   meet at the line, and only the queue tells you where to look.

Threshold sensitivity is built in: the run re-counts name mismatches at
35/45/55 m windows and outliers at 60/90/120 m, so you can see how much of
the headline rate is threshold choice and how much is data.

## Method notes

- **Distances are ground metres under a local equirectangular scaling**
  anchored at 36.02° N: `x = lon · 111320 · cos(36.02°)`,
  `y = lat · 110540`. This is *not* Web Mercator (EPSG:3857) — at this
  latitude Mercator inflates distances by ~24% and would silently shrink
  every threshold. Over a single county the local scaling's distortion is
  far below the thresholds' sensitivity bands. If you adapt the pipeline
  to another latitude, re-anchor `LAT0`.
- **Name comparison uses a candidate set**, never the nearest segment
  alone — at intersections and in grid street patterns the nearest
  centerline is often the cross street, and a nearest-only comparison
  manufactures false mismatches.
- **Repair pass is reported, not hidden.** The script applies scripted
  candidate fixes (name alignment to the nearest candidate centerline,
  snapping outliers ≤ 150 m, dropping duplicate extras), re-runs the
  checks, and reports counts *before and after*, side by side.
- Everything is deterministic and reproducible: fixed inputs in, same
  exception report out.

## Inputs

Bring your own layers — **the reference run's county data is not
included** (the talk anonymizes the county, and this repo does the same).
Place these in a `data/` folder beside the script:

| File | Contents |
|---|---|
| `data/addresses.geojson` | Site-address points (GeoJSON, WGS84 lon/lat) with properties `PREDIR`, `NAME`, `TYPE`, `SUFDIR`, `STNUM`, `ESN` |
| `data/roads.geojson` | Road centerlines (GeoJSON) with properties `PREFIX`, `NAME`, `TYPE`, `L_F_ADD`, `R_F_ADD`, `L_T_ADD`, `R_T_ADD` |
| `data/county_boundary.zip` | Census cartographic boundary shapefile (zipped) covering your county and its neighbor |

Then set `COUNTY_A` and `COUNTY_B` at the top of `run_spatial_qa.py` to
the two counties' FIPS GEOIDs (state + county FIPS, as in the Census
file's `GEOID` field). Field names follow a typical county E-911 /
NG911 GIS layer export; rename your fields to match, or adjust the
property lookups — they are the only schema assumptions in the code.

Published-version note: input filenames and the county GEOIDs were
generalized for publication; no check logic or threshold differs from
the reference run.

## Run

```bash
pip install -r requirements.txt
python3 run_spatial_qa.py
```

One command, end to end: load → five checks → sensitivity recount →
boundary/seam analysis → repair + re-run → outputs.

## Outputs

- `results/exception_report.csv` — the ranked exception report: one row
  per flagged record (`class, street, house_number, distance_m`),
  grouped by class, classes ranked by flag count. `distance_m` is the
  distance to the nearest centerline for outliers and to the county
  boundary for the seam queue; blank where distance isn't the flagging
  measure.
- `results/stats.json` — full counts before/after the repair pass, the
  threshold-sensitivity table, and seam exposure.
- `figures/` — parity, outlier/duplicate, seam, and before/after charts.

## Reference run (TNGIC 2026)

One East Tennessee county's own layers: 26,792 address points against
5,365 centerline segments. First pass flagged **9,302 records (34.7%)**
across the name-mismatch, outlier, duplicate, and range classes —
before the separate 250 m seam queue. Counts by class and the
sensitivity table are in
[`results/real-stats.json`](results/real-stats.json). (Street-level
detail from that run is withheld, consistent with the county's
anonymization.)

## License

MIT — © 2026 Nicholas J. Sokol. See [LICENSE](LICENSE).
