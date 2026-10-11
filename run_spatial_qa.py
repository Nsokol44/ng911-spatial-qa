#!/usr/bin/env python3
"""NG911 spatial QA pipeline — five spatial checks on a county's own E-911 layers.

Addresses: county site-address points (26,792 in the reference run).
Centerlines: county-published roads with address ranges + ESN L/R fields
(5,365 segments). Boundary: Census cartographic boundary shared with the
neighboring county. The reference run's county is anonymized (an East
Tennessee county); its data is not distributed with this code.

Published-version note: input paths were generalized for publication. The
pipeline reads data/addresses.geojson, data/roads.geojson, and
data/county_boundary.zip (the Census cartographic boundary file) from a
data/ folder beside this script — substitute your own layers, and set
COUNTY_A / COUNTY_B below to your county pair's FIPS GEOIDs.
"""
import csv
import json
import math
import os
import re
import zipfile
import io
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import shapefile
from shapely.geometry import LineString, MultiLineString, Point, shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "figures")
RESULTS = os.path.join(HERE, "results")
os.makedirs(FIG, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

# FIPS GEOIDs (state + county FIPS) of the study county and the neighboring
# county that shares the seam under test. Generalized for publication —
# replace with your own county pair before running the seam section.
COUNTY_A = "00000"
COUNTY_B = "00001"

LAT0 = 36.02
KX = 111320 * math.cos(math.radians(LAT0))
KY = 110540.0


def to_m(lon, lat):
    return lon * KX, lat * KY


SUFFIX = {"STREET": "ST", "AVENUE": "AVE", "ROAD": "RD", "DRIVE": "DR",
          "LANE": "LN", "COURT": "CT", "CIRCLE": "CIR", "HIGHWAY": "HWY",
          "PARKWAY": "PKWY", "PLACE": "PL", "TRAIL": "TRL", "BOULEVARD": "BLVD"}


def norm_name(s):
    s = re.sub(r"\s+", " ", re.sub(r"[.,']", "", (s or "").upper().strip()))
    return " ".join(SUFFIX.get(t, t) for t in s.split(" "))


def street_name(predir, name, typ, sufdir):
    return norm_name(" ".join(x for x in [predir or "", name or "", typ or "", sufdir or ""] if x))


# ---------- county centerlines ----------
roads_gj = json.load(open(os.path.join(DATA, "roads.geojson")))
segments = []
for f in roads_gj["features"]:
    p = f["properties"]
    g = f.get("geometry")
    if not g:
        continue
    lines = []
    if g["type"] == "LineString":
        lines = [g["coordinates"]]
    elif g["type"] == "MultiLineString":
        lines = g["coordinates"]
    for coords in lines:
        pts = [to_m(x, y) for x, y in coords]
        if len(pts) >= 2:
            segments.append({
                "geom": LineString(pts),
                "name": street_name(p.get("PREFIX"), p.get("NAME"), p.get("TYPE"), None),
                "lo": min([v for v in (p.get("L_F_ADD"), p.get("R_F_ADD")) if v] or [0]),
                "hi": max([v for v in (p.get("L_T_ADD"), p.get("R_T_ADD")) if v] or [0]),
            })
named = [s for s in segments if s["name"]]
seg_geoms = [s["geom"] for s in segments]
named_geoms = [s["geom"] for s in named]
named_names = [s["name"] for s in named]
all_tree = STRtree(seg_geoms)
named_tree = STRtree(named_geoms)

# ---------- county address points ----------
gj = json.load(open(os.path.join(DATA, "addresses.geojson")))
addresses, skipped = [], 0
for f in gj["features"]:
    p = f["properties"]
    g = f.get("geometry")
    if not g or not g.get("coordinates"):
        skipped += 1
        continue
    street = street_name(p.get("PREDIR"), p.get("NAME"), p.get("TYPE"), p.get("SUFDIR"))
    num = (p.get("STNUM") or "").strip()
    if not street or not num:
        skipped += 1
        continue
    x, y = to_m(*g["coordinates"])
    m = re.match(r"\d+", num)
    addresses.append({"pt": Point(x, y), "street": street, "num": num,
                      "numval": int(m.group()) if m else None,
                      "esn": (p.get("ESN") or "").strip()})
N = len(addresses)
print(f"segments={len(segments)} named={len(named)} addresses={N} skipped={skipped}")


def run_checks(addrs, parity_r=45, outlier_d=60, do_range=True):
    parity, outliers, range_fail = [], [], []
    for i, rec in enumerate(addrs):
        d_any = rec["pt"].distance(seg_geoms[all_tree.nearest(rec["pt"])])
        if d_any > outlier_d:
            outliers.append((i, d_any))
            continue
        near = np.atleast_1d(named_tree.query(rec["pt"].buffer(parity_r), predicate="intersects"))
        cand = [k for k in near if named_names[k] == rec["street"]]
        if not cand:
            parity.append(i)
        elif do_range and rec["numval"] is not None:
            ok = False
            for k in cand:
                s = named[k]
                if s["lo"] <= rec["numval"] <= s["hi"]:
                    ok = True
                    break
            if not ok:
                range_fail.append(i)
    by_key = defaultdict(list)
    for i, rec in enumerate(addrs):
        by_key[(rec["street"], rec["num"])].append(i)
    dup_groups = []
    for key, idxs in by_key.items():
        if len(idxs) > 1:
            base = addrs[idxs[0]]["pt"]
            close = [j for j in idxs[1:] if base.distance(addrs[j]["pt"]) < 5]
            if close:
                dup_groups.append((key, [idxs[0]] + close))
    return parity, outliers, dup_groups, range_fail


P, O, D, RF = run_checks(addresses)
dup_extra = sum(len(g) - 1 for _, g in D)

sens = {}
for r_ in (35, 45, 55):
    sens[f"parity_r{r_}"] = len(run_checks(addresses, parity_r=r_, do_range=False)[0])
for t in (60, 90, 120):
    sens[f"outlier_t{t}"] = len(run_checks(addresses, outlier_d=t, do_range=False)[1])

# ---------- seam: real county line vs neighboring county ----------
with zipfile.ZipFile(os.path.join(DATA, "county_boundary.zip")) as z:
    base = [n[:-4] for n in z.namelist() if n.endswith(".shp")][0]
    rdr = shapefile.Reader(shp=io.BytesIO(z.read(base + ".shp")),
                           dbf=io.BytesIO(z.read(base + ".dbf")))
    fields = [f[0] for f in rdr.fields[1:]]
    gi = fields.index("GEOID")
    county = {}
    for rec, shp in zip(rdr.iterRecords(), rdr.iterShapes()):
        if rec[gi] in (COUNTY_A, COUNTY_B):
            county[rec[gi]] = shape(shp.__geo_interface__)


def proj(geom):
    from shapely.geometry import Polygon, MultiPolygon
    if geom.geom_type == "Polygon":
        return Polygon([(x * KX, y * KY) for x, y in geom.exterior.coords],
                       [[(x * KX, y * KY) for x, y in r.coords] for r in geom.interiors])
    return MultiPolygon([proj(g) for g in geom.geoms])


boundary = proj(county[COUNTY_A]).boundary.intersection(proj(county[COUNTY_B]).boundary)
seam_exposed = [i for i, rec in enumerate(addresses) if rec["pt"].distance(boundary) <= 50]

# ---------- repair + re-run ----------
fixed = [dict(r) for r in addresses]
for i in P:
    near = np.atleast_1d(named_tree.query(fixed[i]["pt"].buffer(45), predicate="intersects"))
    if len(near):
        j = min(near, key=lambda k: fixed[i]["pt"].distance(named_geoms[k]))
        fixed[i]["street"] = named_names[j]
for i, d_any in O:
    if d_any <= 150:
        j = all_tree.nearest(fixed[i]["pt"])
        fixed[i]["pt"] = nearest_points(seg_geoms[j], fixed[i]["pt"])[0]
drop = set()
for _, grp in D:
    drop.update(grp[1:])
kept = [r for i, r in enumerate(fixed) if i not in drop]
P2, O2, D2, RF2 = run_checks(kept)

stats = dict(
    n_addresses=N, n_segments=len(segments), n_named_segments=len(named),
    parity=len(P), parity_pct=round(100 * len(P) / N, 2),
    dup_groups=len(D), dup_extra=dup_extra,
    outliers=len(O), outlier_pct=round(100 * len(O) / N, 2),
    range_fail=len(RF), range_pct=round(100 * len(RF) / N, 2),
    sensitivity=sens,
    seam_within_50m=len(seam_exposed),
    boundary_km=round(boundary.length / 1000, 1),
    after=dict(parity=len(P2), outliers=len(O2),
               dup_extra=sum(len(g) - 1 for _, g in D2), range_fail=len(RF2), n=len(kept)),
    top_parity_streets=Counter(addresses[i]["street"] for i in P).most_common(8),
)
stats["total_exceptions"] = len(P) + len(O) + dup_extra + len(RF)
stats["exception_pct"] = round(100 * stats["total_exceptions"] / N, 2)
aft = stats["after"]
stats["after_total"] = aft["parity"] + aft["outliers"] + aft["dup_extra"] + aft["range_fail"]
stats["after_pct"] = round(100 * stats["after_total"] / aft["n"], 2)
json.dump(stats, open(os.path.join(RESULTS, "stats.json"), "w"), indent=2, default=str)
print(json.dumps(stats, indent=2, default=str))

# ---------- ranked exception report (CSV) ----------
# One row per flagged record, grouped by class, classes ranked by flag count.
# distance_m = distance to nearest centerline (outliers) or to the county
# boundary (seam queue); blank where distance is not the flagging measure.
seam_queue = [i for i, rec in enumerate(addresses) if rec["pt"].distance(boundary) <= 250]
rows = []
for i in P:
    rows.append(("name_mismatch", addresses[i]["street"], addresses[i]["num"], ""))
for i in RF:
    rows.append(("range_parity", addresses[i]["street"], addresses[i]["num"], ""))
for i, d_any in O:
    rows.append(("outlier", addresses[i]["street"], addresses[i]["num"], round(d_any, 1)))
for _, grp in D:
    for j in grp[1:]:
        rows.append(("duplicate", addresses[j]["street"], addresses[j]["num"], ""))
for i in seam_queue:
    rows.append(("seam_queue", addresses[i]["street"], addresses[i]["num"],
                 round(addresses[i]["pt"].distance(boundary), 1)))
class_counts = Counter(r[0] for r in rows)
rank = sorted(class_counts, key=lambda c: -class_counts[c])
rows.sort(key=lambda r: rank.index(r[0]))
with open(os.path.join(RESULTS, "exception_report.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["class", "street", "house_number", "distance_m"])
    w.writerows(rows)
print("ranked exception report:", dict(class_counts))

# ---------- figures ----------
NAVY, TEAL, RED = "#12304f", "#0e7c7b", "#c0392b"
CAP = "Real county E-911 layers, anonymized excerpt · county addresses vs county centerlines"
pts = np.array([(r["pt"].x, r["pt"].y) for r in addresses])
EXT = (pts[:, 0].min() - 400, pts[:, 0].max() + 400, pts[:, 1].min() - 400, pts[:, 1].max() + 400)


def draw_ctx(ax, boundary_too=False):
    for s in segments:
        x, y = s["geom"].xy
        ax.plot(x, y, color="#b9c5cd", lw=0.7, zorder=1)
    if boundary_too and not boundary.is_empty:
        for g in (boundary.geoms if hasattr(boundary, "geoms") else [boundary]):
            x, y = g.xy
            ax.plot(x, y, color=NAVY, lw=1.8, ls=(0, (5, 3)), zorder=5)
    ax.set_aspect("equal"); ax.axis("off")


pf = np.array([addresses[i]["pt"].coords[0] for i in P])
ol = np.array([addresses[i]["pt"].coords[0] for i, _ in O])

fig, ax = plt.subplots(figsize=(10.5, 7.4), dpi=200)
draw_ctx(ax)
ax.scatter(pts[:, 0], pts[:, 1], s=4, color=TEAL, alpha=0.5, zorder=3, label="Address points (pass)")
ax.scatter(pf[:, 0], pf[:, 1], s=26, color=RED, zorder=6, edgecolor="white", lw=0.5,
           label="Parity flags — no same-name street within 45 m")
ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
ax.legend(loc="upper left", fontsize=15, framealpha=1.0, markerscale=1.6).set_zorder(20)
fig.tight_layout(); fig.savefig(f"{FIG}/fig1_parity.png"); plt.close(fig)

cx, cy = (np.median(ol[:, 0]), np.median(ol[:, 1])) if len(ol) else (np.median(pts[:, 0]), np.median(pts[:, 1]))
fig, ax = plt.subplots(figsize=(10.5, 7.4), dpi=200)
draw_ctx(ax)
m = (abs(pts[:, 0] - cx) < 2600) & (abs(pts[:, 1] - cy) < 2600)
ax.scatter(pts[m, 0], pts[m, 1], s=9, color=TEAL, alpha=0.5, zorder=3, label="Address points")
ax.scatter(ol[:, 0], ol[:, 1], s=95, facecolor="none", edgecolor=RED, lw=1.9, zorder=6,
           label="Outliers (>60 m from nearest centerline)")
shown = 0
for _, grp in D:
    ps = [addresses[j]["pt"] for j in grp]
    mx = np.mean([p.x for p in ps]); my = np.mean([p.y for p in ps])
    if abs(mx - cx) < 2600 and abs(my - cy) < 2600:
        ax.add_patch(plt.Circle((mx, my), 95, fill=False, ec="#8e44ad", lw=1.8, zorder=6)); shown += 1
ax.plot([], [], color="#8e44ad", lw=2.2, label=f"Duplicate clusters in view ({shown})")
ax.set_xlim(cx - 2600, cx + 2600); ax.set_ylim(cy - 2600, cy + 2600)
ax.legend(loc="upper left", fontsize=15, framealpha=1.0, markerscale=1.4).set_zorder(20)
fig.tight_layout(); fig.savefig(f"{FIG}/fig2_dup_out.png"); plt.close(fig)

seam_wide = [i for i, rec in enumerate(addresses) if rec["pt"].distance(boundary) <= 250]
sp_all = np.array([addresses[i]["pt"].coords[0] for i in seam_wide])
if len(sp_all):
    dens = [np.sum(np.abs(sp_all[:, 0] - x) + np.abs(sp_all[:, 1] - y) < 2000) for x, y in sp_all]
    bx_, by_ = sp_all[int(np.argmax(dens))]
else:
    bx_, by_ = cx, cy
fig, ax = plt.subplots(figsize=(10.5, 7.4), dpi=200)
draw_ctx(ax, boundary_too=True)
m = (abs(pts[:, 0] - bx_) < 2400) & (abs(pts[:, 1] - by_) < 2400)
ax.scatter(pts[m, 0], pts[m, 1], s=10, color=TEAL, alpha=0.55, zorder=3, label="Address points")
if len(sp_all):
    inm = (abs(sp_all[:, 0] - bx_) < 2400) & (abs(sp_all[:, 1] - by_) < 2400)
    ax.scatter(sp_all[inm, 0], sp_all[inm, 1], s=70, color="#e67e22", zorder=6, edgecolor="white", lw=0.7,
               label="Within 250 m of the county line")
ax.set_xlim(bx_ - 2400, bx_ + 2400); ax.set_ylim(by_ - 2400, by_ + 2400)
ax.legend(loc="upper left", fontsize=15, framealpha=1.0, markerscale=1.4).set_zorder(20)
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_seam.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(10.6, 5.9), dpi=200)
cats = ["Parity", "Out-of-range\nnumbers", "Duplicate\nextras", "Outliers"]
before = [len(P), len(RF), dup_extra, len(O)]
after = [aft["parity"], aft["range_fail"], aft["dup_extra"], aft["outliers"]]
x = np.arange(4); w = 0.36
b1 = ax.bar(x - w / 2, before, w, color=RED, label="First run")
b2 = ax.bar(x + w / 2, after, w, color=TEAL, label="After scripted fixes + re-run")
for bars in (b1, b2):
    for r in bars:
        ax.annotate(f"{int(r.get_height()):,}", (r.get_x() + r.get_width() / 2, r.get_height()),
                    ha="center", va="bottom", fontsize=14, weight="bold", color=NAVY)
ax.set_xticks(x, cats, fontsize=14.5)
ax.set_ylabel("Flagged records", fontsize=14)
ax.tick_params(axis="y", labelsize=12.5)
ax.set_ylim(0, 6400)
ax.legend(fontsize=13.5, loc="upper left"); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_before_after.png"); plt.close(fig)
print("figures done")
