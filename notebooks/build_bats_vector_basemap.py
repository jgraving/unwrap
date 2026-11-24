"""Build the bat figure's pinned vector basemap material.

Downloads the Geofabrik israel-and-palestine free-shp extract, clips
water/roads/buildings to the composite map extents, and writes
notebooks/assets/goldshtein_bats_figure/basemap_osm.geojson.gz.  The
renderer draws that extract as the vector basemap; update its SHA-256
in goldshtein_bats_figure_core.ASSET_SHA256 after re-running.

Data © OpenStreetMap contributors (ODbL), via Geofabrik.

Usage (requires network):
    curl -LO https://download.geofabrik.de/asia/israel-and-palestine/israel-and-palestine-latest-free.shp.zip
    python notebooks/build_bats_vector_basemap.py   # expects israel-and-palestine-latest-free.shp.zip in CWD
"""
"""Clip the Geofabrik free-shp extract to the figure's map extent and write
the pinned vector-basemap material (buildings/roads/water) as GeoJSON gz."""
import gzip
import json
from pathlib import Path

import geopandas as gpd
from pyproj import Transformer

SRC = Path("ishp")  # unzipped Geofabrik free-shp extract
OUT = (Path(__file__).resolve().parent / "assets" / "goldshtein_bats_figure"
     / "basemap_osm.geojson.gz")

# map extents in EPSG:32636 (composite a/d maps and notebook maps) with margin
transform = Transformer.from_crs(32636, 4326, always_xy=True)
xs = (738300 - 200, 746300 + 200)
ys = (3663600 - 200, 3672700 + 200)
corners = [transform.transform(x, y) for x in xs for y in ys]
west = min(c[0] for c in corners); east = max(c[0] for c in corners)
south = min(c[1] for c in corners); north = max(c[1] for c in corners)
bbox = (west, south, east, north)
print("bbox deg:", [round(v, 5) for v in bbox])

features = []

buildings = gpd.read_file(SRC / "gis_osm_buildings_a_free_1.shp", bbox=bbox)
buildings = buildings[~buildings.geometry.is_empty]
print("buildings:", len(buildings), "| classes:", buildings.get("fclass", buildings.columns[:3]).value_counts().head(6).to_dict() if "fclass" in buildings else "-")
for _, row in buildings.iterrows():
    geom = row.geometry
    geoms = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for g in geoms:
        features.append({"type": "Feature", "properties": {"kind": "building"},
                         "geometry": g.__geo_interface__})

roads = gpd.read_file(SRC / "gis_osm_roads_free_1.shp", bbox=bbox)
keep = {"primary", "secondary", "tertiary", "secondary_link", "tertiary_link",
        "unclassified", "residential", "road", "living_street"}
roads = roads[roads["fclass"].isin(keep)]
print("roads:", len(roads), roads["fclass"].value_counts().to_dict())
for _, row in roads.iterrows():
    features.append({"type": "Feature",
                     "properties": {"kind": "road", "fclass": row["fclass"]},
                     "geometry": row.geometry.__geo_interface__})

water = gpd.read_file(SRC / "gis_osm_water_a_free_1.shp", bbox=bbox)
water = water[~water.geometry.is_empty]
print("water:", len(water), water["fclass"].value_counts().head(6).to_dict() if "fclass" in water else "-")
for _, row in water.iterrows():
    features.append({"type": "Feature", "properties": {"kind": "water"},
                     "geometry": row.geometry.__geo_interface__})

fc = {"type": "FeatureCollection",
      "name": "osm_vector_basemap",
      "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.1:CRS:84"}},
      "bbox": list(bbox),
      "features": features}
OUT.parent.mkdir(parents=True, exist_ok=True)
with gzip.open(OUT, "wt", encoding="utf-8") as fh:
    json.dump(fc, fh, ensure_ascii=False, separators=(",", ":"))
print("wrote", OUT, OUT.stat().st_size, "bytes")
