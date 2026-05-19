"""
utils.py
========
Spatial utility functions for the malaria spatiotemporal model.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import shape


def build_adjacency_matrix_from_geojson(
    geojson_path: str | Path,
    label: Literal["id", "displayName"] | str = "id",
    output_path: str | Path | None = None,
) -> pd.DataFrame:
    """Parse a GeoJSON FeatureCollection and compute a binary spatial adjacency
    matrix suitable for use as *W* in the CARBayesST spatiotemporal model.

    Two areas are considered adjacent (W[i,j] = 1) if their geometries share
    at least a border segment or vertex (queen contiguity).  The diagonal is
    always zero and the matrix is guaranteed symmetric.

    Parameters
    ----------
    geojson_path:
        Path to the GeoJSON FeatureCollection file.
    label:
        Which value to use as the row/column identifier.

        * ``"id"``          – top-level ``feature["id"]`` (DHIS2 orgunit UID,
                              e.g. ``"V42cbWvTYhu"``).  This matches what
                              chap-core writes into the ``location`` column.
        * ``"displayName"`` – ``feature["properties"]["displayName"]``.  Use
                              this when the training CSV was built from a
                              pipeline that stored human-readable names in
                              ``orgunitname`` (e.g. the legacy uppercase names).
        * Any other string  – reads ``feature["properties"][label]``.

    output_path:
        If provided, the adjacency matrix is written to this path as a CSV.
        The CSV has a header row and an index column of area labels so R can
        read it directly::

            W_df <- read.csv("adjacency.csv", row.names = 1, check.names = FALSE)
            W    <- as.matrix(W_df)

    Returns
    -------
    pd.DataFrame
        Square binary adjacency matrix with shape ``(n_areas, n_areas)``.
        Row and column labels are the chosen orgunit identifiers.
    """
    geojson_path = Path(geojson_path)
    gdf = gpd.read_file(geojson_path)

    # ── Resolve labels ────────────────────────────────────────────────────────
    if label == "id":
        # GeoJSON top-level feature id lands in the 'id' column via read_file
        if "id" not in gdf.columns:
            # Fallback: re-parse to pull feature-level id
            with open(geojson_path) as fh:
                raw = json.load(fh)
            ids = [f.get("id") or f["properties"].get("id") for f in raw["features"]]
            gdf["_label"] = ids
        else:
            gdf["_label"] = gdf["id"].astype(str)
    else:
        prop = label  # e.g. "displayName"
        if prop not in gdf.columns:
            raise KeyError(
                f"Property '{prop}' not found in GeoJSON features. "
                f"Available columns: {list(gdf.columns)}"
            )
        gdf["_label"] = gdf[prop].astype(str)

    # Drop features without valid geometry or label
    gdf = gdf.dropna(subset=["geometry", "_label"]).reset_index(drop=True)

    labels = gdf["_label"].tolist()
    n = len(labels)

    # ── Build adjacency via spatial index ─────────────────────────────────────
    # Use a spatial index (STRtree) to avoid O(n²) geometry intersection tests.
    # Two areas are adjacent if their geometries touch or overlap after buffering
    # by a tiny epsilon to bridge floating-point gaps in shared borders.
    EPSILON = 1e-8  # degrees – negligible for any real coordinate system

    adj = np.zeros((n, n), dtype=np.int32)
    tree = gdf.sindex

    for i, geom_i in enumerate(gdf.geometry):
        if geom_i is None or geom_i.is_empty:
            continue
        # Buffer slightly so polygons that share an exact border are detected
        buffered = geom_i.buffer(EPSILON)
        candidate_idxs = list(tree.intersection(buffered.bounds))
        for j in candidate_idxs:
            if j == i:
                continue
            geom_j = gdf.geometry.iloc[j]
            if geom_j is None or geom_j.is_empty:
                continue
            # Adjacent if they share any part of their boundary
            if geom_i.touches(geom_j) or geom_i.intersects(geom_j):
                adj[i, j] = 1
                adj[j, i] = 1

    # Ensure diagonal is zero
    np.fill_diagonal(adj, 0)

    W = pd.DataFrame(adj, index=labels, columns=labels)

    if output_path is not None:
        W.to_csv(output_path)
        print(f"Wrote {n}×{n} adjacency matrix to {output_path}")

    return W


# ── CLI convenience ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Build a spatial adjacency matrix from a GeoJSON file."
    )
    parser.add_argument("geojson", help="Path to the input GeoJSON file.")
    parser.add_argument(
        "--label",
        default="id",
        help=(
            "Feature field to use as row/column label. "
            "'id' (default) uses the top-level feature id; "
            "'displayName' uses properties.displayName."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Path to write the CSV adjacency matrix (default: print summary only).",
    )
    args = parser.parse_args()

    W = build_adjacency_matrix_from_geojson(
        args.geojson, label=args.label, output_path=args.output
    )
    n_areas = len(W)
    n_edges = int(W.values.sum()) // 2
    avg_neighbors = W.values.sum(axis=1).mean()
    isolated = int((W.values.sum(axis=1) == 0).sum())
    print(f"Areas       : {n_areas}")
    print(f"Edges       : {n_edges}")
    print(f"Avg neighbors: {avg_neighbors:.2f}")
    print(f"Isolated    : {isolated}")
