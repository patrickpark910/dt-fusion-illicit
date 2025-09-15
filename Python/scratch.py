#!/usr/bin/env python3
"""
Tokamak blanket volume using Miller cross-sections and Pappus' centroid theorem.

Method
------
1) Generate a closed polygon of the Miller D-shape cross-section in (R,Z) using
   the provided `miller_points` (returns centimeters).
2) Ensure CCW orientation and compute the inner area's centroid and perimeter
   via the shoelace/Green's theorem formulas.
3) Compute unit outward normals along the polygon (discrete mitered average of
   forward/backward tangents) and offset each vertex by +t (in cm).
4) Recompute area and centroid for the offset polygon.
5) Toroidal volumes from Pappus: V = 2π * (R̄) * A, where R̄ is the centroid's
   R-coordinate (distance to the rotation axis). Blanket volume is V_out - V_in.
6) Report a Steiner–Minkowski approximation as a numerical sanity check.

All code comments are in function docstrings as requested.
"""

import argparse
import math
import numpy as np


def miller_points(R0_m: float, a_m: float, delta: float, kappa: float, n: int = 400) -> np.ndarray:
    """Return Nx2 array of (r,z) points [cm] for a Miller D-shape.
    R0_m : major radius [m]
    a_m  : minor radius [m]
    delta: triangularity delta
    kappa: elongation kappa
    n    : number of points around the boundary
    """
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    R = R0_m + a_m*np.cos(theta + delta*np.sin(theta))
    Z = kappa*a_m*np.sin(theta)
    pts_m = np.column_stack([R, Z])
    pts_cm = 100.0*pts_m
    # Close the loop (Polygon expects a closed path)
    return np.vstack([pts_cm, pts_cm[0]])


def _as_open_polygon(pts: np.ndarray) -> np.ndarray:
    """Return polygon as an open sequence of vertices (drop duplicated final vertex if present)."""
    pts = np.asarray(pts)
    if np.allclose(pts[0], pts[-1]):
        return pts[:-1].copy()
    return pts.copy()


def polygon_area_centroid(pts: np.ndarray):
    """Compute signed area and centroid (Cx,Cy) of a planar polygon (closed or open)."""
    P = _as_open_polygon(pts)
    x = P[:, 0]
    y = P[:, 1]
    x_next = np.roll(x, -1)
    y_next = np.roll(y, -1)
    cross = x * y_next - x_next * y
    A = 0.5 * float(np.sum(cross))
    Cx = (1.0 / (6.0 * A)) * float(np.sum((x + x_next) * cross))
    Cy = (1.0 / (6.0 * A)) * float(np.sum((y + y_next) * cross))
    return A, Cx, Cy


def polygon_perimeter(pts: np.ndarray) -> float:
    """Compute the perimeter length of a closed polygon (closed or open input accepted)."""
    P = _as_open_polygon(pts)
    diffs = np.diff(np.vstack([P, P[0]]), axis=0)
    return float(np.sum(np.hypot(diffs[:, 0], diffs[:, 1])))


def ensure_ccw(pts: np.ndarray) -> np.ndarray:
    """Ensure polygon is CCW-oriented; return as a closed polygon (first point repeated)."""
    P = _as_open_polygon(pts)
    A, _, _ = polygon_area_centroid(P)
    if A < 0.0:
        P = P[::-1]
    return np.vstack([P, P[0]])


def outward_normals_ccw(pts: np.ndarray) -> np.ndarray:
    """Compute unit outward normals for a CCW-oriented polygon using mitered tangents."""
    P = _as_open_polygon(ensure_ccw(pts))
    fwd = np.roll(P, -1, axis=0) - P
    bwd = P - np.roll(P, 1, axis=0)
    fwd_n = fwd / (np.linalg.norm(fwd, axis=1, keepdims=True) + 1e-16)
    bwd_n = bwd / (np.linalg.norm(bwd, axis=1, keepdims=True) + 1e-16)
    tan = fwd_n + bwd_n
    bad = (np.linalg.norm(tan, axis=1) < 1e-12)
    tan[bad] = fwd_n[bad]
    tan /= (np.linalg.norm(tan, axis=1, keepdims=True) + 1e-16)
    normals = np.column_stack([tan[:, 1], -tan[:, 0]])
    normals /= (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-16)
    return normals


def offset_polygon_outward(pts: np.ndarray, t_cm: float) -> np.ndarray:
    """Offset a CCW polygon outward by a normal distance t_cm; returns a closed polygon."""
    P_ccw = ensure_ccw(pts)
    P_open = _as_open_polygon(P_ccw)
    N = outward_normals_ccw(P_open)
    Q_open = P_open + t_cm * N
    return np.vstack([Q_open, Q_open[0]])


def blanket_volume(
    R0_m: float,
    a_m: float,
    delta: float,
    kappa: float,
    t_m: float,
    n: int = 1000,
):
    """Compute plasma, outer, and blanket toroidal volumes for a Miller D-shape.

    Parameters
    ----------
    R0_m : float
        Major radius to plasma axis [m].
    a_m : float
        Plasma minor radius [m].
    delta : float
        Triangularity δ.
    kappa : float
        Elongation κ.
    t_m : float
        Blanket thickness measured normal to plasma boundary [m].
    n : int
        Number of points along the boundary for discretization.

    Returns
    -------
    dict
        Keys:
          - A_in_cm2, C_in_cm (R̄,Z̄), P_in_cm
          - A_out_cm2, C_out_cm (R̄,Z̄)
          - V_plasma_m3, V_total_out_m3, V_blanket_m3
          - V_blanket_steiner_m3  (sanity check)
    """
    if t_m < 0:
        raise ValueError("Blanket thickness t_m must be non-negative.")

    pts_cm = miller_points(R0_m, a_m, delta, kappa, n)
    pts_cm = ensure_ccw(pts_cm)

    A_in, Cx_in, Cy_in = polygon_area_centroid(pts_cm)
    P_in = polygon_perimeter(pts_cm)

    t_cm = 100.0 * t_m
    pts_out_cm = offset_polygon_outward(pts_cm, t_cm)

    A_out, Cx_out, Cy_out = polygon_area_centroid(pts_out_cm)

    V_in_cm3 = 2.0 * math.pi * Cx_in * A_in
    V_out_cm3 = 2.0 * math.pi * Cx_out * A_out
    V_blanket_cm3 = V_out_cm3 - V_in_cm3

    A_delta = P_in * t_cm + math.pi * t_cm * t_cm
    V_steiner_cm3 = 2.0 * math.pi * Cx_in * A_delta

    to_m3 = 1.0e-6
    return {
        "A_in_cm2": A_in,
        "C_in_cm": (Cx_in, Cy_in),
        "P_in_cm": P_in,
        "A_out_cm2": A_out,
        "C_out_cm": (Cx_out, Cy_out),
        "V_plasma_m3": V_in_cm3 * to_m3,
        "V_total_out_m3": V_out_cm3 * to_m3,
        "V_blanket_m3": V_blanket_cm3 * to_m3,
        "V_blanket_steiner_m3": V_steiner_cm3 * to_m3,
    }


def main():
    """CLI for computing tokamak blanket volume from torus dimensions and thickness."""


    out = blanket_volume(4, 1.2, 0.5, 1.5, 1, n=20)

    print("\n=== Plasma cross-section (cm units) ===")
    print(f"Area A_in [cm^2]     : {out['A_in_cm2']:.6f}")
    print(f"Centroid R̄_in [cm]   : {out['C_in_cm'][0]:.6f}")
    print(f"Perimeter P_in [cm]  : {out['P_in_cm']:.6f}")

    print("\n=== Offset (outer) cross-section (cm units) ===")
    print(f"Area A_out [cm^2]    : {out['A_out_cm2']:.6f}")
    print(f"Centroid R̄_out [cm]  : {out['C_out_cm'][0]:.6f}")

    print("\n=== Volumes ===")
    print(f"Plasma volume V_in [m^3]        : {out['V_plasma_m3']:.6f}")
    print(f"Total (outer) volume V_out [m^3]: {out['V_total_out_m3']:.6f}")
    print(f"Blanket volume V_blanket [m^3]  : {out['V_blanket_m3']:.6f}")
    print(f"Steiner approx (check) [m^3]    : {out['V_blanket_steiner_m3']:.6f}")


if __name__ == "__main__":
    main()
