import numpy as np

def miller_model(R0, a, kappa, delta, extrude=0.0, calc_vol=True, n=100):
    """
    Generate an offset Miller contour and (optionally) compute total toroidal volume.

    Args:
        R0 (float): major radius (same units as a, extrude)
        a (float): minor radius
        kappa (float): elongation
        delta (float): triangularity
        extrude (float): outward offset thickness (e.g., for a blanket layer)
        calc_vol (bool): if True, return toroidal volume by Pappus
        n (int): number of parametric points around the contour

    Returns:
        If calc_vol:
            (R_off, Z_off, volume_total)  # arrays closed (first point repeated)
        Else:
            (R_off, Z_off)
    """
    t = np.linspace(0.0, 2.0*np.pi, n, endpoint=False)

    # Miller contour (Miller+98; Ball et al.)
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)

    # Derivatives wrt parameter t
    dR_dt = -a * np.sin(t + delta * np.sin(t)) * (1.0 + delta * np.cos(t))
    dZ_dt =  kappa * a * np.cos(t)

    # Outward normal from (dR, dZ): N = (dZ, -dR), then unit-normalize
    N_R = dZ_dt
    N_Z = -dR_dt
    N_mag = np.hypot(N_R, N_Z)
    N_Ru, N_Zu = N_R / N_mag, N_Z / N_mag

    # Offset curve
    R_off = R + extrude * N_Ru
    Z_off = Z + extrude * N_Zu

    # Close polygon (repeat first point)
    R_off = np.concatenate([R_off, R_off[:1]])
    Z_off = np.concatenate([Z_off, Z_off[:1]])

    if not calc_vol:
        return R_off, Z_off

    A, Cx, _ = _poly_area_centroid(R_off, Z_off)
    volume_total = 2.0 * np.pi * Cx * abs(A)

    return R_off, Z_off, volume_total


def miller_inboard_outboard_volumes(R0, a, kappa, delta, extrude=0.0, n=200):
    """
    Compute **inboard** (R <= R0) and **outboard** (R >= R0) toroidal volumes
    of an offset Miller boundary using Sutherland–Hodgman clipping against x=R0.

    Args:
        R0 (float): major radius (split plane is R=R0)
        a (float): minor radius
        kappa (float): elongation
        delta (float): triangularity
        extrude (float): outward offset thickness
        n (int): contour resolution (more points => higher geometric fidelity)

    Returns:
        dict with:
          - 'RZ_in': (R_in, Z_in) closed polygon (may be empty if degenerate)
          - 'RZ_out': (R_out, Z_out) closed polygon
          - 'V_in': inboard toroidal volume
          - 'V_out': outboard toroidal volume
          - 'V_total': V_in + V_out (≈ total from full shape)
    """
    R_off, Z_off, V_total_direct = miller_model(R0, a, kappa, delta, extrude=extrude, calc_vol=True, n=n)

    # Clip against half-planes R <= R0 (inboard) and R >= R0 (outboard)
    R_in, Z_in = _clip_vertical_halfplane(R_off, Z_off, x0=R0, keep='le')
    R_out, Z_out = _clip_vertical_halfplane(R_off, Z_off, x0=R0, keep='ge')

    # Areas and centroids
    A_in, Cx_in, _ = _poly_area_centroid(R_in, Z_in) if len(R_in) >= 3 else (0.0, 0.0, 0.0)
    A_out, Cx_out, _ = _poly_area_centroid(R_out, Z_out) if len(R_out) >= 3 else (0.0, 0.0, 0.0)

    # Pappus volumes (use |A| but keep centroid's R for the Pappus radius)
    V_in  = 2.0 * np.pi * Cx_in  * abs(A_in)  if abs(A_in)  > 0 else 0.0
    V_out = 2.0 * np.pi * Cx_out * abs(A_out) if abs(A_out) > 0 else 0.0

    # Minor numerical differences can occur; return both split sum and direct total
    return {
        'RZ_in':  (R_in, Z_in),
        'RZ_out': (R_out, Z_out),
        'V_in':   V_in,
        'V_out':  V_out,
        'V_total': V_in + V_out,
        'V_total_direct': V_total_direct
    }


def _poly_area_centroid(R, Z):
    """
    Signed area and centroid (Cx, Cz) of a closed polygon.
    Input may be open or closed; this function will close it if needed.

    Returns:
        (A, Cx, Cz) with A signed (CCW positive). Works for convex or concave polygons.
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    if (R[0] != R[-1]) or (Z[0] != Z[-1]):
        R = np.concatenate([R, R[:1]])
        Z = np.concatenate([Z, Z[:1]])

    x, y = R, Z
    cross = x[:-1]*y[1:] - x[1:]*y[:-1]
    A = 0.5 * np.sum(cross)

    if np.isclose(A, 0.0):
        return 0.0, np.nan, np.nan

    Cx = (1.0 / (6.0 * A)) * np.sum((x[:-1] + x[1:]) * cross)
    Cy = (1.0 / (6.0 * A)) * np.sum((y[:-1] + y[1:]) * cross)
    return A, Cx, Cy


def _clip_vertical_halfplane(R, Z, x0, keep='le'):
    """
    Sutherland–Hodgman clip of polygon (R,Z) against vertical line R = x0.
    keep:
        'le' -> keep R <= x0 (inboard)
        'ge' -> keep R >= x0 (outboard)
    Returns closed polygon (R_clipped, Z_clipped). Empty arrays if fully clipped.
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    if (R[0] != R[-1]) or (Z[0] != Z[-1]):
        R = np.concatenate([R, R[:1]])
        Z = np.concatenate([Z, Z[:1]])

    def inside(x):
        return (x <= x0 + 1e-12) if keep == 'le' else (x >= x0 - 1e-12)

    Rout = []
    Zout = []

    for i in range(len(R) - 1):
        x1, y1 = R[i],   Z[i]
        x2, y2 = R[i+1], Z[i+1]
        i1, i2 = inside(x1), inside(x2)

        if i1 and i2:
            # both inside
            Rout.append(x2); Zout.append(y2)
        elif i1 and not i2:
            # exiting: add intersection
            if not np.isclose(x2, x1):
                t = (x0 - x1) / (x2 - x1)
                yi = y1 + t * (y2 - y1)
                Rout.append(x0); Zout.append(yi)
        elif (not i1) and i2:
            # entering: add intersection, then the endpoint
            if not np.isclose(x2, x1):
                t = (x0 - x1) / (x2 - x1)
                yi = y1 + t * (y2 - y1)
                Rout.append(x0); Zout.append(yi)
            Rout.append(x2); Zout.append(y2)
        else:
            # both outside: add nothing
            pass

    if len(Rout) == 0:
        return np.array([]), np.array([])

    # Close polygon
    if (Rout[0] != Rout[-1]) or (Zout[0] != Zout[-1]):
        Rout.append(Rout[0]); Zout.append(Zout[0])

    return np.array(Rout), np.array(Zout)
