import numpy as np, math

def miller_points(R0_m: float, a_m: float, delta: float, kappa: float, n: int = 400) -> np.ndarray:
    """Return Nx2 array of (R,Z) boundary points in cm for a Miller D-shape (last point repeats first)."""
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    R = R0_m + a_m*np.cos(theta + delta*np.sin(theta))
    Z = kappa*a_m*np.sin(theta)
    pts_cm = 100.0*np.column_stack([R, Z])
    return np.vstack([pts_cm, pts_cm[0]])

def _open(P):
    """Drop duplicate last vertex if closed."""
    P = np.asarray(P)
    return P[:-1] if np.allclose(P[0], P[-1]) else P

def _area_centroid(P):
    """Signed area and centroid of polygon (cm^2, cm)."""
    P = _open(P); x, y = P[:,0], P[:,1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    c = x*y1 - x1*y
    A = 0.5*np.sum(c)
    Cx = np.sum((x+x1)*c)/(6*A)
    Cy = np.sum((y+y1)*c)/(6*A)
    return A, Cx, Cy

def _ensure_ccw(P):
    """CCW orientation."""
    A, _, _ = _area_centroid(P)
    P = _open(P)
    return P if A > 0 else P[::-1]

def _normals(P):
    """Unit outward normals at vertices for CCW polygon."""
    P = _open(P)
    f = np.roll(P, -1, 0) - P
    b = P - np.roll(P, 1, 0)
    f /= (np.linalg.norm(f, axis=1, keepdims=True)+1e-16)
    b /= (np.linalg.norm(b, axis=1, keepdims=True)+1e-16)
    t = f + b
    t /= (np.linalg.norm(t, axis=1, keepdims=True)+1e-16)
    n = np.column_stack([t[:,1], -t[:,0]])
    n /= (np.linalg.norm(n, axis=1, keepdims=True)+1e-16)
    return n

def volumes(R0_m, a_m, delta, kappa, t_m, n=1000):
    """Return (V_plasma, V_outer, V_blanket) in m^3 for blanket thickness t_m (m)."""
    inner = _ensure_ccw(miller_points(R0_m, a_m, delta, kappa, n))
    A_in, Cx_in, _ = _area_centroid(inner)
    outer = inner + (100.0*t_m)*_normals(inner)
    A_out, Cx_out, _ = _area_centroid(outer)
    V_in_cm3  = 2*math.pi*Cx_in *A_in
    V_out_cm3 = 2*math.pi*Cx_out*A_out
    to_m3 = 1e-6
    return V_in_cm3*to_m3, V_out_cm3*to_m3, (V_out_cm3-V_in_cm3)*to_m3

# Example:
R0, a, kappa, delta, t = 3.30, 1.13, 1.84, 0.45, 1.02                       # meters, -, -, meters
V_in, V_out, V_bl = volumes(R0, a, delta, kappa, t, n=1200)
print("V_plasma [m^3] :", V_in)
print("V_outer  [m^3] :", V_out)
print("V_blanket[m^3] :", V_bl)
