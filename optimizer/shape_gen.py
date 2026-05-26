"""
Parametric STL generation for shape optimization.

Shape: 3-D body of revolution (axisymmetric).
  - Defined by radial profile r(x) along the flow axis (x).
  - Length along x is fixed to L_FLOW.
  - Volume = pi * integral(r^2, 0, L) is enforced to equal V_TARGET by uniform scaling.
  - Endpoints (nose/tail) are free parameters: blunt (large r) → high Cd,
    streamlined (r→0) → low Cd.
  - Non-zero endpoints get flat disk caps so the STL is a closed solid.

Parameters vector layout  (length = N_CTRL + 2):
  params[0]          — nose radius (x = 0)
  params[1..N_CTRL]  — interior radii (evenly spaced x)
  params[N_CTRL+1]   — tail radius (x = L)
"""

import struct
import numpy as np
from scipy.interpolate import CubicSpline

# ── defaults (can be overridden by the caller) ──────────────────────────────
L_FLOW   = 0.06   # fixed length along flow direction [m]
V_TARGET = 1.5e-5  # fixed volume [m³]
N_CTRL   = 8       # number of FREE interior control points
N_X      = 80      # longitudinal mesh resolution (STL rings)
N_THETA  = 64      # angular mesh resolution (STL sectors)


# ── geometry ─────────────────────────────────────────────────────────────────

def params_to_radii(
    params: np.ndarray,
    n_x: int = N_X,
    L: float = L_FLOW,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Map parameter vector → continuous radial profile r(x).

    params has length N_CTRL + 2:
      params[0]         = nose radius
      params[1:-1]      = interior radii
      params[-1]        = tail radius
    """
    n_pts = len(params)
    x_ctrl = np.linspace(0.0, L, n_pts)
    r_ctrl = np.abs(params)  # ensure non-negative

    cs = CubicSpline(x_ctrl, r_ctrl, bc_type="not-a-knot")
    x_pts = np.linspace(0.0, L, n_x)
    r_pts = np.clip(cs(x_pts), 0.0, None)
    return x_pts, r_pts


def compute_volume(x_pts: np.ndarray, r_pts: np.ndarray) -> float:
    return float(np.pi * np.trapz(r_pts**2, x_pts))


def scale_to_volume(
    params: np.ndarray,
    V_target: float = V_TARGET,
    L: float = L_FLOW,
) -> np.ndarray:
    """Return params scaled so that the resulting body has volume == V_target."""
    x_pts, r_pts = params_to_radii(params, L=L)
    V = compute_volume(x_pts, r_pts)
    if V < 1e-20:
        return params
    # Volume ∝ r², so scale r by sqrt(V_target/V)
    return params * np.sqrt(V_target / V)


# ── STL generation ────────────────────────────────────────────────────────────

def _triangle(a, b, c) -> bytes:
    """Pack one STL triangle (normal + 3 verts + attribute)."""
    ab = b - a
    ac = c - a
    n = np.cross(ab, ac)
    nlen = np.linalg.norm(n)
    n = n / nlen if nlen > 1e-30 else np.zeros(3)
    data = struct.pack("<fff", *n.astype(np.float32))
    for v in (a, b, c):
        data += struct.pack("<fff", *v.astype(np.float32))
    data += struct.pack("<H", 0)
    return data


def generate_stl_bytes(
    params: np.ndarray,
    V_target: float = V_TARGET,
    L: float = L_FLOW,
    n_x: int = N_X,
    n_theta: int = N_THETA,
) -> bytes:
    """
    Build a watertight binary STL for the body of revolution.

    Flat disk caps are added at x=0 and x=L when the corresponding endpoint
    radius is non-zero, keeping the mesh closed.
    """
    p_scaled = scale_to_volume(params, V_target, L)
    x_pts, r_pts = params_to_radii(p_scaled, n_x=n_x, L=L)

    thetas = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)

    # Vertex grid: verts[i, j] = (x, r*cos θ, r*sin θ)
    verts = np.zeros((n_x, n_theta, 3))
    for i, (x, r) in enumerate(zip(x_pts, r_pts)):
        verts[i, :, 0] = x
        verts[i, :, 1] = r * cos_t
        verts[i, :, 2] = r * sin_t

    triangles: list[bytes] = []

    # ── lateral surface (quad strip) ─────────────────────────────────────
    for i in range(n_x - 1):
        for j in range(n_theta):
            jn = (j + 1) % n_theta
            v00 = verts[i,   j]
            v01 = verts[i,   jn]
            v10 = verts[i+1, j]
            v11 = verts[i+1, jn]
            triangles.append(_triangle(v00, v10, v11))
            triangles.append(_triangle(v00, v11, v01))

    # ── nose cap (x = 0) ─────────────────────────────────────────────────
    r_nose = r_pts[0]
    if r_nose > 1e-10:
        apex = np.array([x_pts[0], 0.0, 0.0])
        ring = verts[0]
        for j in range(n_theta):
            jn = (j + 1) % n_theta
            # normals point toward -x
            triangles.append(_triangle(apex, ring[jn], ring[j]))

    # ── tail cap (x = L) ─────────────────────────────────────────────────
    r_tail = r_pts[-1]
    if r_tail > 1e-10:
        apex = np.array([x_pts[-1], 0.0, 0.0])
        ring = verts[-1]
        for j in range(n_theta):
            jn = (j + 1) % n_theta
            triangles.append(_triangle(apex, ring[j], ring[jn]))

    # ── pack binary STL ───────────────────────────────────────────────────
    header = b"\x00" * 80
    buf = bytearray(header)
    buf += struct.pack("<I", len(triangles))
    for tri in triangles:
        buf += tri
    return bytes(buf)


def save_stl(
    params: np.ndarray,
    path: str,
    V_target: float = V_TARGET,
    L: float = L_FLOW,
    **kwargs,
) -> None:
    data = generate_stl_bytes(params, V_target=V_target, L=L, **kwargs)
    with open(path, "wb") as f:
        f.write(data)


# ── parameter helpers ─────────────────────────────────────────────────────────

def default_params(n_ctrl: int = N_CTRL, V_target: float = V_TARGET, L: float = L_FLOW) -> np.ndarray:
    """
    Initial guess: prolate-spheroid-like shape.
    r(x) = r_max * sin(pi * x / L), which gives V = pi/2 * r_max^2 * L.
    We solve for r_max from V_target and then sample at the control points.
    """
    r_max_sq = 2.0 * V_target / (np.pi * L)
    r_max = np.sqrt(max(r_max_sq, 1e-20))
    n_total = n_ctrl + 2
    x_ctrl = np.linspace(0.0, L, n_total)
    r_ctrl = r_max * np.sin(np.pi * x_ctrl / L)
    return r_ctrl


def param_bounds(n_ctrl: int = N_CTRL, V_target: float = V_TARGET, L: float = L_FLOW):
    """Return (lower, upper) arrays for CMA-ES bounds."""
    r_max_sq = 2.0 * V_target / (np.pi * L)
    r_max = np.sqrt(max(r_max_sq, 1e-20))
    n = n_ctrl + 2
    return np.zeros(n), np.full(n, r_max * 4.0)


if __name__ == "__main__":
    # Quick smoke test: generate default shape and save
    p = default_params()
    save_stl(p, "test_shape.stl")
    x, r = params_to_radii(scale_to_volume(p))
    V = compute_volume(x, r)
    print(f"N_params={len(p)}  V={V:.3e} m³  (target={V_TARGET:.3e})")
    print("Saved test_shape.stl")
