import struct
import numpy as np
from scipy.interpolate import CubicSpline

L_FLOW = 0.06  # flow direction [m]
V_TARGET = 1.5e-5  # target volume [m^3]
N_CTRL = 8
N_X = 80
N_THETA = 64


def params_to_radii(
    params: np.ndarray,
    n_x: int = N_X,
    L: float = L_FLOW,
) -> tuple[np.ndarray, np.ndarray]:
    n_pts = len(params)
    x_ctrl = np.linspace(0.0, L, n_pts)
    r_ctrl = np.abs(params)

    r_min = r_ctrl.max() * 0.02
    cs = CubicSpline(x_ctrl, r_ctrl, bc_type="not-a-knot")
    x_pts = np.linspace(0.0, L, n_x)
    r_pts = np.clip(cs(x_pts), r_min, None)
    if r_ctrl[0] < r_min:
        r_pts[0] = 0.0
    if r_ctrl[-1] < r_min:
        r_pts[-1] = 0.0
    return x_pts, r_pts


def compute_volume(x_pts: np.ndarray, r_pts: np.ndarray) -> float:
    return float(np.pi * np.trapz(r_pts**2, x_pts))


def scale_to_volume(
    params: np.ndarray,
    V_target: float = V_TARGET,
    L: float = L_FLOW,
) -> np.ndarray:
    x_pts, r_pts = params_to_radii(params, L=L)
    V = compute_volume(x_pts, r_pts)
    if V < 1e-20:
        return params
    return params * np.sqrt(V_target / V)


def _triangle(a, b, c) -> bytes:
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
    p_scaled = scale_to_volume(params, V_target, L)
    x_pts, r_pts = params_to_radii(p_scaled, n_x=n_x, L=L)

    thetas = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)

    verts = np.zeros((n_x, n_theta, 3))
    for i, (x, r) in enumerate(zip(x_pts, r_pts)):
        verts[i, :, 0] = x
        verts[i, :, 1] = r * cos_t
        verts[i, :, 2] = r * sin_t

    triangles: list[bytes] = []

    for i in range(n_x - 1):
        for j in range(n_theta):
            jn = (j + 1) % n_theta
            v00 = verts[i, j]
            v01 = verts[i, jn]
            v10 = verts[i + 1, j]
            v11 = verts[i + 1, jn]
            triangles.append(_triangle(v00, v10, v11))
            triangles.append(_triangle(v00, v11, v01))

    r_nose = r_pts[0]
    if r_nose > 1e-10:
        apex = np.array([x_pts[0], 0.0, 0.0])
        ring = verts[0]
        for j in range(n_theta):
            jn = (j + 1) % n_theta
            triangles.append(_triangle(apex, ring[jn], ring[j]))

    r_tail = r_pts[-1]
    if r_tail > 1e-10:
        apex = np.array([x_pts[-1], 0.0, 0.0])
        ring = verts[-1]
        for j in range(n_theta):
            jn = (j + 1) % n_theta
            triangles.append(_triangle(apex, ring[j], ring[jn]))

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


def default_params(n_ctrl: int = N_CTRL, V_target: float = V_TARGET, L: float = L_FLOW) -> np.ndarray:
    r_max_sq = 2.0 * V_target / (np.pi * L)
    r_max = np.sqrt(max(r_max_sq, 1e-20))
    n_total = n_ctrl + 2
    x_ctrl = np.linspace(0.0, L, n_total)
    return r_max * np.sin(np.pi * x_ctrl / L)


def sphere_natural_L(V_target: float = V_TARGET) -> float:
    R = (3.0 * V_target / (4.0 * np.pi)) ** (1.0 / 3.0)
    return 2.0 * R


def sphere_params(n_ctrl: int = N_CTRL, V_target: float = V_TARGET, L: float = None) -> np.ndarray:
    R = (3.0 * V_target / (4.0 * np.pi)) ** (1.0 / 3.0)
    if L is None:
        L = 2.0 * R
    n_total = n_ctrl + 2
    x_ctrl = np.linspace(0.0, L, n_total)
    return np.sqrt(np.maximum(0.0, R**2 - (x_ctrl - L / 2.0) ** 2))


INIT_SHAPES = {
    "football": default_params,
    "sphere": sphere_params,
}


def make_init_params(
    shape: str = "football",
    n_ctrl: int = N_CTRL,
    V_target: float = V_TARGET,
    L: float = L_FLOW,
) -> np.ndarray:
    fn = INIT_SHAPES.get(shape)
    if fn is None:
        raise ValueError(f"Unknown init shape '{shape}'. Choose: {list(INIT_SHAPES)}")
    return fn(n_ctrl=n_ctrl, V_target=V_target, L=L)


def param_bounds(n_ctrl: int = N_CTRL, V_target: float = V_TARGET, L: float = L_FLOW):
    r_max_sq = 2.0 * V_target / (np.pi * L)
    r_max = np.sqrt(max(r_max_sq, 1e-20))
    n = n_ctrl + 2
    return np.zeros(n), np.full(n, r_max * 4.0)


def shape_L(shape: str, V_target: float = V_TARGET, L: float = L_FLOW) -> float:
    if shape == "sphere":
        return sphere_natural_L(V_target)
    return L


if __name__ == "__main__":
    import sys
    shape = sys.argv[1] if len(sys.argv) > 1 else "football"
    L = shape_L(shape)
    p = make_init_params(shape, V_target=V_TARGET, L=L)
    out = f"test_shape_{shape}.stl"
    save_stl(p, out, V_target=V_TARGET, L=L)
    x, r = params_to_radii(p, L=L)
    V = compute_volume(x, r)
    r_max = r.max()
    print(f"shape={shape}  L={L*1000:.1f}mm  D={r_max*2*1000:.1f}mm  L/D={L/(r_max*2+1e-20):.2f}  V={V:.3e}")
    print(f"Saved {out}")
