"""
Straight rectangular aluminum fin, steady 2D conduction with convection.

Finite volume model on the fin's length by thickness cross section (per metre of depth).
Base held at a fixed temperature, convection on both faces and the tip.

Checks:
  1. Heat rate and tip temperature against the classic 1D fin solution with a convective tip.
  2. Grid convergence of the heat rate.
Design study:
  Heat rate and fin efficiency against fin length, to find where extra length stops paying off.
"""
import json
import os

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

K = 205.0       # W/m K, aluminum
H = 40.0        # W/m^2 K, forced air
T_BASE = 80.0   # C
T_INF = 25.0    # C
THICK = 0.002   # m, fin thickness


def solve_2d(L, t=THICK, nx=200, ny=10):
    """Cell centred finite volume solution. Returns temperature field and base heat rate per metre depth."""
    dx, dy = L / nx, t / ny
    n = nx * ny
    idx = lambda i, j: j * nx + i
    rows, cols, vals = [], [], []
    b = np.zeros(n)
    g_x = K * dy / dx           # conductance between x neighbours (per metre depth)
    g_y = K * dx / dy           # conductance between y neighbours
    # convection on a face in series with half cell conduction
    g_face_y = 1.0 / (1.0 / (H * dx) + (dy / 2) / (K * dx))
    g_face_x = 1.0 / (1.0 / (H * dy) + (dx / 2) / (K * dy))
    g_base = K * dy / (dx / 2)

    def add(p, q, g):
        rows.extend([p, p]); cols.extend([p, q]); vals.extend([g, -g])

    for j in range(ny):
        for i in range(nx):
            p = idx(i, j)
            diag_extra = 0.0
            if i > 0:
                add(p, idx(i - 1, j), g_x)
            else:
                diag_extra += g_base; b[p] += g_base * T_BASE
            if i < nx - 1:
                add(p, idx(i + 1, j), g_x)
            else:
                diag_extra += g_face_x; b[p] += g_face_x * T_INF
            if j > 0:
                add(p, idx(i, j - 1), g_y)
            else:
                diag_extra += g_face_y; b[p] += g_face_y * T_INF
            if j < ny - 1:
                add(p, idx(i, j + 1), g_y)
            else:
                diag_extra += g_face_y; b[p] += g_face_y * T_INF
            rows.append(p); cols.append(p); vals.append(diag_extra)

    A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    T = spla.spsolve(A.tocsc(), b).reshape(ny, nx)
    q_base = np.sum(g_base * (T_BASE - T[:, 0]))
    x = (np.arange(nx) + 0.5) * dx
    return T, q_base, x


def fin_1d(L, t=THICK):
    """Incropera Table 3.4 case A, convective tip, per metre depth."""
    P = 2.0
    Ac = t
    m = np.sqrt(H * P / (K * Ac))
    M = np.sqrt(H * P * K * Ac) * (T_BASE - T_INF)
    r = H / (m * K)
    q = M * (np.sinh(m * L) + r * np.cosh(m * L)) / (np.cosh(m * L) + r * np.sinh(m * L))
    theta_tip = (T_BASE - T_INF) / (np.cosh(m * L) + r * np.sinh(m * L))
    return q, T_INF + theta_tip, m


def main():
    res = {}
    L = 0.05
    T, q2d, x = solve_2d(L)
    q1d, t_tip_1d, m = fin_1d(L)
    t_tip_2d = T[:, -1].mean()
    biot = H * (THICK / 2) / K
    res["base_case"] = {
        "length_mm": L * 1000, "thickness_mm": THICK * 1000, "k_W_mK": K, "h_W_m2K": H,
        "T_base_C": T_BASE, "T_air_C": T_INF, "Biot_half_thickness": round(biot, 6),
        "q_2d_W_per_m": round(q2d, 3), "q_1d_W_per_m": round(q1d, 3),
        "q_error_pct": round(100 * (q2d - q1d) / q1d, 3),
        "tip_T_2d_C": round(t_tip_2d, 3), "tip_T_1d_C": round(t_tip_1d, 3),
    }

    conv = []
    for nx in [10, 20, 50, 100, 200, 400]:
        _, q, _ = solve_2d(L, nx=nx, ny=max(4, nx // 20))
        conv.append({"nx": nx, "q_W_per_m": round(q, 4)})
    res["grid_convergence"] = conv

    # 1D analytic temperature profile for comparison
    xs = np.linspace(0, L, 200)
    r = H / (m * K)
    theta = (np.cosh(m * (L - xs)) + r * np.sinh(m * (L - xs))) / (np.cosh(m * L) + r * np.sinh(m * L))
    T1d = T_INF + (T_BASE - T_INF) * theta

    fig, axs = plt.subplots(2, 1, figsize=(7, 5.5), gridspec_kw={"height_ratios": [1, 1.4]})
    im = axs[0].imshow(T, aspect="auto", origin="lower", extent=[0, L * 1000, 0, THICK * 1000], cmap="inferno")
    fig.colorbar(im, ax=axs[0], label="T, C")
    axs[0].set_xlabel("x, mm")
    axs[0].set_ylabel("y, mm")
    axs[0].set_title("2D temperature field (thickness stretched for display)")
    axs[1].plot(xs * 1000, T1d, color="#c0504d", label="1D analytical fin")
    axs[1].plot(x[::10] * 1000, T.mean(axis=0)[::10], "o", ms=4, color="#1f3864", label="2D finite volume")
    axs[1].set_xlabel("Distance from base, mm")
    axs[1].set_ylabel("Temperature, C")
    axs[1].legend()
    axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fin_temperature.png"), dpi=150)
    plt.close(fig)

    # Design study: length sweep
    lengths = np.linspace(0.005, 0.15, 30)
    q_fv, eff = [], []
    for Li in lengths:
        _, qi, _ = solve_2d(Li, nx=150, ny=6)
        q_fv.append(qi)
        q_max = H * (2 * Li + THICK) * (T_BASE - T_INF)
        eff.append(qi / q_max)
    q_fv = np.array(q_fv)
    q_inf = np.sqrt(H * 2 * K * THICK) * (T_BASE - T_INF)
    L95 = float(lengths[np.argmax(q_fv >= 0.95 * q_inf)]) if np.any(q_fv >= 0.95 * q_inf) else None
    res["length_study"] = {
        "q_infinite_fin_W_per_m": round(q_inf, 2),
        "length_for_95pct_of_max_mm": round(L95 * 1000, 1) if L95 else None,
        "mL_for_base_case": round(m * L, 3),
        "efficiency_at_50mm": round(float(np.interp(0.05, lengths, eff)), 3),
    }

    fig, ax1 = plt.subplots(figsize=(6.5, 4))
    ax1.plot(lengths * 1000, q_fv, color="#1f3864", label="Heat rate")
    ax1.axhline(q_inf, color="#1f3864", ls=":", lw=1)
    ax1.set_xlabel("Fin length, mm")
    ax1.set_ylabel("Heat rate, W per m depth", color="#1f3864")
    ax2 = ax1.twinx()
    ax2.plot(lengths * 1000, eff, color="#c0504d", label="Efficiency")
    ax2.set_ylabel("Fin efficiency", color="#c0504d")
    ax2.set_ylim(0, 1.05)
    ax1.set_title("Longer fins give diminishing returns")
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fin_length_study.png"), dpi=150)
    plt.close(fig)

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
