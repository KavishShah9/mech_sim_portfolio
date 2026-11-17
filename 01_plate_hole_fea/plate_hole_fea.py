"""
Plate with a central hole under uniaxial tension.

2D plane stress finite element model written from scratch with NumPy and SciPy.
Quarter symmetry model, 4 node quadrilateral elements, 2x2 Gauss integration,
graded mapped mesh around the hole.

Checks:
  1. Mesh convergence of the stress concentration factor Kt.
  2. Kt against the Heywood fit used in Peterson's Stress Concentration Factors
     for a finite width plate with a central hole.
"""
import json
import os

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

E = 200e3      # MPa, steel
NU = 0.30
S0 = 100.0     # MPa, applied gross stress on the loaded edge
W = 100.0      # mm, plate width (perpendicular to load)
LX = 100.0     # mm, half length in the load direction (plate length 2W)


def build_mesh(a, n_th, n_r, grade=1.8):
    """Mapped quad mesh of the quarter plate [0, LX] x [0, W/2] minus a quarter hole of radius a."""
    ly = W / 2.0
    th_c = np.arctan2(ly, LX)
    n1 = max(2, int(round(n_th * th_c / (np.pi / 2))))
    n2 = n_th - n1
    th = np.concatenate([np.linspace(0, th_c, n1 + 1), np.linspace(th_c, np.pi / 2, n2 + 1)[1:]])
    s = np.linspace(0, 1, n_r + 1) ** grade

    inner = np.c_[a * np.cos(th), a * np.sin(th)]
    outer = np.empty_like(inner)
    for k, t in enumerate(th):
        if t <= th_c + 1e-12:
            outer[k] = [LX, LX * np.tan(t)]
        else:
            outer[k] = [ly / np.tan(t) if t < np.pi / 2 - 1e-12 else 0.0, ly]

    nodes = np.zeros(((n_th + 1) * (n_r + 1), 2))
    idx = lambda i, j: j * (n_th + 1) + i  # i along theta, j along radius
    for j, sj in enumerate(s):
        for i in range(n_th + 1):
            nodes[idx(i, j)] = (1 - sj) * inner[i] + sj * outer[i]

    elems = []
    for j in range(n_r):
        for i in range(n_th):
            elems.append([idx(i, j), idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)])
    elems = np.array(elems)

    info = dict(n1=n1, n_th=n_th, n_r=n_r, idx=idx)
    return nodes, elems, info


def D_plane_stress():
    c = E / (1 - NU ** 2)
    return c * np.array([[1, NU, 0], [NU, 1, 0], [0, 0, (1 - NU) / 2]])


GP = np.array([-1, 1]) / np.sqrt(3)


def shape_derivs(xi, eta):
    dN = 0.25 * np.array([
        [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
        [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)],
    ])
    return dN


def B_matrix(xy, xi, eta):
    dN = shape_derivs(xi, eta)
    J = dN @ xy
    detJ = np.linalg.det(J)
    dNxy = np.linalg.solve(J, dN)
    B = np.zeros((3, 8))
    B[0, 0::2] = dNxy[0]
    B[1, 1::2] = dNxy[1]
    B[2, 0::2] = dNxy[1]
    B[2, 1::2] = dNxy[0]
    return B, detJ


def solve(a, n_th, n_r):
    nodes, elems, info = build_mesh(a, n_th, n_r)
    D = D_plane_stress()
    nn = len(nodes)
    rows, cols, vals = [], [], []
    for e in elems:
        xy = nodes[e]
        ke = np.zeros((8, 8))
        for xi in GP:
            for eta in GP:
                B, detJ = B_matrix(xy, xi, eta)
                if detJ <= 0:
                    raise RuntimeError("Inverted element")
                ke += B.T @ D @ B * detJ
        dofs = np.ravel(np.c_[2 * e, 2 * e + 1])
        rows.append(np.repeat(dofs, 8))
        cols.append(np.tile(dofs, 8))
        vals.append(ke.ravel())
    K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * nn, 2 * nn))

    # Loads: uniform traction S0 on the right edge x = LX (outer row, theta index 0..n1)
    F = np.zeros(2 * nn)
    idx = info["idx"]
    jr = info["n_r"]
    for i in range(info["n1"]):
        n_a, n_b = idx(i, jr), idx(i + 1, jr)
        length = np.linalg.norm(nodes[n_b] - nodes[n_a])
        F[2 * n_a] += S0 * length / 2
        F[2 * n_b] += S0 * length / 2

    # Symmetry: ux = 0 on x = 0, uy = 0 on y = 0
    fixed = []
    for k, (x, y) in enumerate(nodes):
        if abs(x) < 1e-9:
            fixed.append(2 * k)
        if abs(y) < 1e-9:
            fixed.append(2 * k + 1)
    fixed = np.unique(fixed)
    free = np.setdiff1d(np.arange(2 * nn), fixed)
    U = np.zeros(2 * nn)
    U[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])

    # Nodal stresses: extrapolate Gauss point stresses to nodes, then average
    corner = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * np.sqrt(3)
    sig_sum = np.zeros((nn, 3))
    cnt = np.zeros(nn)
    for e in elems:
        xy = nodes[e]
        ue = np.ravel(np.c_[U[2 * e], U[2 * e + 1]])
        sg = []
        for (xi, eta) in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
            B, _ = B_matrix(xy, xi / np.sqrt(3), eta / np.sqrt(3))
            sg.append(D @ B @ ue)
        sg = np.array(sg)
        for k, (r, q) in enumerate(corner):
            N = 0.25 * np.array([(1 - r) * (1 - q), (1 + r) * (1 - q), (1 + r) * (1 + q), (1 - r) * (1 + q)])
            sig_sum[e[k]] += N @ sg
            cnt[e[k]] += 1
    sig = sig_sum / cnt[:, None]

    top_of_hole = idx(info["n_th"], 0)  # node at (0, a)
    kt = sig[top_of_hole, 0] / S0
    return dict(nodes=nodes, elems=elems, U=U, sig=sig, kt=kt, n_elem=len(elems), n_dof=2 * nn)


def kt_peterson(d_over_w):
    r = 1 - d_over_w
    ktn = 2 + 0.284 * r - 0.600 * r ** 2 + 1.32 * r ** 3
    return ktn / r  # gross section Kt


def main():
    results = {}

    # 1. Mesh convergence at d/W = 0.2
    a = 0.2 * W / 2
    conv = []
    for n in [8, 12, 16, 24, 32, 48, 64]:
        r = solve(a, n, n)
        conv.append((r["n_elem"], r["kt"]))
    ref = kt_peterson(0.2)
    results["convergence"] = [{"elements": n, "kt": round(k, 4)} for n, k in conv]
    results["kt_reference_dW_0.2"] = round(ref, 4)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([c[0] for c in conv], [c[1] for c in conv], "o-", color="#1f3864", label="FEA, top of hole")
    ax.axhline(ref, color="#c0504d", ls="--", label=f"Peterson (Heywood) = {ref:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("Number of elements (quarter model)")
    ax.set_ylabel("Kt (gross section)")
    ax.set_title("Mesh convergence, d/W = 0.2")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "convergence.png"), dpi=150)
    plt.close(fig)

    # 2. Sweep d/W with a converged mesh
    sweep = []
    for dw in [0.1, 0.2, 0.3, 0.4, 0.5]:
        r = solve(dw * W / 2, 48, 48)
        ref = kt_peterson(dw)
        sweep.append({"d_over_W": dw, "kt_fea": round(r["kt"], 4), "kt_peterson": round(ref, 4),
                      "error_pct": round(100 * (r["kt"] - ref) / ref, 2)})
    results["sweep"] = sweep

    fig, ax = plt.subplots(figsize=(6, 4))
    dws = np.linspace(0.05, 0.6, 100)
    ax.plot(dws, kt_peterson(dws), color="#c0504d", label="Peterson (Heywood fit)")
    ax.plot([s["d_over_W"] for s in sweep], [s["kt_fea"] for s in sweep], "o", color="#1f3864", label="FEA")
    ax.set_xlabel("Hole diameter / plate width, d/W")
    ax.set_ylabel("Kt (gross section)")
    ax.set_title("Stress concentration vs hole size")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "kt_vs_dW.png"), dpi=150)
    plt.close(fig)

    # 3. Stress field plot at d/W = 0.2
    r = solve(0.2 * W / 2, 48, 48)
    nodes, elems, sig = r["nodes"], r["elems"], r["sig"]
    tris = np.vstack([elems[:, [0, 1, 2]], elems[:, [0, 2, 3]]])
    tri = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    cf = ax.tricontourf(tri, sig[:, 0], levels=30, cmap="viridis")
    fig.colorbar(cf, ax=ax, label="Normal stress in load direction, MPa")
    ax.set_aspect("equal")
    ax.set_xlabel("x, mm")
    ax.set_ylabel("y, mm")
    ax.set_title(f"Quarter model, d/W = 0.2, applied {S0:.0f} MPa, peak {sig[:, 0].max():.0f} MPa")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_field.png"), dpi=150)
    plt.close(fig)
    results["peak_stress_MPa_dW_0.2"] = round(float(sig[:, 0].max()), 1)
    results["model"] = {"E_MPa": E, "nu": NU, "applied_stress_MPa": S0, "plate_width_mm": W,
                        "plate_length_mm": 2 * LX, "element": "Q4 plane stress, 2x2 Gauss"}

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
