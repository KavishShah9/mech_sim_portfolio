"""
Parametric L bracket in CadQuery, sized with beam theory.

The bracket bolts to a wall through two M6 holes in its flange and carries a tip load on its arm.
For each candidate thickness the script builds the solid, reads mass properties from the
OpenCascade kernel, checks root bending stress and tip deflection, and keeps the lightest
design that meets the safety factor. The chosen design is exported as STEP, STL and SVG.

Checks:
  1. CadQuery volume against a closed form volume of the same profile, holes and fillet.
  2. Stress and deflection from Euler Bernoulli beam theory with the arm root as the fixed end.
"""
import json
import os

import numpy as np
import cadquery as cq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

# Geometry (mm)
WIDTH = 40.0        # bracket width
ARM = 90.0          # overall arm length from the wall face
FLANGE = 60.0       # flange height
FILLET = 6.0        # inside corner fillet radius
BOLT_D = 6.6        # M6 clearance
TIP_D = 8.4         # M8 clearance for the load pin
TIP_X = ARM - 10.0  # load pin position from the wall face

# Material and load: 6061 T6 aluminum
E = 68.9e3          # MPa
SY = 276.0          # MPa
RHO = 2.70e-3       # g per mm^3
F = 500.0           # N at the load pin
FOS = 2.0


def build(t):
    pts = [(0, 0), (ARM, 0), (ARM, t), (t, t), (t, FLANGE), (0, FLANGE)]
    body = cq.Workplane("XZ").polyline(pts).close().extrude(WIDTH / 2, both=True)
    body = body.edges("|Y").edges(cq.selectors.NearestToPointSelector((t + 0.1, 0, t + 0.1))).fillet(FILLET)
    bolts = (cq.Workplane("YZ").pushPoints([(-WIDTH / 4, FLANGE * 0.62), (WIDTH / 4, FLANGE * 0.62)])
             .circle(BOLT_D / 2).extrude(t))
    pin = cq.Workplane("XY").pushPoints([(TIP_X, 0)]).circle(TIP_D / 2).extrude(t)
    return body.cut(bolts).cut(pin)


def volume_closed_form(t):
    area = ARM * t + (FLANGE - t) * t + (FILLET ** 2 - np.pi * FILLET ** 2 / 4)
    holes = 2 * np.pi * (BOLT_D / 2) ** 2 * t + np.pi * (TIP_D / 2) ** 2 * t
    return area * WIDTH - holes


def beam_check(t):
    lever = TIP_X - t - FILLET           # from end of fillet to the load pin
    M = F * lever                        # N mm
    I = WIDTH * t ** 3 / 12
    sigma = M * (t / 2) / I              # MPa
    defl = F * lever ** 3 / (3 * E * I)  # mm
    return sigma, defl


def main():
    rows = []
    for t in np.arange(4.0, 10.5, 0.5):
        part = build(t)
        solid = part.val()
        v_cq = solid.Volume()
        v_cf = volume_closed_form(t)
        sigma, defl = beam_check(t)
        rows.append({
            "t_mm": float(t), "volume_mm3": round(v_cq, 1), "volume_closed_form_mm3": round(v_cf, 1),
            "volume_error_pct": round(100 * (v_cq - v_cf) / v_cf, 4),
            "mass_g": round(v_cq * RHO, 1), "root_stress_MPa": round(sigma, 1),
            "tip_deflection_mm": round(defl, 3), "safety_factor": round(SY / sigma, 2),
            "passes": bool(SY / sigma >= FOS),
        })

    passing = [r for r in rows if r["passes"]]
    best = min(passing, key=lambda r: r["mass_g"])
    t = best["t_mm"]
    part = build(t)
    cq.exporters.export(part, os.path.join(OUT, "bracket.step"))
    cq.exporters.export(part, os.path.join(OUT, "bracket.stl"))
    cq.exporters.export(part, os.path.join(OUT, "bracket_iso.svg"), opt={
        "width": 520, "height": 360, "marginLeft": 40, "marginTop": 30,
        "projectionDir": (1.0, -1.4, 0.9), "showAxes": False, "strokeWidth": 0.4,
        "strokeColor": (31, 56, 100), "hiddenColor": (180, 180, 180), "showHidden": True,
    })
    bb = part.val().BoundingBox()
    com = cq.Shape.centerOfMass(part.val())

    ts = [r["t_mm"] for r in rows]
    fig, ax1 = plt.subplots(figsize=(6.5, 4))
    ax1.plot(ts, [r["root_stress_MPa"] for r in rows], "o-", color="#1f3864", label="Root bending stress")
    ax1.axhline(SY / FOS, color="#1f3864", ls="--", lw=1, label=f"Allowable, yield / {FOS:.0f}")
    ax1.set_xlabel("Thickness, mm")
    ax1.set_ylabel("Stress, MPa", color="#1f3864")
    ax2 = ax1.twinx()
    ax2.plot(ts, [r["mass_g"] for r in rows], "s-", color="#c0504d", label="Mass")
    ax2.set_ylabel("Mass, g", color="#c0504d")
    ax1.axvline(t, color="gray", ls=":", lw=1)
    ax1.set_title(f"Lightest passing design: t = {t:.1f} mm, {best['mass_g']:.0f} g")
    ax1.grid(alpha=0.3)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper center")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "thickness_study.png"), dpi=150)
    plt.close(fig)

    res = {
        "inputs": {"width_mm": WIDTH, "arm_mm": ARM, "flange_mm": FLANGE, "fillet_mm": FILLET,
                   "load_N": F, "material": "6061 T6", "yield_MPa": SY, "required_safety_factor": FOS},
        "chosen": best,
        "bounding_box_mm": [round(bb.xlen, 2), round(bb.ylen, 2), round(bb.zlen, 2)],
        "center_of_mass_mm": [round(com.x, 2), round(com.y, 2), round(com.z, 2)],
        "max_volume_error_pct": max(abs(r["volume_error_pct"]) for r in rows),
        "sweep": rows,
    }
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps({k: v for k, v in res.items() if k != "sweep"}, indent=2))


if __name__ == "__main__":
    main()
