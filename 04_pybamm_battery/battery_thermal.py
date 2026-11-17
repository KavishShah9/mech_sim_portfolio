"""
Electrothermal discharge study of a 5 Ah cylindrical lithium ion cell (LG M50, Chen 2020 parameters) in PyBaMM.

Couples the electrochemistry to a lumped thermal model and discharges the cell at several C rates.

Checks:
  1. Capacity delivered at C/10 against the 5 Ah nominal capacity.
  2. Single Particle Model with electrolyte (SPMe) against the full Doyle Fuller Newman (DFN) model at 1C.
  3. Energy balance: heat generated in the cell against heat stored plus heat lost to the surroundings.
"""
import json
import os

import numpy as np
import pybamm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

NOMINAL_AH = 5.0
V_CUT = 2.5


def run(model_cls, c_rate, thermal=True):
    opts = {"thermal": "lumped"} if thermal else {}
    model = model_cls(options=opts)
    params = pybamm.ParameterValues("Chen2020")
    params.update({"Lower voltage cut-off [V]": V_CUT}, check_already_exists=False)
    exp = pybamm.Experiment([f"Discharge at {c_rate}C until {V_CUT} V"])
    sim = pybamm.Simulation(model, parameter_values=params, experiment=exp,
                            solver=pybamm.IDAKLUSolver())
    sol = sim.solve()
    out = {
        "t_s": sol["Time [s]"].entries,
        "V": sol["Voltage [V]"].entries,
        "Q_Ah": sol["Discharge capacity [A.h]"].entries,
    }
    if thermal:
        out["T_C"] = sol["X-averaged cell temperature [K]"].entries - 273.15
        out["T_amb_C"] = float(params["Ambient temperature [K]"]) - 273.15
        # energy balance terms for the lumped thermal model
        t = out["t_s"]
        q_gen = sol["Total heating [W]"].entries
        q_cool = sol["Surface total cooling [W]"].entries
        rho_cp = sol["Volume-averaged effective heat capacity [J.K-1.m-3]"].entries
        vol = float(params["Cell volume [m3]"])
        dTdt = np.gradient(out["T_C"], t)
        out["E_gen_J"] = float(np.trapezoid(q_gen, t))
        out["E_cool_J"] = float(np.trapezoid(q_cool, t))
        out["E_stored_J"] = float(np.trapezoid(rho_cp * vol * dTdt, t))
    return out, sol, params


def main():
    res = {}

    # 1. Capacity check at C/10
    low, _, _ = run(pybamm.lithium_ion.SPMe, 0.1, thermal=False)
    q_low = float(low["Q_Ah"][-1])
    res["capacity_check"] = {"C_rate": 0.1, "delivered_Ah": round(q_low, 3), "nominal_Ah": NOMINAL_AH,
                             "difference_pct": round(100 * (q_low - NOMINAL_AH) / NOMINAL_AH, 2)}

    # 2. SPMe against DFN at 1C
    spme, _, _ = run(pybamm.lithium_ion.SPMe, 1.0)
    dfn, _, _ = run(pybamm.lithium_ion.DFN, 1.0)
    q_grid = np.linspace(0, min(spme["Q_Ah"][-1], dfn["Q_Ah"][-1]), 200)
    v_s = np.interp(q_grid, spme["Q_Ah"], spme["V"])
    v_d = np.interp(q_grid, dfn["Q_Ah"], dfn["V"])
    rmse_mv = float(np.sqrt(np.mean((v_s - v_d) ** 2)) * 1000)
    res["model_comparison_1C"] = {
        "SPMe_capacity_Ah": round(float(spme["Q_Ah"][-1]), 3), "DFN_capacity_Ah": round(float(dfn["Q_Ah"][-1]), 3),
        "voltage_RMSE_mV": round(rmse_mv, 1),
        "SPMe_peak_T_C": round(float(spme["T_C"].max()), 2), "DFN_peak_T_C": round(float(dfn["T_C"].max()), 2),
    }

    # 3. C rate sweep with DFN and lumped thermal
    sweep = {}
    for c in [0.5, 1.0, 2.0, 3.0]:
        o, sol, params = run(pybamm.lithium_ion.DFN, c)
        o["energy_Wh"] = float(np.trapezoid(o["V"] * c * NOMINAL_AH, o["t_s"]) / 3600)
        sweep[c] = o

    res["c_rate_sweep"] = [{
        "C_rate": c, "current_A": c * NOMINAL_AH, "time_min": round(float(o["t_s"][-1] / 60), 1),
        "delivered_Ah": round(float(o["Q_Ah"][-1]), 3), "energy_Wh": round(o["energy_Wh"], 2),
        "temperature_rise_C": round(float(o["T_C"].max() - o["T_amb_C"]), 2),
        "heat_generated_J": round(o["E_gen_J"], 0),
        "heat_to_surroundings_J": round(o["E_cool_J"], 0),
        "heat_stored_J": round(o["E_stored_J"], 0),
        "energy_balance_residual_pct": round(100 * (o["E_gen_J"] + o["E_cool_J"] - o["E_stored_J"]) / o["E_gen_J"], 2),
    } for c, o in sweep.items()]

    colors = {0.5: "#9dc3e6", 1.0: "#2e75b6", 2.0: "#1f3864", 3.0: "#c0504d"}
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for c, o in sweep.items():
        axs[0].plot(o["Q_Ah"], o["V"], color=colors[c], label=f"{c:g}C")
        axs[1].plot(o["t_s"] / 60, o["T_C"], color=colors[c], label=f"{c:g}C")
    axs[0].set_xlabel("Discharge capacity, Ah")
    axs[0].set_ylabel("Terminal voltage, V")
    axs[0].set_title("Voltage vs capacity (DFN)")
    axs[1].set_xlabel("Time, min")
    axs[1].set_ylabel("Cell temperature, C")
    axs[1].set_title("Lumped cell temperature")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "c_rate_study.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(spme["Q_Ah"], spme["V"], color="#c0504d", label="SPMe")
    ax.plot(dfn["Q_Ah"], dfn["V"], "--", color="#1f3864", label="DFN")
    ax.set_xlabel("Discharge capacity, Ah")
    ax.set_ylabel("Terminal voltage, V")
    ax.set_title(f"1C discharge, SPMe vs DFN, RMSE {rmse_mv:.0f} mV")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "spme_vs_dfn.png"), dpi=150)
    plt.close(fig)

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
