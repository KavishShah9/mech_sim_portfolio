"""Rerun every study in this repository."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    "01_plate_hole_fea/plate_hole_fea.py",
    "02_fin_heat_transfer/fin_heat_transfer.py",
    "03_cadquery_bracket/bracket.py",
    "04_pybamm_battery/battery_thermal.py",
]

for s in SCRIPTS:
    print(f"\n=== {s} ===", flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, s)], check=True)
print("\nAll studies finished.")
