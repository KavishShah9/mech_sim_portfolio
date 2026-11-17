# Mechanical Engineering Simulation Portfolio

Kavish Shah. BEng Mechanical Engineering, MEng Industrial Engineering.

I built four engineering studies in Python with coding agents. Each script sets up the physics, runs the model, checks the answer against theory or a second model, and writes its plots and a `results.json` file.

## Projects

1. **[Plate with a hole under tension](01_plate_hole_fea)**: a 2D plane stress finite element solver written from scratch with NumPy and SciPy. Kt matches the Peterson chart within 1.2% for hole sizes from 10% to 50% of plate width.
2. **[Aluminum cooling fin](02_fin_heat_transfer)**: a 2D finite volume conduction and convection model. The heat rate matches the 1D fin solution within 0.01%, and a length study shows where extra fin length stops paying off.
3. **[Parametric L bracket in CadQuery](03_cadquery_bracket)**: builds the bracket as a solid at every thickness from 4 to 10 mm, sizes it with beam theory and exports the lightest passing design as STEP, STL and SVG. The kernel volume matches the closed form volume within 0.0001%.
4. **[Battery cell electrothermal discharge](04_pybamm_battery)**: PyBaMM DFN and SPMe models of a 5 Ah LG M50 cell with lumped thermal coupling at four C rates. The energy balance closes within 0.03%.

## Run it

```bash
pip install -r requirements.txt
python run_all.py
```

Each study also runs on its own, for example `python 01_plate_hole_fea/plate_hole_fea.py`.

## How I check a model

1. I write the physics first: geometry, materials, loads, boundary conditions and the answer I expect from a hand calculation.
2. A coding agent drafts the script from that spec. I read every line, with extra attention on element ordering, units and sign conventions.
3. I compare the result with theory, a second model or a conservation law, and I run a convergence study before I use any number.
4. Inputs, results and plots live next to the code so anyone can rerun the study.
