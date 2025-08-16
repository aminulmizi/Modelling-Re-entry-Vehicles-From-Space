# Modelling Re-entry Vehicles from Space

**Repository:** `Modelling-Re-entry-Vehicles-From-Space`  
**Author / Student ID:** Aminul Mizi  
**Supervisor / School:** Wolfson School of Mechanical, Electrical and Manufacturing Engineering - Loughborough University

---

## Overview

This repository contains the code and supporting material for a final-year project investigating the aerodynamics and trajectories of atmospheric re-entry vehicles. The work combines CFD analyses (STAR-CCM+) of a Falcon-9-inspired booster (including grid-fin redesign studies) with custom Python trajectory simulators (3-DOF and 6-DOF variants). The full project report (`Modelling-reentry-vehicles-from-space-report.pdf`) describing methods, results and discussion is included.

---

## Contents

- `Modelling-reentry-vehicles-from-space-report.pdf`: Full project report describing methodology, CFD setup and results, grid-fin redesign, and trajectory simulation.  
- `3DOF-TRAJECTORY.py`: 3-degree-of-freedom trajectory solver (Runge-Kutta integration, atmospheric model, aerodynamic data interpolation).  
- `6D0F-TRAJECTORY.py`: Intended 6-degree-of-freedom dynamics solver (check header for usage and dependencies).  
- `IN-VACUO.py`: Helper/test script for vacuum or high-altitude initial conditions.  
- `YplusCalc.m`: MATLAB script to compute first-cell height and prism layer settings for a target y+ (mesh guidance).

---

## Quick summary of the work

- **CFD validation:** sphere validation cases were used to select turbulence models; k–ω SST was chosen for best agreement with reference drag data.  
- **Vehicle analysis:** a scaled Falcon-9 inspired geometry was simulated across transonic to hypersonic points. Drag coefficients and flow fields were extracted for use in trajectory code.  
- **Grid-fin redesign:** larger contoured cells and vortex generators reduced the transonic drag rise, improving flow through the grid cells; drag reductions of ≈1.6–5.2% were observed in studied cases.  
- **Trajectory modelling:** a Python 3-DOF solver uses interpolated aerodynamic data and an atmospheric model to project flight paths and compare design variants.

(Full details, tables and plots are in `Modelling-reentry-vehicles-from-space-report.pdf`.)

---

## Dependencies

> Check the top of each script for exact version requirements.

**Python (recommended)**  
- Python 3.8+ (3.10 recommended)  
- `numpy`, `scipy`, `pandas`, `matplotlib`  
- `tqdm` (optional)

Install with pip:
```bash
pip install numpy scipy pandas matplotlib tqdm
```

## Y+ and mesh guidance

The included `YplusCalc.m` MATLAB script computes the first-cell height and prism layer parameters for a target **y+**. It was used to guide prism layer settings for k–ω SST simulations in STAR-CCM+.  
Run the script in MATLAB and supply the required inputs: **density**, **viscosity**, **reference length**, **velocity**, and **target y+**.

---

## Reproducing / extending the CFD work

The STAR-CCM+ project files and meshes are **not** part of this repository. The report documents:

- geometry & scaling approach  
- mesh strategy (polyhedral + prism layers + AMR controls)  
- solver settings (numerical schemes, turbulence model selection, boundary conditions)  
- altitude-dependent gas properties and case points

To reproduce the CFD results you will need access to STAR-CCM+ and the geometry/mesh files (or you will need to recreate them). Consult the report for detailed solver and mesh setup information.

---

## Results & key takeaways

- **k–ω SST** matched sphere validation best and was used for re-entry analyses.  
- The Falcon-9 inspired vehicle exhibits a **transonic drag rise**, with grid-fins a significant contributor.  
- Redesigned grid-fins (contoured cells + vortex generators) reduced transonic drag and improved through-cell flow.  
- Translating CFD improvements into 3-DOF trajectories yielded **modest gains** — ballistic coefficient and integrated vehicle aerodynamics remain dominant factors.

---

## Notes, limitations & next steps

- STAR-CCM+ project and NX CAD files are **not included**; these are required for full reproduction.  
- The 6-DOF solver may require additional input formatting or dependencies — check its header comments.  
- Suggested future work: real-time CFD coupling, extended dynamic modelling, wind-tunnel or flight testing, and structural/manufacturability evaluation of redesigned fins.

---

## Contributing

Want to help? Open a pull request or issue. Useful contributions include clearer run instructions, a sample `aero_table.csv`, unit tests, or exported STAR-CCM+ case metadata.

---

## License

This repository currently has **no license file**. If you want a permissive license, consider adding an `LICENSE` (e.g. MIT).

