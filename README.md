# DRAGON-DOCK: AA279D Control of Distributed Space Systems

## Overview
DRAGON-DOCK is a project for Stanford's AA279D course on Spacecraft Formation-Flying and Rendezvous. It develops a framework for spacecraft rendezvous and docking in low Earth orbit, integrating navigation (using an Unscented Kalman Filter) and continuous control laws for a chief-deputy spacecraft pair. The project includes MATLAB simulations and a comprehensive final report.

## Repository Structure
- **control_helpers/**: Helper functions for control algorithms.
- **propagators/**: Orbital propagation functions (e.g., Keplerian, J2-perturbed).
- **utils/**: Utility functions for conversions (e.g., ECI to RTN, orbital elements).
- **Figures/**: Plots and visualizations generated from simulations.
- **DRAGON-DOCK - AA 279D Final Report.pdf**: Final project report detailing methodology, results, and analysis.
- **ps[1-9]_main.m**: Main scripts for Problem Sets 1-9, implementing simulations and analyses.
- **control_results.mat**: Saved MATLAB results for control simulations.

## Setup and Usage
1. **Requirements**: MATLAB (R2020a or later recommended).
2. **Clone the Repository**:
   ```bash
   git clone https://github.com/b9check/DRAGON-DOCK-Stanford-AA279D-Control-of-Distributed-Space-Systems.git
   ```
3. **Run Scripts**: Open MATLAB, navigate to the repository directory, and execute `ps[1-9]_main.m` scripts for each problem set. Ensure `control_helpers`, `propagators`, and `utils` folders are in the MATLAB path.
4. **View Results**: Outputs include plots (saved in `Figures/`) and numerical results (e.g., `control_results.mat`).

## Final Report
The final report (`DRAGON-DOCK - AA 279D Final Report.pdf`) provides a detailed explanation of the project, including mission objectives, dynamics modeling, navigation, control strategies, and results. Refer to it for in-depth methodology and conclusions.

## How to Cite
If you use this code in your work, please cite it as follows:

Carlhammar, A., & Check, B. (2025). DRAGON-DOCK: Code for Spacecraft Rendezvous and Docking in AA279D Distributed Space Systems (Version 1.0) [Computer software]. Stanford University. https://github.com/b9check/DRAGON-DOCK-Stanford-AA279D-Control-of-Distributed-Space-Systems

For BibTeX users:
```bibtex
@software{Carlhammar_Check_2025_DRAGON_DOCK,
  author = {Alexandre Carlhammar and Brian Check},
  title = {{DRAGON-DOCK: Code for Spacecraft Rendezvous and Docking in AA279D Distributed Space Systems}},
  year = {2025},
  publisher = {Stanford University},
  url = {https://github.com/b9check/DRAGON-DOCK-Stanford-AA279D-Control-of-Distributed-Space-Systems},
  note = {Developed for AA279D: Distributed Space Systems, taught by Prof. Simone D'Amico},
  version = {1.0}
}
```

## Contributors
- Alexandre Carlhammar ([@alex-crlhmmr](https://github.com/alex-crlhmmr), acarlham@stanford.edu)
- Brian Check ([@b9check](https://github.com/b9check), bcheck@stanford.edu)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details (if added).
