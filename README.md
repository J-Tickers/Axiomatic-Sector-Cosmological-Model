
# Axiomatic Sector Cosmological Model (ASCM)
A Boundary-Coupled Dual-Manifold Effective Field Theory

This repository contains the complete analytical proofs, multi-epoch Radau ODE integrations, and dynamic Fisher-Laplace forecasting architecture for the ASCM framework. The codebase mathematically verifies the 88-to-10 mode Dirac-Bergmann reduction, the scale-dependent Kadanoff block coarse-graining, and the Real-Intermediate-State (RIS) subtracted baryogenesis transport network.

## Execution and Verification
To reproduce the mathematical integrity audits and Stage-IV empirical forecasts, clone this repository and execute the diagnostic suite.
* Run `python3 unified_constraint_surface_audit.py` to verify the dimensionless constraint norm remains strictly $\mathcal{E}_C < 10^{-10}$ across the $N \in [-60.0, 0.0]$ integration span.
* Run `python3 bulk_spectral_bvp_solver.py` to isolate the massless graviton zero mode ($m_0^2 = 0.0000$) and confirm the positive Kaluza-Klein mass gap.
* Run `python3 phase4_fisher_laplace_evidence.py` to evaluate the dynamic Jacobian and generate the Stage-IV $+20.54$ Laplace log-evidence.

## Heavy Data Artifacts
All symbolic computer-algebra proofs, stiff Radau boundary-value and Boltzmann transport integrators, and dynamic Fisher-matrix and Laplace-evidence likelihood pipelines are permanently archived and publicly available at GitHub:
* **GitHub Repository:** https://github.com/J-Tickers/Axiomatic-Sector-Cosmological-Model

The 200x300 PM grid kernels, 2LPT transfer tables, and Gadget-4 TreePM numerical artifacts required for hydrodynamic scaling are hosted externally via Figshare to preserve repository performance:
* **Numerical Artifacts DOI:** 10.6084/m9.figshare.34734825
* **Numerical Artifacts Live URL:** https://doi.org/10.6084/m9.figshare.34734825
* **Journal Manuscript and Companion Manual Artifacts DOI:** 10.6084/m9.figshare.34746873
* **Official ORCID Identifier:** 0009-0005-7943-2176