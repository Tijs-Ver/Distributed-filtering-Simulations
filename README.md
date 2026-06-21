# Distributed-filtering-Simulations
MATLAB implementation of centralized Kalman filtering, Covariance Intersection (CI), and Split Covariance Intersection (SCI) for distributed state estimation, including Monte Carlo simulations and thesis figures.
Thesis Code Repository

This repository contains all MATLAB code developed and used during the thesis project on distributed filtering for ultra large-scale multi-agent
systems

The repository includes implementations of:

Centralized Kalman Filter (CKF)
Covariance Intersection (CI)
Split Covariance Intersection (SCI)

Applied to a problem of 1D heat conduction.
The repository contains:
Plot generation scripts
Monte Carlo simulation scripts
Monte Carlo error processing scripts
Folders containing plots used in the thesis
Folders containing data used to obtain plots

The repository is structured as:
Centralized_Kalman/
├── Main simulation script
├── Filter implementation
└── Supporting functions

CI/
├── Main simulation script
├── Filter implementation
└── Supporting functions

SCI/
├── Main simulation script
├── Filter implementation
└── Supporting functions

Monte_Carlo_Data/
├── CI/
├── SCI/
└── CKF/

Plots/
├── Thesis figures
└── Plot generation scripts
Purpose

The code in this repository was developed to evaluate the consistency and performance of centralized and decentralized state estimation algorithms. The provided scripts reproduce the simulations, Monte Carlo analyses, and figures presented in the thesis.

Author

Tijs Verhagen

Thesis

Master's Thesis, Mechanical Engineering, Eindhoven University of Technology, 2026.
