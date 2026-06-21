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

Centralized Kalman filter/
-  KF_Plotting (Single run creating all plots saved to plots)
-  KF_Monte_Carlo (Monte Carlo simulation saving the error to /data)
-  KF_Error_Processing (Extracts information from MC and saves to /data)
-  plots folders/ (contains plots used)
-  data folder/ (contains data used)
Distributed Kalman Filter CI/
-  CI_Plotting (Single run creating all plots saved to plots)
-  CI_Monte_Carlo (Monte Carlo simulation saving the error to /data)
-  CI_Error_Processing (Extracts information from MC and saves to /data)
-  plots folders/ (contains plots used)
-  data folder/ (contains data used)
-  OverlappingCI (function script found at: https://github.com/decenter2021/OCI)
SCI/
-  SCI_Plotting (Single run creating all plots saved to plots)
-  SCI_Monte_Carlo (Monte Carlo simulation saving the error to /data)
-  SCI_Error_Processing (Extracts information from MC and saves to /data)
-  plots folders/ (contains plots used)
-  data folder/ (contains data used)
-  OverlappingCI (function script found at: https://github.com/decenter2021/OCI)

Purpose

The code in this repository was developed to evaluate the performance of different estimation algorithms. The provided scripts reproduce the simulations, Monte Carlo analyses, and figures presented in the thesis.

Author

Tijs Verhagen

Thesis

Bachelor's Thesis, Mechanical Engineering, Eindhoven University of Technology, 2026.
