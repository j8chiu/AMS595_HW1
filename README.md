# Monte Carlo Estimation of π

This repository contains MATLAB code for **AMS 595 / DCS 525 Project 1**: estimating π using Monte Carlo methods.  
The project demonstrates different approaches to random sampling, statistical uncertainty, and precision control.

---

## Contents
- `ams595_hw1.m` – main script that runs:
  - **Task 1** – fixed-N for-loop estimation with convergence plots and precision vs. cost analysis.
  - **Task 2** – while-loop estimation with a confidence-interval-based stopping rule to achieve target significant figures.
  - **Task 3** – reusable function `mcpi_precision_plot` with live plotting, precision control, and plot saving.

---

## Requirements
- MATLAB R2022a or later (earlier versions may also work).
- No external toolboxes required.

---

## User Guide

### 1. Getting Started
1. Download or clone this repository.
2. Open MATLAB and set the working directory to the folder containing the files.
3. Run the script:
   ```matlab
   ams595_hw1
