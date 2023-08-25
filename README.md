# Integrated Population Model for WA coastal steelhead

## Summary
This repository contains code to fit a hierarchical Integrated Population Model (IPM) to multiple populations of wild winter steelhead (_Oncorhynchus mykiss_) on the Washington coast. The IPM is a statistical population dynamics models that integrates information on spawner abundances, total harvest, and spawner age structure into a combined run-reconstruction and spawner-recruitment model. It accounts for iteroparity by distinguishing maiden from repeat spawners and estimates time-varying kelt survival rates. The model also estimates population parameters such as productivity and capacity as well as time-varying recruitment residuals.

## Repository Structure 
The repository contains three sub-directories as described below:

| Directory   | Description                                           |
| ----------- | ----------------------------------------------------- |
| `code`      | Contains code to run the analyses and produce figures |
| `data`      | Contains fish and covariate data used in the analysis |
| `functions` | Contains functions for data preparation and sourcing  |

The `code` directory contains scripts to be run in the following order:

| File                 | Description                                  |
| -------------------- | -------------------------------------------- |
| `IPM_sthd_multipop`  | Code to load the fish data and fit the IPM*  |
| `posthoc_analysis`   | Post-hoc analysis for covariate selection**  |
| `manuscript_figures` | Code to produce the main results figures     |

* IPM to be fit twice: using covar_effects=TRUE and covar_effects=FALSE
** Uses IPM output without covariates to perform a covariate selection

## Dependencies
The model fitting relies on the R package 'salmonIPM'. This package can be made available upon request to the main developer Eric Buhle via GitHub (https://github.com/ebuhle) or email (<ebuhle@gmail.com>).

## Output
The code described above produces the main results figures shown in Ohlberger et al. (unpublished).
