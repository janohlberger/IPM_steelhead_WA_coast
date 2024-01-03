# Integrated Population Model for WA coastal steelhead

## Summary
This repository contains code to fit a hierarchical Integrated Population Model (IPM) to multiple populations of wild winter steelhead (_Oncorhynchus mykiss_) on the Washington coast. The IPM is a statistical population dynamics models that integrates information on spawner abundances, total harvest, and spawner age structure into a combined run-reconstruction and spawner-recruitment model. It accounts for iteroparity by distinguishing maiden from repeat spawners and estimates time-varying kelt survival rates. The model also estimates time-varying recruitment residuals and population parameters such as productivity and capacity.

## Repository Structure 
The repository contains the following sub-directories:

| Directory   | Description                                           |
| ----------- | ----------------------------------------------------- |
| `code`      | Contains code to run the analyses and produce figures |
| `data`      | Contains fish and covariate data used in the analysis |

An `output` folder is created when running the first analysis script.

## Analyses
The `code` directory contains the following scripts:

| File                    | Description                               |
| ----------------------- | ----------------------------------------- |
| `IPM_sthd_multipop.R`   | Fits integrated population models (IPMs)  |
| `posthoc_analysis.R`    | Post-hoc analysis for covariate selection |
| `manuscript_figures.R`  | Produces the main text figures            |
| `supplementary_tables.R`| Produces the supplementary tables         |
| `year_effects.R`        | Compares models with/without year effects |

To produce the full set of results presented in the manuscript, four different versions of the IPM need to be fit - using combinations of covar_effects=TRUE or covar_effects=FALSE and year_effect=TRUE or year_effect=FALSE. The main covariate model presented in the manuscript figures uses covar_effects=TRUE and year_effect=FALSE.

To reproduce the full analysis, run scripts in the following order::

1. `IPM_sthd_multipop.R` using covar_effects=FALSE and year_effect=FALSE
2. `posthoc_analysis.R` to identify covariates to be included in the IPM
3. `IPM_sthd_multipop.R` using covar_effects=TRUE and year_effect=FALSE
4. `manuscript_figures.R` to produce the main text figures
5. `supplementary_tables.R` to produce the supplementary tables
6. `IPM_sthd_multipop.R` using year_effect=TRUE with/without covariates
7. `year_effects.R` to compare year effects with vs. without covariates

## Dependencies
The model fitting relies on the R package 'salmonIPM'. This package can be made available upon request to the main developer Eric Buhle via GitHub (https://github.com/ebuhle) or email (<eric.buhle@noaa.gov>).

## Output
The code described above produces the main results figures shown in Ohlberger et al. (unpublished manuscript).
