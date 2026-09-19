# Need → Access → Timely Response → Health Outcome

Reproducible R code for the v6 policy-process and target trial emulation analyses using HRS and CLHLS derived person-trial inputs.

The repository contains analysis code only. Raw data, participant-level derived data, fitted RDS objects, local paths and logs are excluded from version control.

## Analysis structure

1. HRS 2 × 2 Need/outcome construct decomposition.
2. Need-component and stable-health sensitivity analyses.
3. Need-severity and Need-standardized Access equity analyses.
4. Access-stratified target trial emulation on the additive scale.
5. HRS/CLHLS common-construct transportability analysis.
6. Post-response health-system outcomes.
7. Observational policy-coverage scenarios.
8. Nature-style figures with auditable source-data exports.

CLHLS is used only as an external replication cohort. Cross-cohort differences are descriptive heterogeneity and are not interpreted as causal health-system effects.

## Required inputs

Two previously derived RDS inputs are required:

- HRS all person-trial candidates.
- CLHLS all person-trial candidates.

Supply paths through environment variables; do not place restricted data in the repository.

```powershell
$env:NATR_REPO_ROOT='path/to/Need-Access-TimelyResponse-TTE-Frailty'
$env:NATR_OUTPUT_ROOT='path/to/private/output'
$env:NATR_HRS_CANDIDATES='path/to/02_all_person_trial_candidates.rds'
$env:NATR_CLHLS_CANDIDATES='path/to/02_clhls_all_trial_candidates.rds'
Rscript run_all.R
```

## Software

R 4.4 or later is recommended. Required packages are checked in `R/config.R`: `dplyr`, `tidyr`, `purrr`, `readr`, `survey`, `sandwich`, `ggplot2`, `patchwork`, `scales`, `svglite` and `ragg`.

## Outputs

The private output directory contains analysis tables, model objects, validation records, figure source data and figures in SVG, PDF, 600-dpi TIFF and PNG. Figure contracts are documented in `docs/figure_contracts.md`.

The policy-coverage module is an observational policy-impact scenario analysis, not a policy implementation trial.
