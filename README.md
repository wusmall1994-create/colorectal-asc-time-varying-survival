# Time-varying mortality and conditional survival in colorectal adenosquamous carcinoma

R code for the population-based SEER analysis of colorectal adenosquamous carcinoma (ASC) versus adenocarcinoma not otherwise specified (NOS). The repository contains statistical code only. It contains no SEER records, patient identifiers, intermediate patient-level objects, manuscript files, or result files.

## Requirements

- R 4.6.1 or a compatible recent R release
- Packages: `data.table`, `survival`, `splines`, `parallel`, `ggplot2`
- Authorised access to the SEER Research Data used in the study

## Input files

Export the required SEER variables as described in `R/01_prepare_data.R`. Set these environment variables to the local export paths:

```text
SEER_ASC_FILE       colorectal ASC, ICD-O-3 morphology 8560
SEER_AC_NOS_FILE    colorectal adenocarcinoma NOS, morphology 8140
SEER_ALL_HIST_FILE  all eligible colorectal histologies for the expanded comparator
```

`SEER_ALL_HIST_FILE` is needed only for scripts 07–09. The scripts never transmit SEER data and do not write patient identifiers to the public repository.

## Run order

Run from the repository root:

```r
source("R/01_prepare_data.R")
source("R/02_bootstrap_conditional_survival.R")
source("R/03_sensitivity_diagnostics.R")
source("R/04_time_interaction.R")
source("R/05_figure_conditional_survival.R")
source("R/06_figure_time_specific_hr.R")
source("R/07_expanded_comparator_audit.R")
source("R/08_compare_8140_exports.R")
source("R/09_expanded_comparator_analysis.R")
source("R/10_group_bootstrap_intervals.R")
source("R/11_transparency_analyses.R")
```

The bootstrap uses 2,000 samples drawn independently within histology. The main conditional-survival model permits separate histology-specific baseline hazards and shares measured covariate effects. Alternative landmark models, proportional-hazards diagnostics, expanded-comparator analyses, and a smooth scaled-Schoenfeld-residual display assess model and time-interval sensitivity.

## Data governance

SEER data are governed by the SEER Research Data Agreement and cannot be redistributed here. The `.gitignore` excludes common raw-data, patient-level, result, archive, image, and manuscript formats. Users remain responsible for complying with their own SEER data agreement.
