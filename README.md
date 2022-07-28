# Group-3-Pulmonary-Hypertension
Pulmonary hypertension (PH) is a life-threatening complication of interstitial lung disease (ILD-PH). We investigate whether treatment with phosphodiesterase 5 inhibitors (PDE5i) in patients with ILD-PH was associated with improved survival. We present a longitudinal cohort (n=128) of incident group 3 PH patients with over 220 patient.years' follow up to understand the associations of treatment with survival. Code used in the analysis is given below.

## Software and hardware requirements
All statistical analysis was performed using R version 4.0.2 (R Foundation for Statistical Computing, Vienna, Austria; www.r-project.org) using RStudio version 1.3.1073 (Boston, Mass).

## Files
### 1. Functions.R
#### Packages and hand-written functions required for analysis
Loads all libraries and bespoke functions necessary for the code to run. Some functions relate to Bayesian analysis, some to pooling of the analysed imputed datasets and some are merely convenience for the write-up.

### 2. Setup.R
#### Imputation of Missing Data
Uploads the clinical data from an Excel sheet (not included) into R, reformats into longitudinal time-stamped format and imputes to create multiply imputed final datasets for analysis.

### 3. Bayesian Survival Model.R
#### Bayesian Survival Analysis
Survival analysis to understand the associations of treatment with other clinical surrogates and survival.

### 4. Bayesian LMM.R
#### Bayesian linear mixed-models
Mixed model analysis with pooling of datasets to understand the association of treatment with haemodynamic, echocardiographic and functional surrogates of clinical status.

### 5. Figure x.R
#### R code to draw the figures in the manuscript

### 6. Supplementary Figure x.R
#### R code to draw the supplementary figures in the manuscript

### Reference:
Dawes et al, Phosphodiesterase 5 inhibitor treatment and survival in interstitial lung disease pulmonary hypertension: A Bayesian retrospective observational cohort study



