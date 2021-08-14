# Group-3-Pulmonary-Hypertension
Pulmonary hypertension (PH) in chronic lung disease is common, linked to poor prognosis and has no long-term data to understand treatment strategies. We present longitudinal cohort (n=183) of incident group 3 PH patients with over 380 patient.years follow up to understand the associations of treatment with survival. Code used in the analysis is given below.

## Software and hardware requirements
All statistical analysis was performed using R version 4.0.2 (R Foundation for Statistical Computing, Vienna, Austria; www.r-project.org) using RStudio version 1.3.1073 (Boston, Mass). Analysis was completed on a 2017 iMac with 4.2 GHz Quad-Core Intel Core i7 and 64Gb RAM.

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

### 5. (Supplementary) Figure x.R
#### R code to draw the figures in the manuscript

### Reference:
Dawes et al, Pulmonary vasodilator treatment and survival in group 3 pulmonary hypertension: an observational study



