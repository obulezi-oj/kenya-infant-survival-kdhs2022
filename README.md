# Born Unequal: Structural Stratification, Institutional Exclusion, and Infant Survival Across Kenyan Counties

This repository contains the complete R analysis pipeline for evaluating geographic deprivation, non-facility childbirth, and infant mortality ($<12$ months) across all 47 Kenyan counties using the 2022 Kenya Demographic and Health Survey (KDHS).

## Overview

Analyzing $N = 77,381$ unweighted live births ($70,933.9$ survey-weighted births) from the 2022 KDHS Births Recode (`KEBR8CFL`), this project combines complex survey econometrics with cluster-grouped supervised and unsupervised learning across four stages:

1. **Demographic & County Classification (Part A):** Computes survey-weighted sample distributions, DHS synthetic-cohort childhood mortality rates (NNMR, PNNMR, IMR, CMR, U5MR), a 47-county deprivation table, and a $k$-means / Principal Component Analysis (PCA) classification grouping counties into four structural vulnerability regimes.
2. **Design-Based Survey Regressions (Part B):** Estimates Taylor-linearized quasibinomial logistic regressions (`svyglm`) accounting for DHS weights, enumeration area clustering (`v001`), and urban–rural county strata (`v022`) to model non-facility delivery and infant mortality.
3. **Cluster-Grouped Predictive Modeling & Stacking (Part C):** Partitions the 12–59 month survival cohort ($75\%/25\%$) at the Primary Sampling Unit (PSU) level to prevent cluster leakage, benchmarking five survey-weighted models—Logistic Regression, Elastic Net, Random Forest, XGBoost, and a three-layer Deep Multi-Layer Perceptron ($64 \rightarrow 32 \rightarrow 16 \rightarrow 1$)—combined via a binomial Stacked Super-Learner Ensemble.
4. **Model Interpretation & Risk Concentration (Part D):** Generates out-of-sample permutation feature importance, counterfactual mortality predictions across a $40$-cell intersectional grid ($\text{Residence} \times \text{Wealth} \times \text{Education}$), maternal age partial dependence curves, and risk-decile concentration metrics.

## Data Access & Privacy Notice

Raw microdata files from the 2022 KDHS (`KEBR8CFL.DTA`, `.SAV`, or `.DAT`) are not included in this repository in accordance with Demographic and Health Surveys (DHS) Program data-use terms.

* To reproduce the analysis, register and download the **2022 Kenya DHS Births Recode (`KEBR8CFL`)** from [The DHS Program](https://dhsprogram.com/).
* Place the downloaded `.zip` archive or extracted dataset files inside your local working directory; the script automatically detects and loads `.dta`, `.sav`, or `.dat`/`.dct` formats.
* Add `*.dta`, `*.sav`, `*.dat`, `*.zip`, and `extracted_data/` to your `.gitignore` file before pushing changes.

## Requirements

The script automatically checks for and installs any missing CRAN packages on first run:

* **Data & Survey Design:** `haven`, `tidyverse`, `survey`, `srvyr`, `DHS.rates`, `broom`, `scales`, `forcats`
* **Statistical & Machine Learning:** `glmnet`, `ranger`, `xgboost`

## Running the Pipeline

1. Clone this repository to your local machine.
2. Place the `KEBR8CFL` dataset file(s) in the project folder.
3. Update `work_dir` at the top of the R script to match your local directory path.
4. Run the script end-to-end (`Rscript analysis.R` or via RStudio). All tables and 300-DPI figures are saved automatically to the `outputs/` directory.

## Generated Outputs

### Tables (`outputs/*.csv`)
* `Table01_Sample_Characteristics_Weighted.csv` — Unweighted and survey-weighted sample breakdown.
* `Table02_National_and_Subgroup_Mortality.csv` — Five-year synthetic-cohort mortality rates by wealth and residence.
* `Table03_County_Structural_Inequality_47.csv` — Survey-weighted deprivation, care exclusion, and IMR across all 47 counties.
* `Table04_County_Structural_Clusters_ML.csv` — PCA scores and $k$-means structural regime assignments.
* `Table05_Survey_Logistic_HomeBirth_aOR.csv` — Adjusted odds ratios for non-facility home delivery.
* `Table06_Survey_Logistic_InfantDeath_aOR.csv` — Adjusted odds ratios for infant mortality ($<12$ months).
* `Table07_ML_DL_Ensemble_Model_Comparison.csv` — Out-of-sample ROC-AUC, LogLoss, Brier Score, and top-decile death capture.
* `Table08_Ensemble_Feature_Importance_XAI.csv` — Cross-model permutation feature importance and consensus scores.
* `Table09_Intersectional_Risk_Matrix.csv` — Stacked Ensemble predicted IMR across 40 residence, wealth, and education strata.
* `Table10_Risk_Decile_Calibration_Concentration.csv` — Observed vs. predicted mortality and cumulative death capture by risk decile.

### Figures (`outputs/*.png`)
* `Figure01_Class_Gradient_Exclusion_Mortality.png` — Non-facility delivery rates across household wealth quintiles.
* `Figure02_County_Divide_Top_vs_Bottom.png` — Comparison of the 10 highest-exclusion vs. 10 highest-access counties.
* `Figure03_County_Structural_Typology_PCA.png` — PCA projection of Kenya's 47 counties colored by $k$-means regime.
* `Figure04_Cluster_Structural_Profiles.png` — Heatmap of mean deprivation and mortality indicators across county regimes.
* `Figure05_Forest_Plot_Adjusted_Odds_Ratios.png` — Survey-weighted adjusted odds ratios ($95\%$ CI) for home birth and infant death.
* `Figure06_ML_DL_Ensemble_ROC_Curves.png` — Out-of-sample weighted ROC curves across all six predictive models.
* `Figure07_XAI_Feature_Importance_Comparison.png` — Normalized permutation feature importance across Random Forest, XGBoost, and Deep MLP.
* `Figure08_Intersectional_Heat_Matrix.png` — Predicted IMR across urban/rural residence, wealth, and maternal education.
* `Figure09_LifeCourse_NonLinear_Risk_Curve.png` — Partial dependence curve of infant mortality across maternal age (15–44 years).
* `Figure10_Criminological_Risk_Concentration_Lift.png` — Cumulative share of weighted infant deaths captured across predicted risk deciles.
