# Manuscript 1 — Psychometric Validation of the EVAPSI Stress Scale

## Overview

This folder contains the full R analysis pipeline for the first EVAPSI manuscript, which covers:

- **Multivariate normality** testing (Mardia's test)
- **Confirmatory Factor Analysis (CFA)** of the original 7-factor model
- **Exploratory Graph Analysis (EGA)** with bootstrapping and dimension stability
- **Exploratory Factor Analysis (EFA)** using polychoric correlations and oblimin rotation
- **CFA of the EFA-derived model** on a held-out half-sample (cross-validation)
- **Model comparison** (original vs. EFA-derived structure)
- **Reliability** (omega) for both models
- **Measurement invariance (MI)** testing across:
  - Biological sex (`sexo_bio`)
  - Gender identity (`id_genero`)
  - Study level (`study_level`)
- **Factor score extraction** and **Spearman correlations** with DASS-21 and PWB subscales

---

## File

| File | Description |
|------|-------------|
| `evapsi_stress.R` | Full analysis script (sequential, top to bottom) |

---

## Data Required

Place the following file in this folder (or adjust the path in the script) before running:

| File | Description |
|------|-------------|
| `evapsi_data.csv` | Main dataset. Columns 7–35: EVAPSI items (EU1–EU29); columns 36–39: DASS-21 subscales; columns 40–45: PWB subscales. Grouping variables: `sexo_bio`, `id_genero`, `study_level`, `facultad`, `edad`. |

---

## Analysis Workflow

```
1. Descriptive table (table1)
2. Split-half sampling (p1 = EFA sample, p2 = CFA sample)
3. Multivariate normality (MVN)
4. CFA — original model (fit_og) on full sample
5. UVA — unique variable analysis (EGAnet)
6. EGA + bootEGA on p1
7. EFA (oblimin, 7 factors, polychoric) on p1
8. CFA — EFA-derived model (fit_efa) on p2
9. Model comparison + reliability
10. Measurement invariance — biological sex
    └── Configural → Partial thresholds → Partial loadings
11. Measurement invariance — gender identity
    └── Configural → Partial thresholds → Partial loadings
12. Measurement invariance — study level
    └── Configural → Thresholds → Loadings (full)
13. Factor score extraction (EBM method)
14. Spearman correlations: EU×EU, EU×DASS, EU×PWB
```

---

## Key Decisions & Notes

- **Estimator:** WLSM / WLSMV (robust, for ordered categorical data)
- **Parameterization:** Delta (for MI testing with `measEq.syntax`)
- **Identification:** Wu & Estabrook (2016) for categorical MI
- **Partial invariance search:** custom `search_ninv()` function — iteratively frees the parameter with the largest LRT chi-square difference
- **Non-invariant thresholds freed:**
  - Sex: EU16.thr1, EU3.thr3
  - Gender: EU16.thr2, EU3.thr1
- **Non-invariant loadings freed:**
  - Sex: lambda.20_5, lambda.9_3, lambda.15_4
  - Gender: lambda.3_1, lambda.20_5, lambda.9_3
- **Study level:** Full threshold and loading invariance supported

---

## Required Packages

```r
install.packages(c("lavaan", "semTools", "EGAnet", "psych",
                   "MVN", "dplyr", "table1", "knitr", "readxl"))
```

---

## Output

The script produces:

- Fit index tables (`results_og`, `results_efa`, `comp_models`)
- Standardized factor loadings
- MI fit tables (`MI_sex`, `MI_gen`, `MI_lvl`)
- Spearman correlation data frames (`EU_EU_results`, `EU_DASS_results`, `EU_PWB_results`)
- Saved factor scores: `EU_scores` (RDS file)
