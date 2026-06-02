# =============================================================================
# EVAPSI — Manuscript 1
# Psychometric Validation and Measurement Invariance of the Academic Stress Scale
#
# Authors:   EVAPSI Research Group
# Affiliation: Universidad Pública de Navarra (UPNA)
# Last updated: 2026
#
# Description:
#   Full analysis pipeline including:
#     1. Descriptive statistics
#     2. Multivariate normality
#     3. Confirmatory Factor Analysis (CFA) — original model
#     4. Exploratory Graph Analysis (EGA) with bootstrapping
#     5. Exploratory Factor Analysis (EFA) — polychoric + oblimin
#     6. CFA — EFA-derived model (cross-validation)
#     7. Model comparison and reliability
#     8. Measurement invariance (MI) by sex, gender identity, and study level
#     9. Factor score extraction and convergent validity correlations
#
# Data:  evapsi_data.csv
#        Columns 7–35:  EVAPSI items (EU1–EU29)
#        Columns 36–39: DASS-21 subscales
#        Columns 40–45: PWB subscales
#
# Estimator: WLSMV (robust, ordered-categorical data)
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(lavaan)
library(MVN)
library(semTools)
library(knitr)
library(EGAnet)
library(psych)
library(table1)


# -----------------------------------------------------------------------------
# 1. HELPER FUNCTION — Iterative search for non-invariant parameters
# -----------------------------------------------------------------------------
#
# search_ninv() identifies the single constrained parameter (loading or
# threshold) whose release produces the largest chi-square improvement in a
# likelihood-ratio test (Satorra & Bentler, 2000). It is called iteratively
# after each partial-invariance step.
#
# Arguments:
#   data  : data frame used in the CFA
#   model : lavaan model syntax string (with equality constraints)
#   fit   : fitted lavaan object corresponding to `model`
#   group : character; grouping variable name
#   par   : "loadings" or "thresholds"
#
# Returns a list with:
#   chi_diff     : vector of chi-square differences for each candidate
#   max_chi_diff : maximum chi-square difference
#   label        : label of the most non-invariant parameter

search_ninv <- function(data, model, fit,
                        group,
                        par = c("loadings", "thresholds")) {

  par <- match.arg(par)
  pe  <- parameterestimates(fit)

  # Extract candidate parameter labels from group 1
  if (par == "loadings") {
    all_labels <- subset(pe, op == "=~" & group == 1)$label
  } else {
    all_labels <- subset(pe, op == "|"  & group == 1)$label
  }

  # Keep only labels still constrained in the model syntax
  labels <- all_labels[
    sapply(all_labels, function(lab) {
      grepl(paste0("c\\(", lab, ",\\s*", lab, "\\)"), model)
    })
  ]

  chi_diff <- rep(NA_real_, length(labels))

  for (i in seq_along(labels)) {

    pattern     <- paste0("c\\(", labels[i], ",\\s*", labels[i], "\\)")
    partial_mod <- gsub(pattern, "c(NA, NA)", model, perl = TRUE)

    if (identical(partial_mod, model)) next   # no change — skip

    partial_fit <- cfa(
      model            = partial_mod,
      data             = data,
      estimator        = "WLSMV",
      ordered          = TRUE,
      group            = group,
      parameterization = "delta"
    )

    lrt <- tryCatch(
      lavTestLRT(partial_fit, fit, method = "satorra.2000"),
      error = function(e) NULL
    )

    if (is.null(lrt)) next

    df_diff <- lrt["Df diff"][[1]][2]
    if (df_diff == 0) next

    chi_diff[i] <- lrt["Chisq diff"][[1]][2]
  }

  list(
    chi_diff     = chi_diff,
    max_chi_diff = max(chi_diff, na.rm = TRUE),
    label        = if (all(is.na(chi_diff))) NA else labels[which.max(chi_diff)]
  )
}


# -----------------------------------------------------------------------------
# 2. DATA
# -----------------------------------------------------------------------------

evapsi <- read.csv("evapsi_data.csv")


# -----------------------------------------------------------------------------
# 3. DESCRIPTIVE STATISTICS
# -----------------------------------------------------------------------------

table1(
  ~ sexo_bio + id_genero + edad + study_level + facultad,
  data = evapsi
)


# -----------------------------------------------------------------------------
# 4. SPLIT-HALF SAMPLING (EFA sample = p1; CFA sample = p2)
# -----------------------------------------------------------------------------

set.seed(123)
s <- sample(nrow(evapsi))
f <- floor(nrow(evapsi) / 2)

evapsi_p1 <- evapsi[s[1:f],        7:35]
evapsi_p2 <- evapsi[s[-(1:f)],     7:35]

# Remove "EU" prefix for EGAnet / psych functions
colnames(evapsi_p1) <- gsub("EU", "", colnames(evapsi_p1))


# -----------------------------------------------------------------------------
# 5. MULTIVARIATE NORMALITY (Mardia's test)
# -----------------------------------------------------------------------------

normality <- mvn(evapsi[, 7:35], mvn_test = "mardia")
normality$multivariate_normality
kable(normality$descriptives)


# -----------------------------------------------------------------------------
# 6. CFA — ORIGINAL 7-FACTOR MODEL (full sample)
# -----------------------------------------------------------------------------

mod_og <- '
  F1 =~ EU1  + EU2  + EU3  + EU4  + EU5
  F2 =~ EU6  + EU7  + EU8  + EU9
  F3 =~ EU10 + EU11 + EU12 + EU13
  F4 =~ EU14 + EU15 + EU16 + EU17
  F5 =~ EU18 + EU19 + EU20 + EU21
  F6 =~ EU22 + EU23 + EU24 + EU25
  F7 =~ EU26 + EU27 + EU28 + EU29
'

fit_og <- cfa(mod_og, data = evapsi, ordered = TRUE, estimator = "WLSM")
summary(fit_og, fit.measures = TRUE)

# Fit indices
names_fit <- c(
  "chisq.scaled", "df.scaled", "pvalue.scaled",
  "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled",
  "cfi.scaled", "tli.scaled", "srmr"
)
results_og <- round(data.matrix(fitmeasures(fit_og, fit.measures = names_fit)), 3)
print(results_og)

# Standardized factor loadings
std_all      <- standardizedSolution(fit_og)
std_loadings <- subset(std_all, op == "=~")
std_loadings[, c("lhs", "rhs", "est.std", "pvalue")]


# -----------------------------------------------------------------------------
# 7. UNIQUE VARIABLE ANALYSIS (UVA)
# -----------------------------------------------------------------------------

evapsi_uva <- UVA(data = evapsi_p1)
evapsi_uva$keep_remove


# -----------------------------------------------------------------------------
# 8. EXPLORATORY GRAPH ANALYSIS (EGA)
# -----------------------------------------------------------------------------

evapsi_ega <- EGA(evapsi_p1)
summary(evapsi_ega)

# Bootstrap EGA
evapsi_boot <- bootEGA(data = evapsi_p1, seed = 1)

# Visual comparison: empirical vs. bootstrap
evapsi_compare <- compare.EGA.plots(
  evapsi_ega, evapsi_boot,
  labels = c("Empirical", "Bootstrap")
)

summary(evapsi_boot)
dimensionStability(evapsi_boot)


# -----------------------------------------------------------------------------
# 9. EXPLORATORY FACTOR ANALYSIS (EFA)
#    Polychoric correlations, 7 factors, oblimin rotation, PA extraction
# -----------------------------------------------------------------------------

evapsi_p1_poly <- psych::polychoric(evapsi_p1)$rho
evapsi_efa     <- fa(evapsi_p1_poly, nfactors = 7, rotate = "oblimin", fm = "pa")
print(evapsi_efa)

# Items with loadings >= .40 per factor
efa_loadings     <- as.data.frame(unclass(evapsi_efa$loadings))
items_per_factor <- lapply(
  colnames(efa_loadings),
  function(f) rownames(efa_loadings)[abs(efa_loadings[[f]]) >= 0.4]
)
names(items_per_factor) <- colnames(efa_loadings)
items_per_factor


# -----------------------------------------------------------------------------
# 10. CFA — EFA-DERIVED MODEL (held-out p2 sample)
# -----------------------------------------------------------------------------

mod_efa <- '
  F1 =~ EU1  + EU2  + EU3  + EU4  + EU5
  F2 =~ EU6  + EU7  + EU8
  F3 =~ EU9  + EU10 + EU11 + EU12 + EU13
  F4 =~ EU14 + EU15 + EU16 + EU17
  F5 =~ EU18 + EU19 + EU20
  F6 =~ EU22 + EU23 + EU24
  F7 =~ EU26 + EU27 + EU28 + EU29
'

fit_efa <- cfa(mod_efa, data = evapsi_p2, ordered = TRUE, estimator = "WLSM")
summary(fit_efa, fit.measures = TRUE)

results_efa <- round(data.matrix(fitmeasures(fit_efa, fit.measures = names_fit)), 3)
print(results_efa)

# Standardized loadings
std_all      <- standardizedSolution(fit_efa)
std_loadings <- subset(std_all, op == "=~")
std_loadings[, c("lhs", "rhs", "est.std", "pvalue")]


# -----------------------------------------------------------------------------
# 11. MODEL COMPARISON
# -----------------------------------------------------------------------------

lavTestLRT(fit_og, fit_efa)

comp_models        <- matrix(NA, nrow = 2, ncol = 9)
comp_models[1, ]   <- round(data.matrix(fitmeasures(fit_og,  fit.measures = names_fit)), 3)
comp_models[2, ]   <- round(data.matrix(fitmeasures(fit_efa, fit.measures = names_fit)), 3)
colnames(comp_models) <- c("chisq", "df", "pvalue", "rmsea",
                           "rmsea.ci.lower", "rmsea.ci.upper",
                           "cfi", "tli", "srmr")
rownames(comp_models) <- c("Original model", "EFA-derived model")
print(comp_models)


# -----------------------------------------------------------------------------
# 12. RELIABILITY (omega)
# -----------------------------------------------------------------------------

models <- list(
  model_og  = semTools::reliability(fit_og),
  model_efa = semTools::reliability(fit_efa)
)
for (name in names(models)) {
  cat("\n", strrep("=", 40), "\n")
  cat("  ", name, "\n")
  cat(strrep("=", 40), "\n")
  print(models[[name]])
}


# =============================================================================
# 13. MEASUREMENT INVARIANCE (MI)
#     Approach: Wu & Estabrook (2016), delta parameterization
#     Sequence: configural → partial thresholds → partial loadings
#     Non-invariant parameters freed iteratively using search_ninv()
# =============================================================================

# Redefine mod_efa on full sample for MI analyses
mod_efa <- '
  F1 =~ EU1  + EU2  + EU3  + EU4  + EU5
  F2 =~ EU6  + EU7  + EU8
  F3 =~ EU9  + EU10 + EU11 + EU12 + EU13
  F4 =~ EU14 + EU15 + EU16 + EU17
  F5 =~ EU18 + EU19 + EU20
  F6 =~ EU22 + EU23 + EU24
  F7 =~ EU26 + EU27 + EU28 + EU29
'


# --- 13a. MI by Biological Sex -----------------------------------------------

# Exclude "Otro" (n too small for two-group CFA)
table(evapsi$sexo_bio)
evapsi_sex <- evapsi[evapsi$sexo_bio != "Otro", ]


# Configural model
mod_conf_sex <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi_sex,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "sexo_bio",
    group.equal      = "configural"
  )
)

fit_conf_sex <- cfa(
  mod_conf_sex,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
summary(fit_conf_sex, fit.measures = TRUE)


# Threshold invariance model
mod_thres_sex <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi_sex,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "sexo_bio",
    group.equal      = "thresholds"
  )
)

fit_thres_sex <- cfa(
  mod_thres_sex,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
summary(fit_thres_sex, fit.measures = TRUE)
lavTestLRT(fit_conf_sex, fit_thres_sex)   # Full thresholds significantly worse


# Iterative partial threshold search
# sex_thres <- search_ninv(
#   data  = evapsi_sex,
#   model = mod_thres_sex,
#   fit   = fit_thres_sex,
#   group = "sexo_bio",
#   par   = "thresholds"
# )
# sex_thres   # -> EU16.thr1

# Free EU16.thr1
sex_thres_flagged1 <- paste0("c\\(", "EU16.thr1", ", ", "EU16.thr1", "\\)")
mod_thres_sex2     <- gsub(sex_thres_flagged1, "c(NA, NA)", mod_thres_sex, perl = TRUE)

fit_thres_sex2 <- cfa(
  mod_thres_sex2,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
lavTestLRT(fit_conf_sex, fit_thres_sex2)


# sex_thres2 <- search_ninv(
#   data  = evapsi_sex,
#   model = mod_thres_sex2,
#   fit   = fit_thres_sex2,
#   group = "sexo_bio",
#   par   = "thresholds"
# )
# sex_thres2   # -> EU3.thr3

# Free EU3.thr3
sex_thres_flagged2 <- paste0("c\\(", "EU3.thr3", ", ", "EU3.thr3", "\\)")
mod_thres_sex3     <- gsub(sex_thres_flagged2, "c(NA, NA)", mod_thres_sex2, perl = TRUE)

fit_thres_sex3 <- cfa(
  mod_thres_sex3,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
lavTestLRT(fit_conf_sex, fit_thres_sex3)   # Partial threshold model accepted


# Metric invariance (loadings), carrying forward freed thresholds
mod_load_sex <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi_sex,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "sexo_bio",
    group.equal      = c("thresholds", "loadings")
  )
)
mod_load_sex <- gsub(sex_thres_flagged1, "c(NA, NA)", mod_load_sex, perl = TRUE)
mod_load_sex <- gsub(sex_thres_flagged2, "c(NA, NA)", mod_load_sex, perl = TRUE)

fit_load_sex <- cfa(
  mod_load_sex,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
summary(fit_load_sex, fit.measures = TRUE)
lavTestLRT(fit_thres_sex3, fit_load_sex)   # Full loadings significantly worse


# Iterative partial loading search
# sex_load <- search_ninv(
#   data  = evapsi_sex,
#   model = mod_load_sex,
#   fit   = fit_load_sex,
#   group = "sexo_bio",
#   par   = "loadings"
# )
# sex_load   # -> lambda.20_5

sex_load_flagged1 <- "c\\(lambda\\.20_5, lambda\\.20_5\\)"
mod_load_sex2     <- gsub(sex_load_flagged1, "c(NA, NA)", mod_load_sex, perl = TRUE)

fit_load_sex2 <- cfa(
  mod_load_sex2,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
lavTestLRT(fit_thres_sex3, fit_load_sex2)


# sex_load2 <- search_ninv(
#   data  = evapsi_sex,
#   model = mod_load_sex2,
#   fit   = fit_load_sex2,
#   group = "sexo_bio",
#   par   = "loadings"
# )
# sex_load2   # -> lambda.9_3

sex_load_flagged2 <- "c\\(lambda\\.9_3, lambda\\.9_3\\)"
mod_load_sex3     <- gsub(sex_load_flagged2, "c(NA, NA)", mod_load_sex2, perl = TRUE)

fit_load_sex3 <- cfa(
  mod_load_sex3,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
lavTestLRT(fit_thres_sex3, fit_load_sex3)


# sex_load3 <- search_ninv(
#   data  = evapsi_sex,
#   model = mod_load_sex3,
#   fit   = fit_load_sex3,
#   group = "sexo_bio",
#   par   = "loadings"
# )
# sex_load3   # -> lambda.15_4

sex_load_flagged3 <- "c\\(lambda\\.15_4, lambda\\.15_4\\)"
mod_load_sex4     <- gsub(sex_load_flagged3, "c(NA, NA)", mod_load_sex3, perl = TRUE)

fit_load_sex4 <- cfa(
  mod_load_sex4,
  data    = evapsi_sex,
  ordered = TRUE,
  group   = "sexo_bio"
)
lavTestLRT(fit_thres_sex3, fit_load_sex4)   # Partial loading model accepted


# MI fit summary table — Biological Sex
MI_sex        <- matrix(NA, nrow = 3, ncol = 9)
MI_sex[1, ]   <- round(data.matrix(fitmeasures(fit_conf_sex,   fit.measures = names_fit)), 3)
MI_sex[2, ]   <- round(data.matrix(fitmeasures(fit_thres_sex3, fit.measures = names_fit)), 3)
MI_sex[3, ]   <- round(data.matrix(fitmeasures(fit_load_sex4,  fit.measures = names_fit)), 3)
colnames(MI_sex) <- c("chisq", "df", "pvalue", "rmsea",
                      "rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(MI_sex) <- c("Configural", "Partial thresholds", "Partial loadings")
print(MI_sex)


# --- 13b. MI by Gender Identity ----------------------------------------------

table(evapsi$id_genero)
evapsi_gen <- evapsi[evapsi$id_genero != "Otro", ]


# Configural model
mod_conf_gen <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi_gen,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "id_genero",
    group.equal      = "configural"
  )
)

fit_conf_gen <- cfa(
  mod_conf_gen,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
summary(fit_conf_gen, fit.measures = TRUE)


# Threshold invariance model
mod_thres_gen <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi_gen,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "id_genero",
    group.equal      = "thresholds"
  )
)

fit_thres_gen <- cfa(
  mod_thres_gen,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
summary(fit_thres_gen, fit.measures = TRUE)
lavTestLRT(fit_conf_gen, fit_thres_gen)


# Iterative partial threshold search
# gen_thres <- search_ninv(
#   data  = evapsi_gen,
#   model = mod_thres_gen,
#   fit   = fit_thres_gen,
#   group = "id_genero",
#   par   = "thresholds"
# )
# gen_thres   # -> EU16.thr2

gen_thres_flagged1 <- paste0("c\\(", "EU16.thr2", ", ", "EU16.thr2", "\\)")
mod_thres_gen2     <- gsub(gen_thres_flagged1, "c(NA, NA)", mod_thres_gen, perl = TRUE)

fit_thres_gen2 <- cfa(
  mod_thres_gen2,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
lavTestLRT(fit_conf_gen, fit_thres_gen2)


# gen_thres2 <- search_ninv(
#   data  = evapsi_gen,
#   model = mod_thres_gen2,
#   fit   = fit_thres_gen2,
#   group = "id_genero",
#   par   = "thresholds"
# )
# gen_thres2   # -> EU3.thr1

gen_thres_flagged2 <- paste0("c\\(", "EU3.thr1", ", ", "EU3.thr1", "\\)")
mod_thres_gen3     <- gsub(gen_thres_flagged2, "c(NA, NA)", mod_thres_gen2, perl = TRUE)

fit_thres_gen3 <- cfa(
  mod_thres_gen3,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
lavTestLRT(fit_conf_gen, fit_thres_gen3)   # Partial threshold model accepted


# Metric invariance (loadings)
mod_load_gen <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi_gen,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "id_genero",
    group.equal      = c("thresholds", "loadings")
  )
)
mod_load_gen <- gsub(gen_thres_flagged1, "c(NA, NA)", mod_load_gen, perl = TRUE)
mod_load_gen <- gsub(gen_thres_flagged2, "c(NA, NA)", mod_load_gen, perl = TRUE)

fit_load_gen <- cfa(
  mod_load_gen,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
summary(fit_load_gen, fit.measures = TRUE)
lavTestLRT(fit_thres_gen3, fit_load_gen)


# Iterative partial loading search
# gen_load <- search_ninv(
#   data  = evapsi_gen,
#   model = mod_load_gen,
#   fit   = fit_load_gen,
#   group = "id_genero",
#   par   = "loadings"
# )
# gen_load   # -> lambda.3_1

gen_load_flagged1 <- "c\\(lambda\\.3_1, lambda\\.3_1\\)"
mod_load_gen2     <- gsub(gen_load_flagged1, "c(NA, NA)", mod_load_gen, perl = TRUE)

fit_load_gen2 <- cfa(
  mod_load_gen2,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
lavTestLRT(fit_thres_gen3, fit_load_gen2)


# gen_load2 <- search_ninv(
#   data  = evapsi_gen,
#   model = mod_load_gen2,
#   fit   = fit_load_gen2,
#   group = "id_genero",
#   par   = "loadings"
# )
# gen_load2   # -> lambda.20_5

gen_load_flagged2 <- "c\\(lambda\\.20_5, lambda\\.20_5\\)"
mod_load_gen3     <- gsub(gen_load_flagged2, "c(NA, NA)", mod_load_gen2, perl = TRUE)

fit_load_gen3 <- cfa(
  mod_load_gen3,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
lavTestLRT(fit_thres_gen3, fit_load_gen3)


# gen_load3 <- search_ninv(
#   data  = evapsi_gen,
#   model = mod_load_gen3,
#   fit   = fit_load_gen3,
#   group = "id_genero",
#   par   = "loadings"
# )
# gen_load3   # -> lambda.9_3

gen_load_flagged3 <- "c\\(lambda\\.9_3, lambda\\.9_3\\)"
mod_load_gen4     <- gsub(gen_load_flagged3, "c(NA, NA)", mod_load_gen3, perl = TRUE)

fit_load_gen4 <- cfa(
  mod_load_gen4,
  data    = evapsi_gen,
  ordered = TRUE,
  group   = "id_genero"
)
lavTestLRT(fit_thres_gen3, fit_load_gen4)   # Partial loading model accepted


# MI fit summary table — Gender Identity
MI_gen        <- matrix(NA, nrow = 3, ncol = 9)
MI_gen[1, ]   <- round(data.matrix(fitmeasures(fit_conf_gen,   fit.measures = names_fit)), 3)
MI_gen[2, ]   <- round(data.matrix(fitmeasures(fit_thres_gen3, fit.measures = names_fit)), 3)
MI_gen[3, ]   <- round(data.matrix(fitmeasures(fit_load_gen4,  fit.measures = names_fit)), 3)
colnames(MI_gen) <- c("chisq", "df", "pvalue", "rmsea",
                      "rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(MI_gen) <- c("Configural", "Partial thresholds", "Partial loadings")
print(MI_gen)


# --- 13c. MI by Study Level --------------------------------------------------

table(evapsi$study_level)


# Configural model
mod_conf_lvl <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "study_level",
    group.equal      = "configural"
  )
)

fit_conf_lvl <- cfa(
  mod_conf_lvl,
  data    = evapsi,
  ordered = TRUE,
  group   = "study_level"
)
summary(fit_conf_lvl, fit.measures = TRUE)


# Threshold invariance model
mod_thres_lvl <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "study_level",
    group.equal      = "thresholds"
  )
)

fit_thres_lvl <- cfa(
  mod_thres_lvl,
  data    = evapsi,
  ordered = TRUE,
  group   = "study_level"
)
summary(fit_thres_lvl, fit.measures = TRUE)
lavTestLRT(fit_conf_lvl, fit_thres_lvl)   # Full thresholds accepted


# Metric invariance (loadings)
mod_load_lvl <- as.character(
  measEq.syntax(
    mod_efa,
    data             = evapsi,
    ordered          = TRUE,
    parameterization = "delta",
    ID.fac           = "std.lv",
    ID.cat           = "Wu.Estabrook.2016",
    group            = "study_level",
    group.equal      = c("thresholds", "loadings")
  )
)

fit_load_lvl <- cfa(
  mod_load_lvl,
  data    = evapsi,
  ordered = TRUE,
  group   = "study_level"
)
summary(fit_load_lvl, fit.measures = TRUE)
lavTestLRT(fit_thres_lvl, fit_load_lvl)   # Full loadings accepted


# MI fit summary table — Study Level
MI_lvl        <- matrix(NA, nrow = 3, ncol = 9)
MI_lvl[1, ]   <- round(data.matrix(fitmeasures(fit_conf_lvl,  fit.measures = names_fit)), 3)
MI_lvl[2, ]   <- round(data.matrix(fitmeasures(fit_thres_lvl, fit.measures = names_fit)), 3)
MI_lvl[3, ]   <- round(data.matrix(fitmeasures(fit_load_lvl,  fit.measures = names_fit)), 3)
colnames(MI_lvl) <- c("chisq", "df", "pvalue", "rmsea",
                      "rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(MI_lvl) <- c("Configural", "Full thresholds", "Full loadings")
print(MI_lvl)


# =============================================================================
# 14. FACTOR SCORES & CONVERGENT VALIDITY
#     EBM (Empirical Bayes Modal) method via lavPredict()
# =============================================================================

# Fit EFA-derived model on full sample for factor score extraction
full_fit_efa <- cfa(mod_efa, data = evapsi, ordered = TRUE, estimator = "WLSM")

idx        <- lavInspect(full_fit_efa, "case.idx")
EU_Fscores <- lavPredict(full_fit_efa, type = "lv", method = "EBM")

# Append factor scores to main data frame
for (g in seq_along(EU_Fscores)) {
  for (fs in colnames(EU_Fscores[[g]])) {
    evapsi[idx[[g]], fs] <- EU_Fscores[[g]][, fs]
  }
}

EU_Fscores <- as.data.frame(EU_Fscores)
saveRDS(EU_Fscores, file = "EU_scores.rds")


# --- 14a. Inter-factor correlations (EU × EU) --------------------------------

EU_EU_pairs <- expand.grid(
  EU_factor1 = colnames(EU_Fscores),
  EU_factor2 = colnames(EU_Fscores),
  stringsAsFactors = FALSE
)

EU_EU_cortests <- apply(EU_EU_pairs, 1, function(x) {
  cor.test(
    EU_Fscores[[x["EU_factor1"]]],
    EU_Fscores[[x["EU_factor2"]]],
    method = "spearman", exact = FALSE
  )
})

EU_EU_results <- do.call(rbind, lapply(seq_along(EU_EU_cortests), function(i) {
  data.frame(
    EU_factor1 = EU_EU_pairs$EU_factor1[i],
    EU_factor2 = EU_EU_pairs$EU_factor2[i],
    rho        = unname(EU_EU_cortests[[i]]$estimate),
    p_value    = EU_EU_cortests[[i]]$p.value
  )
}))
print(EU_EU_results)


# --- 14b. Correlations with DASS-21 subscales (EU × DASS) -------------------

EU_DASS_pairs <- expand.grid(
  EU_factor   = colnames(EU_Fscores),
  DASS_factor = colnames(evapsi[, 36:39]),
  stringsAsFactors = FALSE
)

EU_DASS_cortests <- apply(EU_DASS_pairs, 1, function(x) {
  cor.test(
    EU_Fscores[[x["EU_factor"]]],
    evapsi[[x["DASS_factor"]]],
    method = "spearman", exact = FALSE
  )
})

EU_DASS_results <- do.call(rbind, lapply(seq_along(EU_DASS_cortests), function(i) {
  data.frame(
    EU_factor   = EU_DASS_pairs$EU_factor[i],
    DASS_factor = EU_DASS_pairs$DASS_factor[i],
    rho         = unname(EU_DASS_cortests[[i]]$estimate),
    p_value     = EU_DASS_cortests[[i]]$p.value
  )
}))
print(EU_DASS_results)


# --- 14c. Correlations with PWB subscales (EU × PWB) ------------------------

EU_PWB_pairs <- expand.grid(
  EU_factor  = colnames(EU_Fscores),
  PWB_factor = colnames(evapsi[, 40:45]),
  stringsAsFactors = FALSE
)

EU_PWB_cortests <- apply(EU_PWB_pairs, 1, function(x) {
  cor.test(
    EU_Fscores[[x["EU_factor"]]],
    evapsi[[x["PWB_factor"]]],
    method = "spearman", exact = FALSE
  )
})

EU_PWB_results <- do.call(rbind, lapply(seq_along(EU_PWB_cortests), function(i) {
  data.frame(
    EU_factor  = EU_PWB_pairs$EU_factor[i],
    PWB_factor = EU_PWB_pairs$PWB_factor[i],
    rho        = unname(EU_PWB_cortests[[i]]$estimate),
    p_value    = EU_PWB_cortests[[i]]$p.value
  )
}))
print(EU_PWB_results)

# =============================================================================
# END OF SCRIPT
# =============================================================================
