## ----setup, include=FALSE---------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----message=FALSE, warning=FALSE-------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(lavaan)
library(MVN)
library(semTools)
library(knitr)
library(EGAnet)
library(psych)
library(table1)


## ---------------------------------------------------------------------------------------------------------
search_ninv <- function(data, model, fit, group, 
                        par = c("loadings", "thresholds")) {

  par <- match.arg(par)
  pe  <- parameterestimates(fit)

  ## extract candidate labels from group 1
  if (par == "loadings") {
    all_labels <- subset(pe, op == "=~" & group == 1)$label
  } else if (par == "thresholds") {
    all_labels <- subset(pe, op == "|" & group == 1)$label
  }

  ## keep only labels that are STILL constrained in the model syntax
  labels <- all_labels[
    sapply(all_labels, function(lab) {
      grepl(paste0("c\\(", lab, ",\\s*", lab, "\\)"), model)
    })
  ]

  chi_diff <- rep(NA_real_, length(labels))

  for (i in seq_along(labels)) {

    pattern <- paste0("c\\(", labels[i], ",\\s*", labels[i], "\\)")
    partial_mod <- gsub(pattern, "c(NA, NA)", model, perl = TRUE)

    ## skip if freeing this parameter does not change the model
    if (identical(partial_mod, model)) next

    partial_fit <- cfa(
      model = partial_mod,
      data  = data,
      estimator = "WLSMV",
      ordered = TRUE,
      group = group,
      parameterization = "delta"
    )

    ## run LRT safely
    lrt <- tryCatch(
      lavTestLRT(partial_fit, fit, method = "satorra.2000"),
      error = function(e) NULL
    )

    ## skip failed or non-testable comparisons
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


## ---------------------------------------------------------------------------------------------------------
evapsi = read.csv("evapsi_data.csv")


## ---------------------------------------------------------------------------------------------------------
table1(~ sexo_bio +
         id_genero +
         edad +
         study_level +
         facultad,
       data = evapsi)


## ---------------------------------------------------------------------------------------------------------
set.seed(123)
s = sample(nrow(evapsi))
f = floor(nrow(evapsi)/2)

evapsi_p1 = evapsi[s[1:f],7:35]
evapsi_p2 = evapsi[s[-(1:f)],7:35]

colnames(evapsi_p1) <- gsub("EU", "", colnames(evapsi_p1))


## ---------------------------------------------------------------------------------------------------------
normality = mvn(evapsi[,7:35],   mvn_test = "mardia")
normality$multivariate_normality
kable(normality$descriptives)


## ---------------------------------------------------------------------------------------------------------
mod_og = 'F1 =~ EU1 + EU2 + EU3 + EU4 + EU5
          F2 =~ EU6 + EU7 + EU8 + EU9
          F3 =~ EU10 + EU11 + EU12 + EU13
          F4 =~ EU14 + EU15 + EU16 + EU17
          F5 =~ EU18 + EU19 + EU20 + EU21
          F6 =~ EU22 + EU23 + EU24 + EU25
          F7 =~ EU26 + EU27 + EU28 + EU29'


## ---------------------------------------------------------------------------------------------------------
fit_og = cfa(mod_og, data=evapsi, ordered=TRUE, estimator="WLSM")
check_og = summary(fit_og, fit.measures = TRUE)


## ---------------------------------------------------------------------------------------------------------
names_fit = c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "cfi.scaled", "tli.scaled", "srmr")
results_og = round(data.matrix(fitmeasures(fit_og, fit.measures = names_fit)), 3)
results_og


## ---------------------------------------------------------------------------------------------------------
std_all = standardizedSolution(fit_og)
std_loadings = subset(std_all, op=="=~")
std_loadings[, c("lhs", "rhs", "est.std", "pvalue")]


## ---------------------------------------------------------------------------------------------------------
evapsi_uva = UVA(data=evapsi_p1)
evapsi_uva$keep_remove


## ---------------------------------------------------------------------------------------------------------
evapsi_ega = EGA(evapsi_p1)


## ---------------------------------------------------------------------------------------------------------
summary(evapsi_ega)


## ---------------------------------------------------------------------------------------------------------
evapsi_boot = bootEGA(data=evapsi_p1, seed=1)


## ---------------------------------------------------------------------------------------------------------
evapsi_compare <- compare.EGA.plots(evapsi_ega, evapsi_boot,
                                    labels = c("Empirical", "Bootstrap"))


## ---------------------------------------------------------------------------------------------------------
summary(evapsi_boot)
dimensionStability(evapsi_boot)


## ---------------------------------------------------------------------------------------------------------
evapsi_p1_poly = psych::polychoric(evapsi_p1)$rho
evapsi_efa = fa(evapsi_p1_poly, nfactors=7, rotate="oblimin", fm="pa")
evapsi_efa


## ---------------------------------------------------------------------------------------------------------
efa_loadings = as.data.frame(unclass(evapsi_efa$loadings))

items_per_factor = lapply(colnames(efa_loadings), function(f) {rownames(efa_loadings)[abs(efa_loadings[[f]]) >= 0.4]})
names(items_per_factor) = colnames(efa_loadings)
items_per_factor


## ---------------------------------------------------------------------------------------------------------
mod_efa = 'F1 =~ EU1 + EU2 + EU3 + EU4 + EU5
           F2 =~ EU6 + EU7 + EU8
           F3 =~ EU9 + EU10 + EU11 + EU12 + EU13
           F4 =~ EU14 + EU15 + EU16 + EU17
           F5 =~ EU18 + EU19 + EU20
           F6 =~ EU22 + EU23 + EU24
           F7 =~ EU26 + EU27 + EU28 + EU29'


## ---------------------------------------------------------------------------------------------------------
fit_efa = cfa(mod_efa, data=evapsi_p2, ordered=TRUE, estimator="WLSM")
check_efa = summary(fit_efa, fit.measures = TRUE)


## ---------------------------------------------------------------------------------------------------------
names_fit = c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "cfi.scaled", "tli.scaled", "srmr")
results_efa = round(data.matrix(fitmeasures(fit_efa, fit.measures = names_fit)), 3)
results_efa


## ---------------------------------------------------------------------------------------------------------
std_all = standardizedSolution(fit_efa)
std_loadings = subset(std_all, op=="=~")
std_loadings[, c("lhs", "rhs", "est.std", "pvalue")]


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_og, fit_efa)


## ---------------------------------------------------------------------------------------------------------
# Prepare matrix
names_fit = c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "cfi.scaled", "tli.scaled", "srmr")
comp_models = matrix(NA, nrow = 2, ncol = 9)

# Add results to matrix
comp_models[1,] = round(data.matrix(fitmeasures(fit_og,
                                                fit.measures = names_fit)), 3)
comp_models[2,] = round(data.matrix(fitmeasures(fit_efa,
                                                fit.measures = names_fit)), 3)

# Print matrix
colnames(comp_models) = c("chisq","df","pvalue", "rmsea","rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(comp_models) = c("original", "from efa")
print(comp_models)


## ----message=FALSE----------------------------------------------------------------------------------------
models <- list(
  model_og = semTools::reliability(fit_og),
  model_efa = semTools::reliability(fit_efa))
for (name in names(models)) {
  cat("\n====================\n")
  cat("  ", name, "\n")
  cat("====================\n")
  print(models[[name]])}


## ---------------------------------------------------------------------------------------------------------
mod_efa = 'F1 =~ EU1 + EU2 + EU3 + EU4 + EU5
           F2 =~ EU6 + EU7 + EU8
           F3 =~ EU9 + EU10 + EU11 + EU12 + EU13
           F4 =~ EU14 + EU15 + EU16 + EU17
           F5 =~ EU18 + EU19 + EU20
           F6 =~ EU22 + EU23 + EU24
           F7 =~ EU26 + EU27 + EU28 + EU29'


## ---------------------------------------------------------------------------------------------------------
table(evapsi$sexo_bio)

evapsi_sex = evapsi[evapsi$sexo_bio!="Otro", ]


## ----results='hide'---------------------------------------------------------------------------------------
mod_conf_sex = as.character(measEq.syntax(mod_efa,
                                data = evapsi_sex,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "sexo_bio",
                                group.equal = "configural")) 
cat(as.character(mod_conf_sex))


## ----results='hide'---------------------------------------------------------------------------------------
fit_conf_sex = cfa(mod_conf_sex,  
                           data = evapsi_sex,            
                           ordered = T,             
                           group = "sexo_bio")
summary (fit_conf_sex, fit.measures = T)


## ----results='hide'---------------------------------------------------------------------------------------
mod_thres_sex = as.character(measEq.syntax(mod_efa,
                                data = evapsi_sex,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "sexo_bio",
                                group.equal = "thresholds")) 
cat(as.character(mod_thres_sex))


## ----results='hide'---------------------------------------------------------------------------------------
fit_thres_sex = cfa(mod_thres_sex,  
                           data = evapsi_sex,            
                           ordered = T,             
                           group = "sexo_bio")
summary (fit_thres_sex, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_sex, fit_thres_sex)


## ----eval=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------
## sex_thres = search_ninv(data = evapsi_sex,
##             model = mod_thres_sex,
##             fit = fit_thres_sex,
##             group = "sexo_bio",
##             par = "thresholds")
## sex_thres


## ----results='hide'---------------------------------------------------------------------------------------
sex_thres_flagged = paste0("c\\(", "EU16.thr1", ", ", "EU16.thr1", "\\)")

mod_thres_sex2 = gsub(sex_thres_flagged, "c(NA, NA)", mod_thres_sex, perl=TRUE)

cat(mod_thres_sex2)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_thres_sex2 = cfa(mod_thres_sex2,  
                           data = evapsi_sex,            
                           ordered = T,             
                           group = "sexo_bio")
summary (fit_thres_sex2, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_sex, fit_thres_sex2)


## ----eval=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------------
## sex_thres2 = search_ninv(data = evapsi_sex,
##             model = mod_thres_sex2,
##             fit = fit_thres_sex2,
##             group = "sexo_bio",
##             par = "thresholds")
## sex_thres2


## ----results='hide'---------------------------------------------------------------------------------------
sex_thres_flagged2 = paste0("c\\(", "EU3.thr3", ", ", "EU3.thr3", "\\)")

mod_thres_sex3 = gsub(sex_thres_flagged2, "c(NA, NA)", mod_thres_sex2, perl=TRUE)

cat(mod_thres_sex3)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_thres_sex3 = cfa(mod_thres_sex3,  
                           data = evapsi_sex,            
                           ordered = T,             
                           group = "sexo_bio")
summary (fit_thres_sex3, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_sex, fit_thres_sex3)


## ----results='hide'---------------------------------------------------------------------------------------
mod_load_sex = as.character(measEq.syntax(mod_efa,
                                data = evapsi_sex,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "sexo_bio",
                                group.equal = c("thresholds", "loadings"))) 
cat(as.character(mod_load_sex))


## ----results='hide'---------------------------------------------------------------------------------------
mod_load_sex = gsub(sex_thres_flagged, "c(NA, NA)", mod_load_sex, perl=TRUE)
mod_load_sex = gsub(sex_thres_flagged2, "c(NA, NA)", mod_load_sex, perl=TRUE)

cat(mod_load_sex)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_sex = cfa(mod_load_sex,  
                           data = evapsi_sex,            
                           ordered = T,             
                           group = "sexo_bio")
summary (fit_load_sex, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_sex3, fit_load_sex)


## ----eval=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------
## sex_load = search_ninv(data = evapsi_sex,
##             model = mod_load_sex,
##             fit = fit_load_sex,
##             group = "sexo_bio",
##             par = "loadings")
## sex_load


## ----results='hide'---------------------------------------------------------------------------------------
sex_load_flagged = "c\\(lambda\\.20_5, lambda\\.20_5\\)"

mod_load_sex2 = gsub(sex_load_flagged, "c(NA, NA)", mod_load_sex, perl=TRUE)

cat(mod_load_sex2)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_sex2 = cfa(mod_load_sex2,
                    data = evapsi_sex,  
                    ordered = T,  
                    group = "sexo_bio")
summary (fit_load_sex2, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_sex3, fit_load_sex2)


## ----eval=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------------
## sex_load2 = search_ninv(data = evapsi_sex,
##             model = mod_load_sex2,
##             fit = fit_load_sex2,
##             group = "sexo_bio",
##             par = "loadings")
## sex_load2


## ----results='hide'---------------------------------------------------------------------------------------
sex_load_flagged2 = "c\\(lambda\\.9_3, lambda\\.9_3\\)"

mod_load_sex3 = gsub(sex_load_flagged2, "c(NA, NA)", mod_load_sex2, perl=TRUE)

cat(mod_load_sex3)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_sex3 = cfa(mod_load_sex3,
                    data = evapsi_sex,  
                    ordered = T,  
                    group = "sexo_bio")
summary (fit_load_sex3, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_sex3, fit_load_sex3)


## ----eval=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------------
## sex_load3 = search_ninv(data = evapsi_sex,
##             model = mod_load_sex3,
##             fit = fit_load_sex3,
##             group = "sexo_bio",
##             par = "loadings")
## sex_load3


## ----results='hide'---------------------------------------------------------------------------------------
sex_load_flagged3 = "c\\(lambda\\.15_4, lambda\\.15_4\\)"

mod_load_sex4 = gsub(sex_load_flagged3, "c(NA, NA)", mod_load_sex3, perl=TRUE)

cat(mod_load_sex4)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_sex4 = cfa(mod_load_sex4,
                    data = evapsi_sex,  
                    ordered = T,  
                    group = "sexo_bio")
summary (fit_load_sex4, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_sex3, fit_load_sex4)


## ---------------------------------------------------------------------------------------------------------
MI_sex = matrix(NA, nrow = 3, ncol = 9)

MI_sex[1,] = round(data.matrix(fitmeasures(fit_conf_sex, fit.measures = names_fit)), 3)
MI_sex[2,] = round(data.matrix(fitmeasures(fit_thres_sex3, fit.measures = names_fit)), 3)
MI_sex[3,] = round(data.matrix(fitmeasures(fit_load_sex4, fit.measures = names_fit)), 3)


## ---------------------------------------------------------------------------------------------------------
colnames(MI_sex) = c("chisq","df","pvalue", "rmsea","rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(MI_sex) = c("configural", "partial thresholds", "partial loadings")
print(MI_sex)


## ---------------------------------------------------------------------------------------------------------
table(evapsi$id_genero)

evapsi_gen = evapsi[evapsi$id_genero!="Otro", ]


## ----results='hide'---------------------------------------------------------------------------------------
mod_conf_gen = as.character(measEq.syntax(mod_efa,
                                data = evapsi_gen,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "id_genero",
                                group.equal = "configural")) 
cat(as.character(mod_conf_gen))


## ----results='hide'---------------------------------------------------------------------------------------
fit_conf_gen = cfa(mod_conf_gen,  
                           data = evapsi_gen,            
                           ordered = T,             
                           group = "id_genero")
summary (fit_conf_gen, fit.measures = T)


## ----results='hide'---------------------------------------------------------------------------------------
mod_thres_gen = as.character(measEq.syntax(mod_efa,
                                data = evapsi_gen,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "id_genero",
                                group.equal = "thresholds")) 
cat(as.character(mod_thres_gen))


## ----results='hide'---------------------------------------------------------------------------------------
fit_thres_gen = cfa(mod_thres_gen,  
                           data = evapsi_gen,            
                           ordered = T,             
                           group = "id_genero")
summary (fit_thres_gen, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_gen, fit_thres_gen)


## ----eval=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------
## gen_thres = search_ninv(data = evapsi_gen,
##                         model = mod_thres_gen,
##                         fit = fit_thres_gen,
##                         group = "id_genero",
##                         par = "thresholds")
## gen_thres


## ----results='hide'---------------------------------------------------------------------------------------
gen_thres_flagged = paste0("c\\(", "EU16.thr2", ", ", "EU16.thr2", "\\)")

mod_thres_gen2 = gsub(gen_thres_flagged, "c(NA, NA)", mod_thres_gen, perl=TRUE)

cat(mod_thres_gen2)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_thres_gen2 = cfa(mod_thres_gen2,  
                           data = evapsi_gen,            
                           ordered = T,             
                           group = "id_genero")
summary (fit_thres_gen2, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_gen, fit_thres_gen2)


## ----eval=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------
## gen_thres2 = search_ninv(data = evapsi_gen,
##                         model = mod_thres_gen2,
##                         fit = fit_thres_gen2,
##                         group = "id_genero",
##                         par = "thresholds")
## gen_thres2


## ----results='hide'---------------------------------------------------------------------------------------
gen_thres_flagged2 = paste0("c\\(", "EU3.thr1", ", ", "EU3.thr1", "\\)")

mod_thres_gen3 = gsub(gen_thres_flagged2, "c(NA, NA)", mod_thres_gen2, perl=TRUE)

cat(mod_thres_gen3)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_thres_gen3 = cfa(mod_thres_gen3,  
                           data = evapsi_gen,            
                           ordered = T,             
                           group = "id_genero")
summary (fit_thres_gen3, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_gen, fit_thres_gen3)


## ----results='hide'---------------------------------------------------------------------------------------
mod_load_gen = as.character(measEq.syntax(mod_efa,
                                data = evapsi_gen,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "id_genero",
                                group.equal = c("thresholds", "loadings"))) 
cat(as.character(mod_load_gen))


## ----results='hide'---------------------------------------------------------------------------------------
mod_load_gen = gsub(gen_thres_flagged, "c(NA, NA)", mod_load_gen, perl=TRUE)
mod_load_gen = gsub(gen_thres_flagged2, "c(NA, NA)", mod_load_gen, perl=TRUE)

cat(mod_load_gen)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_gen = cfa(mod_load_gen,  
                           data = evapsi_gen,            
                           ordered = T,             
                           group = "id_genero")
summary (fit_load_gen, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_gen3, fit_load_gen)


## ----eval=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------
## gen_load = search_ninv(data = evapsi_gen,
##             model = mod_load_gen,
##             fit = fit_load_gen,
##             group = "id_genero",
##             par = "loadings")
## gen_load


## ----results='hide'---------------------------------------------------------------------------------------
gen_load_flagged = "c\\(lambda\\.3_1, lambda\\.3_1\\)"

mod_load_gen2 = gsub(gen_load_flagged, "c(NA, NA)", mod_load_gen, perl=TRUE)

cat(mod_load_gen2)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_gen2 = cfa(mod_load_gen2,
                    data = evapsi_gen,  
                    ordered = T,  
                    group = "id_genero")
summary (fit_load_gen2, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_gen3, fit_load_gen2)


## ----eval=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------------
## gen_load2 = search_ninv(data = evapsi_gen,
##             model = mod_load_gen2,
##             fit = fit_load_gen2,
##             group = "id_genero",
##             par = "loadings")
## gen_load2


## ----results='hide'---------------------------------------------------------------------------------------
gen_load_flagged2 = "c\\(lambda\\.20_5, lambda\\.20_5\\)"

mod_load_gen3 = gsub(gen_load_flagged2, "c(NA, NA)", mod_load_gen2, perl=TRUE)

cat(mod_load_gen3)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_gen3 = cfa(mod_load_gen3,
                    data = evapsi_gen,  
                    ordered = T,  
                    group = "id_genero")
summary (fit_load_gen3, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_gen3, fit_load_gen3)


## ----eval=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------------
## gen_load3 = search_ninv(data = evapsi_gen,
##             model = mod_load_gen3,
##             fit = fit_load_gen3,
##             group = "id_genero",
##             par = "loadings")
## gen_load3


## ----results='hide'---------------------------------------------------------------------------------------
gen_load_flagged3 = "c\\(lambda\\.9_3, lambda\\.9_3\\)"

mod_load_gen4 = gsub(gen_load_flagged3, "c(NA, NA)", mod_load_gen3, perl=TRUE)

cat(mod_load_gen4)


## ----warning=FALSE, results='hide'------------------------------------------------------------------------
fit_load_gen4 = cfa(mod_load_gen4,
                    data = evapsi_gen,  
                    ordered = T,  
                    group = "id_genero")
summary (fit_load_gen4, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_gen3, fit_load_gen4)


## ---------------------------------------------------------------------------------------------------------
MI_gen = matrix(NA, nrow = 3, ncol = 9)

MI_gen[1,] = round(data.matrix(fitmeasures(fit_conf_gen, fit.measures = names_fit)), 3)
MI_gen[2,] = round(data.matrix(fitmeasures(fit_thres_gen3, fit.measures = names_fit)), 3)
MI_gen[3,] = round(data.matrix(fitmeasures(fit_load_gen4, fit.measures = names_fit)), 3)


## ---------------------------------------------------------------------------------------------------------
colnames(MI_gen) = c("chisq","df","pvalue", "rmsea","rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(MI_gen) = c("configural", "partial thresholds", "partial loadings")
print(MI_gen)


## ---------------------------------------------------------------------------------------------------------
table(evapsi$study_level)


## ----results='hide'---------------------------------------------------------------------------------------
mod_conf_lvl = as.character(measEq.syntax(mod_efa,
                                data = evapsi,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "study_level",
                                group.equal = "configural"))
cat(as.character(mod_conf_lvl))


## ----results='hide'---------------------------------------------------------------------------------------
fit_conf_lvl = cfa(mod_conf_lvl,  
                           data = evapsi,            
                           ordered = T,             
                           group = "study_level")
summary (fit_conf_lvl, fit.measures = T)


## ----results='hide'---------------------------------------------------------------------------------------
mod_thres_lvl = as.character(measEq.syntax(mod_efa,
                                data = evapsi,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "study_level",
                                group.equal = "thresholds")) 
cat(as.character(mod_thres_lvl))


## ----results='hide'---------------------------------------------------------------------------------------
fit_thres_lvl = cfa(mod_thres_lvl,  
                           data = evapsi,            
                           ordered = T,             
                           group = "study_level")
summary (fit_thres_lvl, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_conf_lvl, fit_thres_lvl)


## ----results='hide'---------------------------------------------------------------------------------------
mod_load_lvl = as.character(measEq.syntax(mod_efa,
                                data = evapsi,         
                                ordered = T,           
                                parameterization = "delta",   
                                ID.fac = "std.lv",            
                                ID.cat = "Wu.Estabrook.2016", 
                                group = "study_level",
                                group.equal = c("thresholds", "loadings"))) 
cat(as.character(mod_load_lvl))


## ----results='hide'---------------------------------------------------------------------------------------
fit_load_lvl = cfa(mod_load_lvl,  
                           data = evapsi,            
                           ordered = T,             
                           group = "study_level")
summary (fit_load_lvl, fit.measures = T)


## ---------------------------------------------------------------------------------------------------------
lavTestLRT(fit_thres_lvl, fit_load_lvl)


## ---------------------------------------------------------------------------------------------------------
MI_lvl = matrix(NA, nrow = 3, ncol = 9)

MI_lvl[1,] = round(data.matrix(fitmeasures(fit_conf_lvl, fit.measures = names_fit)), 3)
MI_lvl[2,] = round(data.matrix(fitmeasures(fit_thres_lvl, fit.measures = names_fit)), 3)
MI_lvl[3,] = round(data.matrix(fitmeasures(fit_load_lvl, fit.measures = names_fit)), 3)


## ---------------------------------------------------------------------------------------------------------
colnames(MI_lvl) = c("chisq","df","pvalue", "rmsea","rmsea.ci.lower", "rmsea.ci.upper", "cfi", "tli", "srmr")
rownames(MI_lvl) = c("configural", "partial thresholds", "partial loadings")
print(MI_lvl)


## ----message=FALSE, warning=FALSE-------------------------------------------------------------------------
full_fit_efa = cfa(mod_efa, data=evapsi, ordered=TRUE, estimator="WLSM")
idx <- lavInspect(full_fit_efa, "case.idx")
EU_Fscores <- lavPredict(full_fit_efa, type= "lv", method = "EBM")
for (g in seq_along(EU_Fscores)) {
  for (fs in colnames(EU_Fscores[[g]])) {
    evapsi[ idx[[g]], fs] <- EU_Fscores[[g]][ , fs]}}

EU_Fscores = as.data.frame(EU_Fscores)
saveRDS(EU_Fscores, file="EU_scores")


## ---------------------------------------------------------------------------------------------------------
EU_Fscores = readRDS(file="EU_scores")


## ---------------------------------------------------------------------------------------------------------
EU_EU_pairs <- expand.grid(EU_factor1 = colnames(EU_Fscores),
                           EU_factor2 = colnames(EU_Fscores),
                           stringsAsFactors = FALSE)


## ---------------------------------------------------------------------------------------------------------
EU_EU_cortests <- apply(EU_EU_pairs, 1, function(x) {
  cor.test(EU_Fscores[[x["EU_factor1"]]],
           EU_Fscores[[x["EU_factor2"]]],
           method = "spearman", exact=FALSE)})


## ---------------------------------------------------------------------------------------------------------
EU_EU_results <- do.call(rbind, lapply(seq_along(EU_EU_cortests), function(i)
  {data.frame(EU_factor1 = EU_EU_pairs$EU_factor1[i],
              EU_factor2 = EU_EU_pairs$EU_factor2[i],
              rho = unname(EU_EU_cortests[[i]]$estimate),
              p_value = EU_EU_cortests[[i]]$p.value)}))
View(EU_EU_results)


## ---------------------------------------------------------------------------------------------------------
EU_DASS_pairs <- expand.grid(EU_factor   = colnames(EU_Fscores),
                             DASS_factor = colnames(evapsi[,36:39]),
                             stringsAsFactors = FALSE)


## ---------------------------------------------------------------------------------------------------------
EU_DASS_cortests <- apply(EU_DASS_pairs, 1, function(x) {
  cor.test(EU_Fscores[[x["EU_factor"]]],
           evapsi[[x["DASS_factor"]]],
           method = "spearman", exact=FALSE)})


## ---------------------------------------------------------------------------------------------------------
EU_DASS_results <- do.call(rbind, lapply(seq_along(EU_DASS_cortests), function(i)
  {data.frame(EU_factor = EU_DASS_pairs$EU_factor[i],
              DASS_factor = EU_DASS_pairs$DASS_factor[i],
              rho = unname(EU_DASS_cortests[[i]]$estimate),
              p_value = EU_DASS_cortests[[i]]$p.value)}))
EU_DASS_results


## ---------------------------------------------------------------------------------------------------------
EU_PWB_pairs <- expand.grid(EU_factor   = colnames(EU_Fscores),
                            PWB_factor = colnames(evapsi[,40:45]),
                            stringsAsFactors = FALSE)


## ---------------------------------------------------------------------------------------------------------
EU_PWB_cortests <- apply(EU_PWB_pairs, 1, function(x) {
  cor.test(EU_Fscores[[x["EU_factor"]]],
           evapsi[[x["PWB_factor"]]],
           method = "spearman", exact=FALSE)})


## ---------------------------------------------------------------------------------------------------------
EU_PWB_results <- do.call(rbind, lapply(seq_along(EU_PWB_cortests), function(i)
  {data.frame(EU_factor = EU_PWB_pairs$EU_factor[i],
              PWB_factor = EU_PWB_pairs$PWB_factor[i],
              rho = unname(EU_PWB_cortests[[i]]$estimate),
              p_value = EU_PWB_cortests[[i]]$p.value)}))
EU_PWB_results

