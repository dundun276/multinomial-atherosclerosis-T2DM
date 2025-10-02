# ==============================================================================
# External validation sample size calculation for multinomial prediction model
# ==============================================================================
# Methodology framework:
#   Primary reference: Riley et al. Stat Med. 2021;40:4230-4251. doi:10.1002/sim.9025
#   Code from: Collins GS. Logistic-regression-external-validation-sample-size GitHub repository:
#              https://github.com/gscollins1973/Logistic-regression-external-validation-sample-size
# 
# Key methodological adaption:
#   Extended Riley's binary outcome methodology to multinomial setting by 
#   decomposing into multiple binary submodels
#   Gehringer et al. J Clin Epidemiol. 2024;174:111481. doi:10.1016/j.jclinepi.2024.111481
#   GitHub: https://github.com/celinagehringer/multinomial_risk_prediction_RA
#           (File: Section3.3.3.1_external_sampsize.R)
#
# Net benefit calculation exclusion - methodological justification:
#   Riley's net benefit framework is designed for binary outcomes and cannot be 
#   directly applied to multinomial prediction models without substantial modification
# 
# Adapted by: Yongqi Zheng | Contact: dundun276@126.com
# Usage: For academic non-commercial research only. Citation required.
# ==============================================================================
library(dplyr) # version 1.1.4

### Sample size for O/E ratio
N_OE <- function(phi, SE_OE = 0.051){
  # phi = prevalence
  # SE_OE = standard error of O/E
  # Rearranged to solve for N: N = (1 - phi)/(phi * SE_OE^2)
  round((1 - phi) / (phi * SE_OE^2), 0)
  }



### Sample size for the calibration slope
N_slope <- function(LP, beta0 = 0, beta1 = 1, SE_slope = 0.051){
  # LP = linear predictor (vector)
  # beta0 = calibration intercept 
  # beta1 = calibration slope
  # SE_slope = standard error of the calibration slope
  # Based on Fisher information matrix for logistic regression
  I00 <- mean(exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  I01 <- mean(LP * exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  I11 <- mean(LP * LP * exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  round(I00 / (SE_slope^2 * (I00 * I11 - I01^2)), 0)
  }



### Sample size for the c-statistic
N_C <- function(c.index = 0.7, phi = 0.1, target.se = 0.0255, min.opt = 100, max.opt = 10000){
  # c.index = expected c-statistic
  # phi = prevalence
  # target.se = target standard error of the c-statistic
  # Based on Newcombe's non-parametric variance formula for AUC
  se.c <- function(n, cc = c.index, phi2 = phi, target.se2 = target.se){
    zz <- sqrt((cc * (1 - cc) * (1 + (n/2 - 1) * ((1 - cc)/(2 - cc)) + 
                                   ((n/2 - 1) * cc / (1 + cc)))) / (n^2 * phi2 * (1 - phi2)))
    abs(zz - target.se2)
    }
  round(optimize(se.c, c(min.opt, max.opt), tol = 0.0001)$minimum, 0)
  }



### Sample size for net benefit
# Riley's net benefit framework (N_nb function) is designed specifically for 
# binary outcomes and requires a single risk threshold for clinical decision-making
#
# For multinomial prediction models:
#   No established methodology exists for net benefit calculation
#   Multiple outcome categories complicate benefit-harm tradeoff assessments  
#   Clinical utility requires alternative approaches (e.g., decision curve analysis 
#   adapted for multiple categories, cost-effectiveness analysis)
# 
# N_nb <- function(sens, spec, phi, pt, se.nb = 0.051){
#   # sens = sensitivity
#   # spec = specificity
#   # phi = prevalence
#   # pt = threshold
#   # se.nb = standard error of net benefit
#   w    <- (1 - phi) / phi * (pt / (1 - pt))
#   round((sens * (1-sens)/phi + w^2 * spec * (1-spec)/(1-phi) + 
#            w^2 * (1-spec)^2/(phi * (1-phi))) / se.nb^2, 0)
#   }



### Minimum sample size for accurate estimate of Observed / Expected
# SE = 0.102 ¡ú 95% CI width ¡Ö 0.4 (for rare outcomes, phi < 0.3)
# SE = 0.051 ¡ú 95% CI width ¡Ö 0.2 (for common outcomes, phi > 0.3)

# Submodel 1: prevalence of thickening in TOTAL population = p.2 = 102/1123 = 0.091 ¡ú rare outcome
N_OE_1 <- N_OE(phi = p.2, SE_OE = 0.102) 
N_OE_1 # 962

# Submodel 2: prevalence of plaque in TOTAL population     = p.3 = 554/1123 = 0.493 ¡ú common outcome
N_OE_2 <- N_OE(phi = p.3, SE_OE = 0.051) 
N_OE_2 # 395



### Minimum sample size for accurate estimate of calibration slope
# Different LP distributions for different submodels reflect varying prediction patterns
# LP parameters should be obtained from development data
N <- 1000000  # Large sample for stable estimation of information matrix terms

# Submodel 1: LP ~ N(-0.339, 1.1282)
submodel1_data <- prob_and_lp_means %>% 
  filter(Outcome %in% c("0", "1")) # The codes 0 and 1 in this study represent outcomes 1 and 2
lp_stats_model1 <- submodel1_data %>%
  summarise(
    mean_LP = mean(lp_avg_cat2, na.rm = TRUE),
    sd_LP = sd(lp_avg_cat2, na.rm = TRUE),
    min_LP = min(lp_avg_cat2, na.rm = TRUE),
    max_LP = max(lp_avg_cat2, na.rm = TRUE),
    n = n(),
    n_outcome1 = sum(Outcome == "1"),
    n_outcome2 = sum(Outcome == "2"),
    prevalence_outcome2 = n_outcome2 / n()
    )
round(lp_stats_model1$mean_LP, 3) # -0.339
round(lp_stats_model1$sd_LP, 3)   # 1.128

LP_model1 <- rnorm(N, mean = -0.339, sd = 1.128)
N_slope_1 <- N_slope(LP_model1, SE_slope = 0.072)
N_slope_1 # 1187

# Submodel 2: LP ~ N(0.067, 1.4952)
submodel2_data <- prob_and_lp_means %>% 
  filter(Outcome %in% c("0", "2")) # The codes 0 and 2 in this study represent outcomes 1 and 3
lp_stats_model2 <- submodel2_data %>%
  summarise(
    mean_LP = mean(lp_avg_cat3, na.rm = TRUE),
    sd_LP = sd(lp_avg_cat3, na.rm = TRUE),
    min_LP = min(lp_avg_cat3, na.rm = TRUE),
    max_LP = max(lp_avg_cat3, na.rm = TRUE),
    n = n(),
    n_outcome1 = sum(Outcome == "1"),
    n_outcome3 = sum(Outcome == "3"),
    prevalence_outcome3 = n_outcome3 / n()
    )
round(lp_stats_model2$mean_LP, 3) # -0.067
round(lp_stats_model2$sd_LP, 3)   # 1.495
LP_model2 <- rnorm(N, mean = 0.067, sd = 1.495)
N_slope_2 <- N_slope(LP_model2, SE_slope = 0.072)
N_slope_2 # 923



### Minimum sample size for accurate estimate of concordance statistic
# Riley: "We suggest assuming C-statistics that are the same or about 0.05 lower than 
#         the C-statistic reported for model development"
# SE = 0.0300 ¡ú 95% CI width ¡Ö 0.12 for rare outcomes
# SE = 0.0255 ¡ú 95% CI width ¡Ö 0.10 for common outcomes

# Submodel 1: rare outcome, anticipated C-statistic = 0.74 (internal validation corrected value)
N_C_1 <- N_C(c.index = 0.74, phi = p.2, target.se = 0.0300)
N_C_1 # 819

# Submodel 2: common outcome, anticipated C-statistic = 0.79 (internal validation corrected value)
N_C_2 <- N_C(c.index = 0.79, phi = p.3, target.se = 0.0255)
N_C_2 # 315



### Summary
external_sample_size <- max(N_OE_1, N_slope_1, N_C_1, N_OE_2, N_slope_2, N_C_2)
external_sample_size # 1187

### Minimum sample size: n = 1187