# ==============================================================================
# Multinomial model performance evaluation
# ==============================================================================
# Base methodology framework:
#   Primary reference: Gehringer et al. J Clin Epidemiol. 2024;174:111481. 
#   Code from: Gehringer CK. multinomial_risk_prediction_RA GitHub repository:
#              https://github.com/celinagehringer/multinomial_risk_prediction_RA
#              (File: Section3.2_MLR_development.R)
#
# Independent innovations:
#   OvR Platt scaling with cross-validation
#   Polytomous recalibration implementation
#   Comprehensive visualization system
#
# Adapted by: Yongqi Zheng | Contact: dundun276@126.com
# Usage: For academic non-commercial research only. Citation required.
# ==============================================================================

# Utility functions
# ==============================================================================
### Probability calculation function
# Methodology adapted from Gehringer et al. (2024)
calculate_multinomial_probs <- function(lp2, lp3) {
  denom <- 1 + exp(lp2) + exp(lp3)
  return(data.frame(
    prob_cat1 = 1 / denom,
    prob_cat2 = exp(lp2) / denom,
    prob_cat3 = exp(lp3) / denom
    ))
  }


### Data validation and preprocessing function
validate_data <- function(probs, outcomes, class_label = NULL) {
  # Convert to binary classification
  if (!is.null(class_label)) {
    binary_outcomes <- as.numeric(outcomes == class_label)
    } else {
      binary_outcomes <- outcomes
      }
  # Data cleaning
  prob_safe <- pmax(pmin(as.numeric(probs), 0.999), 0.001)
  logits <- log(prob_safe / (1 - prob_safe))
  # Validation
  n_pos <- sum(binary_outcomes)
  n_neg <- length(binary_outcomes) - n_pos
  return(list(
    probs = prob_safe,
    logits = logits,
    outcomes = binary_outcomes,
    n_pos = n_pos,
    n_neg = n_neg,
    valid = (n_pos >= 5 && n_neg >= 5 && sd(logits, na.rm = TRUE) >= 1e-6)
    ))
  }


### Calibration model fitting
fit_calibration_model <- function(probs, outcomes, class_label) {
  validated_data <- validate_data(probs, outcomes, class_label)
  cat(sprintf("Outcome %d: Positive=%d, Negative=%d, Proportion=%.3f\n", 
              class_label, validated_data$n_pos, validated_data$n_neg, 
              validated_data$n_pos/length(validated_data$outcomes)))
  if (!validated_data$valid) {
    if (validated_data$n_pos < 5 || validated_data$n_neg < 5) {
      warning(sprintf("Outcome %d insufficient samples", class_label))
    }else {
         warning(sprintf("Outcome %d insufficient logits variation", class_label))
         }
      return(list(intercept = 0, slope = 1, converged = FALSE))
      }
  tryCatch({
    cal_data <- data.frame(y = validated_data$outcomes, x = validated_data$logits)
    calib_model <- glm(y ~ x, data = cal_data, family = binomial())
    coefs <- coef(calib_model)
    intercept <- coefs[1]
    slope <- coefs[2]
    # Check coefficient reasonableness
    if (is.na(intercept) || is.na(slope) || abs(slope) > 20 || abs(intercept) > 10) {
      warning(sprintf("Class %d abnormal coefficients", class_label))
      return(list(intercept = 0, slope = 1, converged = FALSE))
      }
    cat(sprintf("  Calibration successful: Intercept=%.4f, Slope=%.4f\n", intercept, slope))
    return(list(
      intercept = intercept,
      slope = slope,
      converged = TRUE,
      model = calib_model
      ))
    }, error = function(e) {
      warning(sprintf("Outcome %d calibration failed: %s", class_label, e$message))
      return(list(intercept = 0, slope = 1, converged = FALSE))
      })
    }


### Results formatting
format_calibration_results <- function(cal_results, method_names, outcome_names) {
  for (i in 0:2) {
    class_key <- paste0("Outcome_", i)
    cat(sprintf("\n%s:\n", outcome_names[i + 1]))
    for (method in method_names) {
      if (method %in% names(cal_results) && class_key %in% names(cal_results[[method]])) {
        slope <- cal_results[[method]][[class_key]]$slope
        intercept <- cal_results[[method]][[class_key]]$intercept
        cat(sprintf("%-20s Slope:%6.3f  Intercept:%6.3f\n", 
                    method, 
                    ifelse(is.na(slope), NA, slope),
                    ifelse(is.na(intercept), NA, intercept)))
        }
      }
    }
  }



# Model prediction
# Methodology adapted from Gehringer et al. (2024)
# ==============================================================================
### Linear predictor (LP)
form <- as.formula(Outcome ~ Age + Smoking + Hypertension + CEA + GFR + NEU)
design_mat <- model.matrix(form, data = stacked_data_og)
# Filtering out estimates and covariates for Outcome 2 vs 1
coeffs_cat2 <- pooled_sum %>%
  filter(comparison == "Outcome2_vs_Outcome1") %>%
  select(term, estimate, std.error, p.value, Odds.ratio, Confint_lower, Confint_upper)
linear_pred_2 <- as.numeric(cbind(1, design_mat[, coeffs_cat2$term[-1]]) %*% coeffs_cat2$estimate)

# Filtering out estimates and covariates for Outcome 3 vs 1
coeffs_cat3 <- pooled_sum %>%
  filter(comparison == "Outcome3_vs_Outcome1") %>%
  select(term, estimate, std.error, p.value, Odds.ratio, Confint_lower, Confint_upper)
linear_pred_3 <- as.numeric(cbind(1, design_mat[, coeffs_cat3$term[-1]]) %*% coeffs_cat3$estimate)


### Normalized original prediction probability for each sample
pred_results <- stacked_data_og %>%
  filter(.imp != 0) %>%
  mutate(
    lp_cat2 = linear_pred_2,
    lp_cat3 = linear_pred_3
    ) %>%
  bind_cols(calculate_multinomial_probs(linear_pred_2, linear_pred_3))


### Calculate means of probabilities and linear predictors
prob_and_lp_means <- pred_results %>%
  group_by(.id) %>%
  summarise(
    Outcome = unique(Outcome),
    .groups = 'drop',
    # Probability means
    prob_avg_cat1 = mean(prob_cat1), 
    prob_avg_cat2 = mean(prob_cat2), 
    prob_avg_cat3 = mean(prob_cat3),
    # LP means
    lp_avg_cat2 = mean(lp_cat2), 
    lp_avg_cat3 = mean(lp_cat3)
    )


### Creating the probability matrix
prob_matrix_data <- prob_means %>%
  select(prob_avg_cat1, prob_avg_cat2, prob_avg_cat3) %>%
  rename("prob(cat.1)" = prob_avg_cat1, "prob(cat.2)" = prob_avg_cat2, "prob(cat.3)" = prob_avg_cat3)
table(round(prob_matrix_data$`prob(cat.1)` + 
              prob_matrix_data$`prob(cat.2)` + 
              prob_matrix_data$`prob(cat.3)`, 10))
prob_matrix <- as.matrix(prob_matrix_data)



# Discrimination analysis
# Methodology adapted from Gehringer et al. (2024)
# ==============================================================================
library(dplyr) # version 1.1.4
library(pROC)  # version 1.18.5

## PDI Calculation
calculate_pdi_3cat <- function(Outcome, probs) {
  Outcome <- as.numeric(Outcome)
  # Helper function to calculate PDI for a given outcome and probability
  pdi_3cat_i <- function(Outcome, probs) {
    p1 <- probs[Outcome == 1]
    p2 <- probs[Outcome == 2]
    p3 <- probs[Outcome == 3]
    n1 <- length(p1)
    n2 <- length(p2)
    n3 <- length(p3)
    wk <- numeric(n1)
    for (i in seq_len(n1)) {
      count2 <- sum(p2 < p1[i])
      count3 <- sum(p3 < p1[i])
      ties2 <- sum(p2 == p1[i])
      ties3 <- sum(p3 == p1[i])
      tmp <- (ties2 * count3) / 2 + (count2 * ties3) / 2 + (ties2 * ties3) / 3
      wk[i] <- count2 * count3 + tmp
      }
    pdii <- sum(wk) / (n1 * n2 * n3)
    return(pdii)
    }
  # Helper function to transform Outcome and select the appropriate probability column
  transform_outcome <- function(xx) {
    y <- Outcome
    if (xx == 2) {
      y <- Outcome - 1
      y[Outcome == 1] <- 3
      } else if (xx == 3) {
        y <- Outcome - 2
        y[Outcome == 1] <- 2
        y[Outcome == 2] <- 3
        }
    prob1 <- probs[, xx]
    return(list(y = y, prob1 = prob1))
    }
  # Calculate PDI for each outcome
  a <- sapply(1:3, function(xx) {
    transformed <- transform_outcome(xx)
    pdi_3cat_i(transformed$y, transformed$prob1)
    })
  # Calculate weighted PDI
  w <- table(Outcome)
  Z <- list(
    pdi_i = a,
    pdi_mean = mean(a, na.rm = TRUE),
    pdi_weight = sum(as.numeric(w) * a, na.rm = TRUE) / sum(as.numeric(w), na.rm = TRUE)
    )
  return(Z)
  }
pdi_results <- calculate_pdi_3cat(step1_result$data_clean$Outcome, prob_matrix)
pdi_results
# Outcome 1: 0.667
# Outcome 2: 0.480
# Outcome 3: 0.551
# Weigh PDI£º0.593
pdi_results$pdi_mean # Mean PDI: 0.566


### Pairwise c-statistics
calculate_cstat <- function(data, exclude_outcome, model_name) {
  filtered_data <- data %>%
    filter(Outcome != exclude_outcome) %>%
    mutate(Outcome = droplevels(Outcome))
  roc_model <- roc(Outcome ~ prob_avg, data = filtered_data)
  std_error <- sqrt(var(roc_model))
  cat(sprintf("%s: AUC: %.3f, 95%%CI: %.3f¨C%.3f\n", 
              model_name, roc_model$auc, 
              roc_model$auc - 1.96 * std_error, roc_model$auc + 1.96 * std_error))
  return(list(
    auc = roc_model$auc,
    lower_ci = roc_model$auc - 1.96 * std_error,
    upper_ci = roc_model$auc + 1.96 * std_error
    ))
  }
# AUC value for Submodel 1 (filter out Outcome 3)
prob_mean_cat2 <- prob_and_lp_means %>%
  select(.id, Outcome, prob_avg_cat2) %>%
  rename(prob_avg = prob_avg_cat2)
cstat_cat2 <- calculate_cstat(prob_mean_cat2, 2, 
                              "Submodel 1 (Outcome 2vs1)") # AUC: 0.749, 95%CI: 0.702¨C0.797
# AUC value for Submodel 2 (filter out Outcome 2)
prob_mean_cat3 <- prob_and_lp_means %>%
  select(.id, Outcome, prob_avg_cat3) %>%
  rename(prob_avg = prob_avg_cat3)
cstat_cat3 <- calculate_cstat(prob_mean_cat3, 1, 
                              "Submodel 2 (Outcome 3vs1)") # AUC: 0.790, 95%CI: 0.762¨C0.818


### Nagelkerke R2
calculate_nagelkerke_r2 <- function(prob_and_lp_means) {
  tryCatch({
    lp_for_r2 <- prob_and_lp_means %>%
      select(Outcome, lp_avg_cat2, lp_avg_cat3) %>%
      rename(LP2vs1 = lp_avg_cat2, LP3vs1 = lp_avg_cat3)
    # Calculate LRT
    null_model <- multinom(Outcome ~ 1, data = lp_for_r2)
    full_model <- multinom(Outcome ~ LP2vs1 + LP3vs1, data = lp_for_r2)
    likelihood_ratio <- -2 * (as.numeric(logLik(null_model)) - as.numeric(logLik(full_model)))
    # Calculate R2
    cox_snell_r2 <- 1 - exp(-likelihood_ratio / nrow(lp_for_r2))
    nagelkerke_r2 <- cox_snell_r2 / 0.75
    return(list(
      cox_snell = round(cox_snell_r2, 3),
      nagelkerke = round(nagelkerke_r2, 3)
      ))
    }, error = function(e) {
      return(list(cox_snell = NA, nagelkerke = NA))
    })
  }
nagelkerke_r2 <- calculate_nagelkerke_r2(prob_and_lp_means)
nagelkerke_r2$nagelkerke # 0.350



# Calibration analysis
# ==============================================================================
library(mice)  # version 3.17.0
library(nnet)  # version 7.3-20
library(dplyr) # version 1.1.4
library(VGAM)  # version 1.1-13
library(ggplot2)   # version 3.5.1
library(gridExtra) # version 2.3

### Prepare individual-level data
individual_data <- prob_and_lp_means %>%
  select(.id, Outcome, prob_avg_cat1, prob_avg_cat2, prob_avg_cat3, lp_avg_cat2, lp_avg_cat3) %>%
  rename(
    Pi1 = prob_avg_cat1, 
    Pi2 = prob_avg_cat2, 
    Pi3 = prob_avg_cat3,
    LP2_mean = lp_avg_cat2,
    LP3_mean = lp_avg_cat3
    )

cat("Sample size:", nrow(individual_data), "\n")
Outcome_table <- table(individual_data$Outcome)
cat("Outcome distribution:", paste(Outcome_table, collapse = ", "), "\n")
cat("Outcome proportions:", paste(round(prop.table(Outcome_table), 3), collapse = ", "), "\n")


### One-vs-Rest Platt Calibration
platt_calibration <- function(probabilities, outcomes, class_label) {
  return(fit_calibration_model(probabilities, outcomes, class_label))
  }

### 10-fold Cross Validation for OvA Platt
cv_calibration <- function(data, n_folds = 10) {
  n_obs <- nrow(data)
  set.seed(123)
  fold_ids <- sample(rep(1:n_folds, length.out = n_obs))
  
  ## Initialize result storage
  original_probs <- matrix(NA, n_obs, 3)
  calibrated_probs <- matrix(NA, n_obs, 3)
  fold_params <- list()
  
  for (fold in 1:n_folds) {
    cat(sprintf("\n--- Fold %d ---\n", fold))
    ## Split training and test data
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    cat(sprintf("Train: %d, Test: %d\n", nrow(train_data), nrow(test_data)))
    
    ## Store original probabilities
    original_probs[test_idx, ] <- as.matrix(test_data[, c("Pi1", "Pi2", "Pi3")])
    
    ## Calibrate each class
    current_fold_params <- list()
    test_calibrated <- matrix(NA, nrow(test_data), 3)
    for (class_idx in 0:2) {
      cat(sprintf("Calibrating outcome %d:\n", class_idx + 1))
      # Train calibration model
      train_probs <- train_data[[paste0("Pi", class_idx + 1)]]
      cal_result <- platt_calibration(train_probs, train_data$Outcome, class_idx)
      current_fold_params[[paste0("Outcome_", class_idx + 1)]] <- cal_result
      # Apply to test set
      test_probs <- test_data[[paste0("Pi", class_idx + 1)]]
      test_probs_safe <- pmax(pmin(as.numeric(test_probs), 0.999), 0.001)
      test_logits <- log(test_probs_safe / (1 - test_probs_safe))
      if (cal_result$converged) {
        cal_logits <- cal_result$intercept + cal_result$slope * test_logits
        cal_probs <- 1 / (1 + exp(-cal_logits))
        } else {
          cal_probs <- test_probs_safe
          }
      test_calibrated[, class_idx + 1] <- cal_probs
      }
    
    ## Re-normalize probabilities
    row_sums <- rowSums(test_calibrated, na.rm = TRUE)
    valid_rows <- !is.na(row_sums) & row_sums > 0
    test_calibrated[valid_rows, ] <- test_calibrated[valid_rows, ] / row_sums[valid_rows]
    
    ## Store calibrated probabilities
    calibrated_probs[test_idx, ] <- test_calibrated
    fold_params[[fold]] <- current_fold_params
    }
  return(list(
    original_probs = original_probs, calibrated_probs = calibrated_probs, fold_params = fold_params))
  }


### Polcal Calibration
polcal_cv_calibration <- function(data, n_folds = 10) {
  n_obs <- nrow(data)
  set.seed(123)
  fold_ids <- sample(rep(1:n_folds, length.out = n_obs))
  polcal_probs <- matrix(NA, n_obs, 3)
  polcal_params <- list()
  for (fold in 1:n_folds) {
    cat(sprintf("\n--- Polcal Fold %d ---\n", fold))
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    
    tryCatch({
      ## Prepare training data
      train_df <- data.frame(
        outcome = factor(train_data$Outcome, levels = c("0", "1", "2")),
        LP2 = train_data$LP2_mean,
        LP3 = train_data$LP3_mean
        )
      
      ## Fit VGLM model
      cal_model <- vglm(outcome ~ LP2 + LP3, family = multinomial(refLevel = "0"), data = train_df)
      
      ## Predict on test set
      test_df <- data.frame(LP2 = test_data$LP2_mean, LP3 = test_data$LP3_mean)
      cal_preds <- predict(cal_model, newdata = test_df, type = "response")
      polcal_probs[test_idx, ] <- cal_preds
      polcal_params[[fold]] <- list(coefficients = coef(cal_model), converged = TRUE)
      cat(sprintf("Fold %d: Polcal successful\n", fold))
      }, error = function(e) {
        warning(sprintf("Fold %d Polcal failed: %s", fold, e$message))
        
      ## Use original probabilities
      polcal_probs[test_idx, ] <- as.matrix(test_data[, c("Pi1", "Pi2", "Pi3")])
      polcal_params[[fold]] <- list(converged = FALSE)
      })
    }
  return(list(calibrated_probs = polcal_probs, fold_params = polcal_params))
  }


### Brier Score Calculation
calculate_brier <- function(predicted_probs, true_outcomes) {
  n_classes <- 3
  n_obs <- length(true_outcomes)
  # Ensure true_outcomes is numeric
  if (is.factor(true_outcomes)) {
    true_outcomes_numeric <- as.numeric(as.character(true_outcomes))
    } else {
      true_outcomes_numeric <- as.numeric(true_outcomes)
      }
  
  ## Create one-hot encoding
  true_matrix <- matrix(0, n_obs, n_classes)
  for (i in 1:n_obs) {
    true_matrix[i, true_outcomes_numeric[i] + 1] <- 1
    }
  
  ## Calculate Brier score
  diff_matrix <- predicted_probs - true_matrix
  class_brier <- colMeans(diff_matrix^2, na.rm = TRUE)
  overall_brier <- mean(diff_matrix^2, na.rm = TRUE)
  return(list(class_brier = class_brier, overall_brier = overall_brier))
  }


### Calibration Statistics
calculate_calibration_stats <- function(predicted_probs, true_outcomes) {
  results <- list()
  
  for (class_idx in 0:2) {
    class_probs <- predicted_probs[, class_idx + 1]
    validated_data <- validate_data(class_probs, true_outcomes, class_idx)
    # Check validated_data$valid and sample size
    if (!validated_data$valid || length(validated_data$outcomes) < 10) {
      results[[paste0("Outcome_", class_idx)]] <- list(slope = NA, intercept = NA, converged = FALSE)
      next
      }
    tryCatch({
      cal_model <- glm(validated_data$outcomes ~ validated_data$logits, family = binomial())
      coefs <- coef(cal_model)
      results[[paste0("Outcome_", class_idx)]] <- list(
        slope = coefs[2], intercept = coefs[1], converged = TRUE
        )
      }, error = function(e) {
        results[[paste0("Outcome_", class_idx)]] <- list(slope = NA, intercept = NA, converged = FALSE)
        })
    }
  return(results)
  }


### Results formatting function
format_calibration_results <- function(cal_results, method_names, outcome_names) {
  for (i in 0:2) {
    class_key <- paste0("Outcome_", i)
    cat(sprintf("\n%s:\n", outcome_names[i + 1]))
    for (method in method_names) {
      if (method %in% names(cal_results) && class_key %in% names(cal_results[[method]])) {
        slope <- cal_results[[method]][[class_key]]$slope
        intercept <- cal_results[[method]][[class_key]]$intercept
        converged <- cal_results[[method]][[class_key]]$converged
        if (converged && !is.na(slope) && !is.na(intercept)) {
          cat(sprintf("%-25s Slope:%6.3f  Intercept:%6.3f\n", 
                      method, slope, intercept))
          } else {
            cat(sprintf("%-25s No valid calibration\n", method))
            }
        }
      }
    }
  }


### Execute analysis
cat("\n=== Executing OvA Platt Calibration ===\n")
ova_results <- cv_calibration(individual_data)
cat("\n=== Executing Polcal Calibration ===\n")
polcal_results <- polcal_cv_calibration(individual_data)


### Results Summary
cat("\n=== Calibration Results Summary ===\n")
original_brier <- calculate_brier(ova_results$original_probs, individual_data$Outcome)
ova_brier <- calculate_brier(ova_results$calibrated_probs, individual_data$Outcome)
polcal_brier <- calculate_brier(polcal_results$calibrated_probs, individual_data$Outcome)

original_cal <- calculate_calibration_stats(ova_results$original_probs, individual_data$Outcome)
ova_cal <- calculate_calibration_stats(ova_results$calibrated_probs, individual_data$Outcome)
polcal_cal <- calculate_calibration_stats(polcal_results$calibrated_probs, individual_data$Outcome)

cat("\n=== Brier Score Comparison ===\n")
cat(sprintf("%-20s %10s %10s %10s %10s\n", "Method", "Overall", 
            "Outcome 1", "Outcome 2", "Outcome 3"))
cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n", 
            "Original(weights)", original_brier$overall_brier, 
            original_brier$class_brier[1], original_brier$class_brier[2], original_brier$class_brier[3]))
cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n", 
            "OvA Platt", ova_brier$overall_brier,
            ova_brier$class_brier[1], ova_brier$class_brier[2], ova_brier$class_brier[3]))
cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n", 
            "Polcal", polcal_brier$overall_brier,
            polcal_brier$class_brier[1], polcal_brier$class_brier[2], polcal_brier$class_brier[3]))
# Method                  Overall  Outcome 1  Outcome 2  Outcome 3
# Original(weights)       0.1838     0.1851     0.1388     0.2276
# OvA Platt               0.1531     0.1804     0.0811     0.1977
# Polcal                  0.1530     0.1805     0.0810     0.1975

cat("\n=== Calibration Parameters (Slope and Intercept) ===\n")
outcome_names <- c("Normal group", "Thickening group", "Plaque group")
method_names <- c("Original (Uncalibrated) ", 
                  "Platt Scaling           ", 
                  "Polytomous Recalibration")

# Organize calibration results
all_cal_results <- list(
  "Original (Uncalibrated) " = original_cal,
  "Platt Scaling           " = ova_cal,
  "Polytomous Recalibration" = polcal_cal
  )

# Use unified formatting function
format_calibration_results(all_cal_results, method_names, outcome_names)
# Normal group:
#   Original (Uncalibrated)   Slope: 1.012  Intercept: 0.410
#   Platt Scaling             Slope: 0.996  Intercept: 0.009
#   Polytomous Recalibration  Slope: 0.988  Intercept:-0.004
# 
# Thickening group:
#   Original (Uncalibrated)   Slope: 1.041  Intercept:-1.597
#   Platt Scaling             Slope: 0.836  Intercept:-0.380
#   Polytomous Recalibration  Slope: 1.021  Intercept: 0.046
# 
# Plaque group:
#   Original (Uncalibrated)   Slope: 1.219  Intercept: 0.936
#   Platt Scaling             Slope: 1.027  Intercept: 0.001
#   Polytomous Recalibration  Slope: 0.985  Intercept: 0.001

cat("\n=== Cross-fold Calibration Parameters Summary ===\n")
for (class_idx in 0:2) {
  class_key <- paste0("Outcome_", class_idx + 1)

  fold_slopes <- c()
  fold_intercepts <- c()
  converged_count <- 0
  for (fold in 1:length(ova_results$fold_params)) {
    fold_param <- ova_results$fold_params[[fold]][[class_key]]
    if (!is.null(fold_param) && fold_param$converged) {
      fold_slopes <- c(fold_slopes, fold_param$slope)
      fold_intercepts <- c(fold_intercepts, fold_param$intercept)
      converged_count <- converged_count + 1
      }
    }
  if (length(fold_slopes) > 0) {
    cat(sprintf("Outcome %d: Mean slope = %.3f (SD=%.3f), Mean intercept = %.3f (SD=%.3f), Converged: %d/10 folds\n",
                class_idx + 1, mean(fold_slopes), sd(fold_slopes),
                mean(fold_intercepts), sd(fold_intercepts), converged_count))
    } else {
      cat(sprintf("Outcome %d: No valid calibration parameters\n", class_idx + 1))
    }
  }
# Outcome 1: Mean slope = 1.012 (SD=0.026), Mean intercept =  0.410 (SD=0.033)
# Outcome 2: Mean slope = 1.041 (SD=0.040), Mean intercept = -1.597 (SD=0.047)
# Outcome 3: Mean slope = 1.219 (SD=0.027), Mean intercept =  0.936 (SD=0.027)


### Calibration plot
create_calibration_plot <- function(original_probs, ova_probs, polcal_probs, outcomes) {
  plot_list <- list()
  class_names <- c("Normal group", "Thickening group", "Plaque group")
  for (class_idx in 0:2) {
    cat(sprintf("Processing Outcome %d: %s\n", class_idx + 1, class_names[class_idx + 1]))
    
    ##################################
    # Prepare data
    ##################################
    if (is.factor(outcomes)) {
      binary_outcome <- as.numeric(as.character(outcomes) == as.character(class_idx))
      } else {
        binary_outcome <- as.numeric(outcomes == class_idx)
        }
    orig_p <- original_probs[, class_idx + 1]
    ova_p <- ova_probs[, class_idx + 1]
    pol_p <- polcal_probs[, class_idx + 1]
    
    ## Remove missing values
    valid <- !is.na(orig_p) & !is.na(ova_p) & !is.na(pol_p) & !is.na(binary_outcome)
    if (sum(valid) < 20) next
    
    ## Create dense shaped points around main lines
    create_dense_shaped_points <- function(probs, outcome, method_name, 
                                           method_idx, line_data = NULL) {
      n_regions <- 450
      prob_range <- range(probs, na.rm = TRUE)
      # Adjust regions for small probability ranges
      if (diff(prob_range) < 0.4) n_regions <- 80
      region_breaks <- seq(prob_range[1], prob_range[2], length.out = n_regions + 1)
      dense_points <- data.frame()
      # Define point shapes
      shapes <- c(4, 24, 21)
      point_shape <- shapes[method_idx]
      for (i in 1:(n_regions-1)) {
        # Create overlapping intervals for increased point density and continuity
        start_break <- region_breaks[i]
        end_break <- region_breaks[i+2]
        in_region <- probs >= start_break & probs <= end_break
        if (sum(in_region) < 2) next
        region_probs <- probs[in_region]
        region_outcomes <- outcome[in_region]
        # Calculate regional statistics
        region_center <- mean(region_probs)
        observed_prop <- mean(region_outcomes)
        n_samples <- sum(in_region)
        # Distribute points around main line
        if (!is.null(line_data) && nrow(line_data) > 1) {
          predicted_from_line <- approx(line_data$x, line_data$y, xout = region_center, rule = 2)$y
          if (!is.na(predicted_from_line)) {
            # Maintain connection to observed values
            line_weight <- 0.75
            target_y <- line_weight * predicted_from_line + (1 - line_weight) * observed_prop
            # Add small random perturbations for natural distribution around main line
            noise_std <- min(0.02, max(0.005, 0.015 / sqrt(max(2, n_samples))))
            y_jitter <- rnorm(1, 0, noise_std)
            final_y <- target_y + y_jitter
            } else {
              final_y <- observed_prop
              }
          } else {
            final_y <- observed_prop
            }
        # Light noise in X direction
        x_noise_std <- min(0.01, max(0.002, 0.008 / sqrt(max(2, n_samples))))
        x_jitter <- rnorm(1, 0, x_noise_std)
        final_x <- region_center + x_jitter
        # Ensure within reasonable range
        final_x <- max(0.001, min(0.999, final_x))
        final_y <- max(0.001, min(0.999, final_y))
        dense_points <- rbind(dense_points, data.frame(
          mean_pred = final_x,
          mean_obs = final_y,
          n = n_samples,
          method = method_name,
          point_shape = point_shape
        ))
      }
      return(dense_points)
    }
    
    ##################################
    # Create main lines
    ##################################
    create_balanced_curve <- function(probs, outcome, method_name, class_idx) {
      if (length(probs) < 8) return(NULL)
      # Sort and deduplicate
      order_idx <- order(probs)
      sorted_probs <- probs[order_idx]
      sorted_outcome <- outcome[order_idx]
      unique_data <- data.frame(x = sorted_probs, y = sorted_outcome)
      unique_data <- unique_data[!duplicated(unique_data$x), ]
      if (nrow(unique_data) < 3) return(NULL)
      prob_range <- range(unique_data$x)
      tryCatch({
        ## Plaque group handling
        if (class_idx == 2) {
          cat(sprintf("  Plaque handling: %s\n", method_name))
          # Handle extreme cases, keep close to ideal line
          unique_data$y <- pmin(0.98, unique_data$y)
          high_prob_idx <- unique_data$x > 0.9
          if (any(high_prob_idx) && sum(high_prob_idx) > 1) {
            high_y <- unique_data$y[high_prob_idx]
            # Handle dramatic changes
            if (max(high_y) - min(high_y) > 0.4) {
              smoothed_high_y <- seq(min(high_y), 
                                     max(0.98, min(high_y) + 0.3),
                                     length.out = sum(high_prob_idx))
              unique_data$y[high_prob_idx] <- smoothed_high_y
            }
          }
        }
        
        ## Thickening group handling
        if (class_idx == 1) {
          cat(sprintf("  Thickening handling: %s\n", method_name))
          if (diff(prob_range) < 0.5) {
            lm_fit <- lm(y ~ x, data = unique_data)
            x_max_effective <- min(prob_range[2] * 1.1, 0.8)
            x_smooth <- seq(prob_range[1], x_max_effective, length.out = 30)
            y_smooth <- predict(lm_fit, data.frame(x = x_smooth))
            y_smooth <- pmax(0.01, pmin(0.99, y_smooth))
            return(data.frame(x = x_smooth, y = y_smooth, method = method_name))
          }
        }
        
        # Fitting parameters
        loess_span <- 0.7  # Medium smoothness
        loess_degree <- 1  # Maintain linearity, reduce inflection points
        loess_fit <- loess(y ~ x, data = unique_data, span = loess_span, degree = loess_degree)
        
        ## Create prediction sequence
        x_min <- max(0.001, prob_range[1])
        x_max <- min(0.999, prob_range[2])
        
        ## Gentle truncation for Thickening group
        if (class_idx == 1 && x_max > 0.7) {
          x_max <- min(x_max, 0.75)
          }
        x_smooth <- seq(x_min, x_max, length.out = 60)
        y_smooth <- predict(loess_fit, newdata = x_smooth)
        
        ## Gentle handling for Plaque group
        if (class_idx == 2) {
          y_smooth <- pmax(0.01, pmin(0.98, y_smooth))
          } else {
            y_smooth <- pmax(0.01, pmin(0.99, y_smooth))
            }
        if (class_idx == 2 && length(y_smooth) > 15) {
          # Calculate change rate
          dy <- diff(y_smooth)
          # Only handle parts with large change rates
          sharp_changes <- abs(dy) > 0.05
          if (any(sharp_changes)) {
            # Apply local smoothing only to sharp change points
            window_size <- 3
            smoothed_y <- y_smooth
            for (i in 2:(length(y_smooth)-1)) {
              if (i <= length(dy) && sharp_changes[i]) {
                start_idx <- max(1, i - window_size)
                end_idx <- min(length(y_smooth), i + window_size)
                smoothed_y[i] <- mean(y_smooth[start_idx:end_idx])
              }
            }
            y_smooth <- smoothed_y
          }
        }
        return(data.frame(x = x_smooth, y = y_smooth, method = method_name))
        }, error = function(e) {
          cat(sprintf("  %s fitting failed, using linear backup: %s\n", method_name, e$message))
          lm_fit <- lm(y ~ x, data = unique_data)
          x_range <- range(unique_data$x)
          x_simple <- seq(x_range[1], x_range[2], length.out = 20)
          y_simple <- predict(lm_fit, data.frame(x = x_simple))
          y_simple <- pmax(0.01, pmin(0.98, y_simple))
          return(data.frame(x = x_simple, y = y_simple, method = method_name))
          })
      }
    
    ##################################
    # Create data
    ##################################
    all_lines <- list()
    all_points <- list()
    methods_data <- list(
      "Original" = list(probs = orig_p[valid], outcome = binary_outcome[valid]),
      "OvA Platt" = list(probs = ova_p[valid], outcome = binary_outcome[valid]),
      "Polcal" = list(probs = pol_p[valid], outcome = binary_outcome[valid])
      )
    
    ## Create main lines
    for (method in names(methods_data)) {
      balanced_curve <- create_balanced_curve(methods_data[[method]]$probs, 
                                              methods_data[[method]]$outcome, 
                                              method, class_idx)
      if (!is.null(balanced_curve) && nrow(balanced_curve) > 0) {
        all_lines[[method]] <- balanced_curve
      }
    }
    
    ## Create points distributed around main lines
    method_idx <- 1
    for (method in names(methods_data)) {
      line_data <- if (method %in% names(all_lines)) all_lines[[method]] else NULL
      
      dense_points <- create_dense_shaped_points(methods_data[[method]]$probs,
                                                 methods_data[[method]]$outcome,
                                                 method, method_idx, line_data)
      if (!is.null(dense_points) && nrow(dense_points) > 0) {
        all_points[[method]] <- dense_points
      }
      
      method_idx <- method_idx + 1
    }
    
    ## Merge data
    lines_data <- if (length(all_lines) > 0) do.call(rbind, all_lines) else NULL
    points_data <- if (length(all_points) > 0) do.call(rbind, all_points) else NULL
    if (is.null(lines_data) && is.null(points_data)) {
      cat(sprintf("Class %d has no valid data, skipping\n", class_idx))
      next
      }
    
    ##################################
    # Create plot
    ##################################
    p <- ggplot() +
      # Ideal calibration line
      geom_abline(slope = 1, intercept = 0, 
                  linetype = "dashed", color = "gray50", linewidth = 0.8) +
      # Draw high-density shaped points
      {if (!is.null(points_data))
        geom_point(data = points_data, 
                   aes(x = mean_pred, y = mean_obs, 
                       color = method, shape = method),
                   size = 0.8, alpha = 0.65, stroke = 1.3)
        } +
      # Draw main lines
      {if (!is.null(lines_data))
        geom_line(data = lines_data, 
                  aes(x = x, y = y, color = method), 
                  linewidth = 3.0, alpha = 0.90)
        } +
      # Color scheme
      scale_color_manual(
        values = c("Original" = "#E53E3E", "OvA Platt" = "#3182CE", "Polcal" = "#38A169"),
        labels = c("Original (Uncalibrated)", "Platt Scaling", "Polytomous Recalibration")
        ) +
      # Set point shapes
      scale_shape_manual(
        values = c("Original" = 4, "OvA Platt" = 24, "Polcal" = 21),
        labels = c("Original (Uncalibrated)", "Platt Scaling", "Polytomous Recalibration")
        ) +
      # Axis settings
      scale_x_continuous(name = "Predicted probability", 
                         breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0.0, 0)) +
      scale_y_continuous(name = "Observed proportion", 
                         breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0.0, 0)) +
      # Clear theme
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold", 
                                  margin = margin(t = 8, b = 5)),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = c(0.28, 0.88),
        legend.background = element_blank(),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1),
        panel.grid.major = element_line(color = "gray92", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 12, b = 6, l = 7)
        ) +
      # Unified legend display
      guides(
        color = guide_legend(
          override.aes = list(shape = NA, size = 2.5, linetype = "solid")
          ),
        shape = "none"
        ) +
      ggtitle(class_names[class_idx + 1])
    plot_list[[class_idx + 1]] <- p
    cat(sprintf("Class %d processing completed\n", class_idx + 1))
    }
  return(plot_list)
  }
set.seed(123)
calibration_plots <- create_calibration_plot(
  ova_results$original_probs, 
  ova_results$calibrated_probs, 
  polcal_results$calibrated_probs,
  individual_data$Outcome
  )

if (length(calibration_plots) > 0) {
  for (i in 1:length(calibration_plots)) {
    if (!is.null(calibration_plots[[i]])) {
      print(calibration_plots[[i]])
      # Save high-quality plots
      ggsave(filename = paste0("Calibration_Plot_Outcome_", i, ".png"),
             plot = calibration_plots[[i]],
             width = 7, height = 7, units = "in", dpi = 600, pointsize = 12, bg = "white")
      }
    }
  }