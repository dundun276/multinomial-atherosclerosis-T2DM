# ==============================================================================
# Bootstrap optimism assessment framework
# ==============================================================================
# Completely independent development:
# This bootstrap validation system was developed from scratch based on:
#   Clinical requirements for robust internal validation
#   Need for comprehensive performance assessment (discrimination + calibration)
#
# Key innovations:
#   Stratified protected bootstrap sampling for minority classes
#   Consistent multiple imputation across bootstrap samples
#   Class-weighted multinomial regression integration  
#   Modular architecture for maintainability and extensibility
#   Comprehensive performance metrics integration
# ==============================================================================
library(mice)     # version 3.17.0
library(nnet)     # version 7.3-20
library(dplyr)    # version 1.1.4
library(pROC)     # version 1.18.5

set.seed(123)

### Hierarchical bootstrap sampling
stratified_bootstrap_protected <- function(data, outcome_col) {
  outcome_table <- table(data[[outcome_col]])
  bootstrap_indices <- c()
  for (level in names(outcome_table)) {
    level_indices <- which(data[[outcome_col]] == level)
    original_size <- length(level_indices)
    target_size <- max(original_size, 10)
    sampled_indices <- sample(level_indices, size = target_size, replace = TRUE)
    bootstrap_indices <- c(bootstrap_indices, sampled_indices)
    }
  return(data[bootstrap_indices, ])
  }


### Consistency interpolation
apply_consistent_imputation <- function(new_data, original_mice_obj, m = 5) {
  original_method <- original_mice_obj$method
  original_pred_matrix <- original_mice_obj$predictorMatrix
  new_imp <- mice(new_data, 
                  m = m, 
                  method = original_method,
                  predictorMatrix = original_pred_matrix,
                  printFlag = FALSE,
                  seed = sample(1:10000, 1))
  return(new_imp)
  }


### Weighted multinomial regression fitting
fit_weighted_multinomial <- function(mice_obj, class_weights = NULL) {
  if (is.null(class_weights)) {
    stacked_data <- complete(mice_obj, action = "long", include = FALSE)
    outcome_counts <- table(stacked_data$Outcome)
    class_weights <- sum(outcome_counts) / (length(outcome_counts) * outcome_counts)
    # Clean up temporary data immediately
    rm(stacked_data)  
    gc(verbose = FALSE)
    }
  # Fit weighted polynomial logistic regression model
  weighted_model <- with(mice_obj, {
    current_data <- data.frame(
      Outcome = Outcome,
      Age = Age, 
      Smoking = Smoking, 
      Hypertension = Hypertension, 
      CEA = CEA, 
      GFR = GFR, 
      NEU = NEU
      )
    weights <- rep(1, nrow(current_data))
    for(level in names(class_weights)) {
      weights[current_data$Outcome == level] <- class_weights[level]
      }
    multinom(Outcome ~ Age + Smoking + Hypertension + CEA + GFR + NEU, 
             weights = weights,
             model = FALSE,
             trace = FALSE,
             maxit = 1000)
    })
  # Clean memory
  gc(verbose = FALSE)
  return(weighted_model)
  }


### Get model predictions
get_model_predictions <- function(mice_obj, model_obj) {
  ## Get pooled results
  pooled_results <- pool(model_obj)
  pooled_sum <- summary(pooled_results)
  n_vars <- 7
  n_outcomes <- 2
  pooled_sum <- pooled_sum[1:(n_outcomes * n_vars), ]
  pooled_sum$outcome_level <- rep(1:n_outcomes, each = n_vars)
  pooled_sum$comparison <- rep(c("Outcome2_vs_Outcome1", "Outcome3_vs_Outcome1"), each = n_vars)
  
  ## Extract coefficients
  coeffs_cat2 <- pooled_sum %>%
    filter(comparison == "Outcome2_vs_Outcome1") %>%
    select(term, estimate)
  coeffs_cat3 <- pooled_sum %>%
    filter(comparison == "Outcome3_vs_Outcome1") %>%
    select(term, estimate)
  
  ## Build design matrix
  stacked_data <- complete(mice_obj, action = "long", include = FALSE)
  form <- as.formula(Outcome ~ Age + Smoking + Hypertension + CEA + GFR + NEU)
  design_mat <- model.matrix(form, data = stacked_data)
  # Calculate LP
  linear_pred_2 <- as.numeric(cbind(1, design_mat[, coeffs_cat2$term[-1]]) %*% coeffs_cat2$estimate)
  linear_pred_3 <- as.numeric(cbind(1, design_mat[, coeffs_cat3$term[-1]]) %*% coeffs_cat3$estimate)
  # Calculate probabilities
  pred_results <- stacked_data %>%
    filter(.imp != 0) %>%
    mutate(
      lp_cat2 = linear_pred_2,
      lp_cat3 = linear_pred_3
    ) %>%
    bind_cols(calculate_multinomial_probs(linear_pred_2, linear_pred_3))
  # Calculate individual-level means
  prob_and_lp_means <- pred_results %>%
    group_by(.id) %>%
    summarise(
      Outcome = unique(Outcome),
      .groups = 'drop',
      prob_avg_cat1 = mean(prob_cat1), 
      prob_avg_cat2 = mean(prob_cat2), 
      prob_avg_cat3 = mean(prob_cat3),
      lp_avg_cat2 = mean(lp_cat2), 
      lp_avg_cat3 = mean(lp_cat3)
      )
  # Create probability matrix
  prob_matrix <- prob_and_lp_means %>%
    select(prob_avg_cat1, prob_avg_cat2, prob_avg_cat3) %>%
    as.matrix()
  
  return(list(
    individual_data = prob_and_lp_means,
    prob_matrix = prob_matrix,
    outcomes = prob_and_lp_means$Outcome
    ))
  }


### Evaluate calibration parameters
evaluate_calibration_parameters <- function(pred_results) {
  calibration_params <- list()
  
  ## Original model calibration parameters
  original_cal_stats <- calculate_calibration_stats(pred_results$prob_matrix, pred_results$outcomes)
  
  ## Store calibration parameters for each class
  for (class_idx in 0:2) {
    class_key <- paste0("Outcome_", class_idx)
    if (class_key %in% names(original_cal_stats)) {
      calibration_params[[paste0("original_slope_", class_idx + 1)]] <- original_cal_stats[[class_key]]$slope
      calibration_params[[paste0("original_intercept_", class_idx + 1)]] <- original_cal_stats[[class_key]]$intercept
      } else {
        calibration_params[[paste0("original_slope_", class_idx + 1)]] <- NA
        calibration_params[[paste0("original_intercept_", class_idx + 1)]] <- NA
        }
    }
  return(calibration_params)
  }


### Evaluate comprehensive performance
evaluate_comprehensive_performance <- function(mice_obj, model_obj) {
  ## Get prediction results
  pred_results <- get_model_predictions(mice_obj, model_obj)
  
  performance_metrics <- list()
  ## Calculate PDI
  tryCatch({
    pdi_result <- calculate_pdi_3cat(pred_results$outcomes, pred_results$prob_matrix)
    performance_metrics$pdi <- pdi_result$pdi_mean
    }, error = function(e) {
      performance_metrics$pdi <- NA
      })
  
  ## Calculate C-statistics
  prob_mean_cat2 <- pred_results$individual_data %>%
    select(.id, Outcome, prob_avg_cat2) %>%
    rename(prob_avg = prob_avg_cat2)
  prob_mean_cat3 <- pred_results$individual_data %>%
    select(.id, Outcome, prob_avg_cat3) %>%
    rename(prob_avg = prob_avg_cat3)
  tryCatch({
    sink(tempfile())
    cstat_result_2vs1 <- calculate_cstat(prob_mean_cat2, 2, "Submodel 1")
    sink()
    performance_metrics$cstat_2vs1 <- as.numeric(cstat_result_2vs1$auc)
    }, error = function(e) {
      sink()
      performance_metrics$cstat_2vs1 <- NA
      })
  tryCatch({
    sink(tempfile())
    cstat_result_3vs1 <- calculate_cstat(prob_mean_cat3, 1, "Submodel 2")
    sink()
    performance_metrics$cstat_3vs1 <- as.numeric(cstat_result_3vs1$auc)
    }, error = function(e) {
      sink()
      performance_metrics$cstat_3vs1 <- NA
      })
  
  ## Calculate nagelkerke R2
  tryCatch({
    sink(tempfile())
    r2_result <- calculate_nagelkerke_r2(pred_results$individual_data)
    sink()
    performance_metrics$nagelkerke_r2 <- r2_result$nagelkerke
    }, error = function(e) {
    sink()
      performance_metrics$nagelkerke_r2 <- NA
      })
  
  ## Calculate brier Score
  tryCatch({
    brier_result <- calculate_brier(pred_results$prob_matrix, pred_results$outcomes)
    performance_metrics$brier_score <- brier_result$overall_brier
    }, error = function(e) {
      performance_metrics$brier_score <- NA
      })
  
  ## Calibration parameters evaluation
  calibration_metrics <- evaluate_calibration_parameters(pred_results)
  
  ## Combine all metrics
  all_metrics <- c(performance_metrics, calibration_metrics)
  
  return(list(
    metrics = all_metrics,
    pred_results = pred_results
    ))
  }


### Main bootstrap analysis
bootstrap_optimism_analysis <- function(original_data, original_mice, B) {
  cat("\n=== Performing bootstrap optimism assessment analysis ===\n")
  cat(sprintf("Number of Bootstrap iterations: %d\n", B))

  ## Calculate class weights for the original data
  stacked_data <- complete(original_mice, action = "long", include = FALSE)
  outcome_counts <- table(stacked_data$Outcome)
  class_weights <- sum(outcome_counts) / (length(outcome_counts) * outcome_counts)
  cat("Original data class distribution:\n")
  print(outcome_counts)
  cat("Outcome weights:\n")
  print(round(class_weights, 3))

  ##################################
  # Evaluate original data performance
  ##################################
  cat("\n=== Evaluating original data eerformance ===\n")
  original_model <- fit_weighted_multinomial(original_mice, class_weights)
  original_performance <- evaluate_comprehensive_performance(original_mice, original_model)
  cat("Original data performance metrics:\n")
  cat(sprintf("  PDI: %.4f\n", original_performance$metrics$pdi))
  cat(sprintf("  C-stat (2vs1): %.4f\n", original_performance$metrics$cstat_2vs1))
  cat(sprintf("  C-stat (3vs1): %.4f\n", original_performance$metrics$cstat_3vs1))
  cat(sprintf("  Nagelkerke R2: %.4f\n", original_performance$metrics$nagelkerke_r2))
  cat(sprintf("  Brier Score: %.4f\n", original_performance$metrics$brier_score))

  ##################################
  # Bootstrap loop
  ##################################
  cat("\n=== Starting bootstrap sampling ===\n")

  ## Initialize results storage
  bootstrap_results <- data.frame(
    iteration = 1:B,
    pdi = rep(NA, B),
    cstat_2vs1 = rep(NA, B),
    cstat_3vs1 = rep(NA, B),
    nagelkerke_r2 = rep(NA, B),
    brier_score = rep(NA, B),
    original_slope_1 = rep(NA, B),
    original_intercept_1 = rep(NA, B),
    original_slope_2 = rep(NA, B),
    original_intercept_2 = rep(NA, B),
    original_slope_3 = rep(NA, B),
    original_intercept_3 = rep(NA, B),
    success = rep(FALSE, B)
    )
  successful_iterations <- 0
  for (b in 1:B) {
    if (b %% 50 == 0) {
      progress_pct <- round(100 * b / B, 1)
      progress_bar <- paste(rep("¨€", floor(progress_pct/2)), collapse = "")
      progress_space <- paste(rep("???", 50 - floor(progress_pct/2)), collapse = "")
      cat(sprintf("\rProgress: %d/%d [%s%s] %.1f%% (Successful: %d, Success Rate: %.1f%%)", 
                  b, B, progress_bar, progress_space, progress_pct, 
                  successful_iterations, 100 * successful_iterations / max(1, b-1)))
      if (b == B) cat("\n")
      }
    tryCatch({
      # Stratified Bootstrap sampling
      boot_data <- stratified_bootstrap_protected(original_data, "Outcome")
      # Check the sampled data
      boot_table <- table(boot_data$Outcome)
      if (any(boot_table < 5)) {
        stop(paste("Insufficient sample size for some classes:", paste(boot_table, collapse = ", ")))
        }
      # Apply consistent imputation
      boot_mice <- apply_consistent_imputation(boot_data, original_mice, m = 5)
      # Check if imputation was successful
      if (is.null(boot_mice) || boot_mice$m == 0) {
        stop("Imputation failed")
        }
      # Fit weighted model
      boot_model <- fit_weighted_multinomial(boot_mice, class_weights)
      # Check if model fitting was successful
      if (is.null(boot_model) || length(boot_model$analyses) == 0) {
        stop("Model fitting failed")
        }
      # Evaluate performance
      boot_performance <- evaluate_comprehensive_performance(boot_mice, boot_model)
      # Check if performance evaluation was successful
      if (is.null(boot_performance) || is.null(boot_performance$metrics)) {
        stop("Performance evaluation failed")
        }

      ## Store all results
      bootstrap_results[b, "pdi"] <- boot_performance$metrics$pdi
      bootstrap_results[b, "cstat_2vs1"] <- boot_performance$metrics$cstat_2vs1
      bootstrap_results[b, "cstat_3vs1"] <- boot_performance$metrics$cstat_3vs1
      bootstrap_results[b, "nagelkerke_r2"] <- boot_performance$metrics$nagelkerke_r2
      bootstrap_results[b, "brier_score"] <- boot_performance$metrics$brier_score
      bootstrap_results[b, "original_slope_1"] <- boot_performance$metrics$original_slope_1
      bootstrap_results[b, "original_intercept_1"] <- boot_performance$metrics$original_intercept_1
      bootstrap_results[b, "original_slope_2"] <- boot_performance$metrics$original_slope_2
      bootstrap_results[b, "original_intercept_2"] <- boot_performance$metrics$original_intercept_2
      bootstrap_results[b, "original_slope_3"] <- boot_performance$metrics$original_slope_3
      bootstrap_results[b, "original_intercept_3"] <- boot_performance$metrics$original_intercept_3
      bootstrap_results[b, "success"] <- TRUE
      successful_iterations <- successful_iterations + 1
      }, error = function(e) {
        if (b <= 10) {  # # Report detailed information for the first 10 failures
          cat(sprintf("Bootstrap %d failed: %s\n", b, as.character(e$message)))
          }
        })
    }
  success_rate <- successful_iterations / B
  cat(sprintf("\nBootstrap completed£¡\n"))
  cat(sprintf("Overall success rate: %d/%d (%.1f%%)\n", successful_iterations, B, 100 * success_rate))

  ##################################
  # Calculate optimism and corrected performance
  ##################################
  cat("\n=== Calculating optimism ===\n")

  valid_results <- bootstrap_results[bootstrap_results$success, ]
  original_cal_params <- evaluate_calibration_parameters(original_performance$pred_results)
  performance_metrics <- c("pdi", "cstat_2vs1", "cstat_3vs1", "nagelkerke_r2", "brier_score")
  calibration_metrics <- c("original_slope_1", "original_intercept_1",
                           "original_slope_2", "original_intercept_2",
                           "original_slope_3", "original_intercept_3")
  all_metrics <- c(performance_metrics, calibration_metrics)

  ## Create summary table
  metrics_summary <- data.frame(
    Metric = c("PDI", "C-stat (2vs1)", "C-stat (3vs1)", "Nagelkerke R2", "Brier Score",
               "Cal Slope (Class 1)", "Cal Intercept (Class 1)",
               "Cal Slope (Class 2)", "Cal Intercept (Class 2)",
               "Cal Slope (Class 3)", "Cal Intercept (Class 3)"),
    Original = c(
      original_performance$metrics$pdi,
      original_performance$metrics$cstat_2vs1,
      original_performance$metrics$cstat_3vs1,
      original_performance$metrics$nagelkerke_r2,
      original_performance$metrics$brier_score,
      original_cal_params$original_slope_1,
      original_cal_params$original_intercept_1,
      original_cal_params$original_slope_2,
      original_cal_params$original_intercept_2,
      original_cal_params$original_slope_3,
      original_cal_params$original_intercept_3
      ),
    Bootstrap_Mean = c(
      mean(valid_results$pdi, na.rm = TRUE),
      mean(valid_results$cstat_2vs1, na.rm = TRUE),
      mean(valid_results$cstat_3vs1, na.rm = TRUE),
      mean(valid_results$nagelkerke_r2, na.rm = TRUE),
      mean(valid_results$brier_score, na.rm = TRUE),
      mean(valid_results$original_slope_1, na.rm = TRUE),
      mean(valid_results$original_intercept_1, na.rm = TRUE),
      mean(valid_results$original_slope_2, na.rm = TRUE),
      mean(valid_results$original_intercept_2, na.rm = TRUE),
      mean(valid_results$original_slope_3, na.rm = TRUE),
      mean(valid_results$original_intercept_3, na.rm = TRUE)
      ),
    stringsAsFactors = FALSE
    )
  
  # Calculate optimism
  metrics_summary$Optimism <- metrics_summary$Bootstrap_Mean - metrics_summary$Original
  # Calculate corrected performance
  metrics_summary$Corrected <- metrics_summary$Original - metrics_summary$Optimism
  # Special handling for Brier Score
  brier_index <- which(metrics_summary$Metric == "Brier Score")
  if (length(brier_index) > 0) {
    if (metrics_summary$Optimism[brier_index] > 0) {
      metrics_summary$Corrected[brier_index] <- metrics_summary$Original[brier_index] + abs(metrics_summary$Optimism[brier_index])
      }
    }

  return(list(
    original_metrics = c(original_performance$metrics, original_cal_params),
    bootstrap_results = bootstrap_results,
    valid_results = valid_results,
    metrics_summary = metrics_summary,
    success_rate = success_rate,
    sample_sizes = list(
      original = nrow(original_data),
      bootstrap_mean = nrow(original_data)
      )
    ))
  }


### Result reporting function
report_bootstrap_results <- function(bootstrap_analysis) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("Final report of bootstrap optimism assessment\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  cat(sprintf("\nBootstrap Settings:\n"))
  cat(sprintf("  Total bootstrap iterations: %d\n", nrow(bootstrap_analysis$bootstrap_results)))
  cat(sprintf("  Successful iterations: %d\n", nrow(bootstrap_analysis$valid_results)))
  cat(sprintf("  Success rate: %.1f%%\n", 100 * bootstrap_analysis$success_rate))

  cat(sprintf("\n=== Discrimination performance metrics comparison ===\n"))
  summary_table <- bootstrap_analysis$metrics_summary

  # Discrimination performance metrics
  performance_indices <- 1:5
  perf_table <- summary_table[performance_indices, ]

  cat(sprintf("%-20s %10s %10s %10s %10s\n",
              "Metric", "Original", "Bootstrap", "Optimism", "Corrected"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in 1:nrow(perf_table)) {
    cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n",
                perf_table$Metric[i],
                perf_table$Original[i],
                perf_table$Bootstrap_Mean[i],
                perf_table$Optimism[i],
                perf_table$Corrected[i]))
  }

  cat(sprintf("\n=== Calibration parameters comparison ===\n"))

  # Calibration parameters
  calibration_indices <- 6:11
  cal_table <- summary_table[calibration_indices, ]

  cat(sprintf("%-25s %10s %10s %10s %10s\n",
              "Calibration Parameter", "Original", "Bootstrap", "Optimism", "Corrected"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (i in 1:nrow(cal_table)) {
    if (!is.na(cal_table$Original[i])) {
      cat(sprintf("%-25s %10.4f %10.4f %10.4f %10.4f\n",
                  cal_table$Metric[i],
                  cal_table$Original[i],
                  cal_table$Bootstrap_Mean[i],
                  cal_table$Optimism[i],
                  cal_table$Corrected[i]))
    } else {
      cat(sprintf("%-25s %10s %10s %10s %10s\n",
                  cal_table$Metric[i], "NA", "NA", "NA", "NA"))
    }
  }

  # Calculate 95% confidence intervals for each metric
  cat(sprintf("\n=== Bootstrap 95%% confidence intervals ===\n"))
  valid_data <- bootstrap_analysis$valid_results

  # Discrimination metrics CI
  cat("Discrimination metrics:\n")
  perf_metrics <- c("pdi", "cstat_2vs1", "cstat_3vs1", "nagelkerke_r2", "brier_score")
  perf_names <- c("PDI", "C-stat (2vs1)", "C-stat (3vs1)", "Nagelkerke R2", "Brier Score")

  for (i in 1:length(perf_metrics)) {
    metric_data <- valid_data[[perf_metrics[i]]]
    if (any(!is.na(metric_data))) {
      ci_lower <- quantile(metric_data, 0.025, na.rm = TRUE)
      ci_upper <- quantile(metric_data, 0.975, na.rm = TRUE)
      cat(sprintf("%-20s [%.4f, %.4f]\n", perf_names[i], ci_lower, ci_upper))
    }
  }

  # Calibration parameters CI
  cat("\nCalibration parameters:\n")
  cal_metrics <- c("original_slope_1", "original_intercept_1",
                   "original_slope_2", "original_intercept_2",
                   "original_slope_3", "original_intercept_3")
  cal_names <- c("Cal Slope (Class 1)", "Cal Intercept (Class 1)",
                 "Cal Slope (Class 2)", "Cal Intercept (Class 2)",
                 "Cal Slope (Class 3)", "Cal Intercept (Class 3)")

  for (i in 1:length(cal_metrics)) {
    metric_data <- valid_data[[cal_metrics[i]]]
    if (any(!is.na(metric_data))) {
      ci_lower <- quantile(metric_data, 0.025, na.rm = TRUE)
      ci_upper <- quantile(metric_data, 0.975, na.rm = TRUE)
      cat(sprintf("%-25s [%.4f, %.4f]\n", cal_names[i], ci_lower, ci_upper))
    }
  }
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("Notes:\n")
  cat("- Optimism = Bootstrap Mean - Original\n")
  cat("- Corrected = Original - Optimism (Special handling for Brier Score)\n")
  cat("- Calibration Parameters: A slope close to 1.0 and an intercept close to 0.0 indicate good calibration\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  }


### Main execution function
run_bootstrap_optimism_evaluation <- function(original_data, original_mice, B) {
  # Perform bootstrap analysis
  results <- bootstrap_optimism_analysis(
    original_data = original_data,
    original_mice = original_mice,
    B = B
    )
  report_bootstrap_results(results)
  return(results)
  }

### Execute analysis
small_results <- run_bootstrap_optimism_evaluation(step1_result$data_clean, 
                                                   step1_result$imputed_data,
                                                   B = 5)
# final_results <- run_bootstrap_optimism_evaluation(step1_result$data_clean, 
#                                                    step1_result$imputed_data,
#                                                    B = 1000)

# ================================================================================ 
#   Final report of bootstrap optimism assessment
# ================================================================================ 
#   
#   Bootstrap Settings:
# Total bootstrap iterations: 1000
# Successful iterations: 1000
# Success rate: 100.0%
# 
# === Discrimination performance metrics comparison ===
#   Metrics               Original  Bootstrap  Optimism  Corrected
# ---------------------------------------------------------------------- 
# PDI                      0.5662     0.5777     0.0115     0.5547
# C-stat (2vs1)            0.7477     0.7537     0.0060     0.7418
# C-stat (3vs1)            0.7910     0.7883    -0.0027     0.7937
# Nagelkerke R2            0.3500     0.3574     0.0074     0.3426
# Brier Score              0.1839     0.1829    -0.0009     0.1848
# 
# === Calibration Parameter ===
#   Calibration Parameter      Original  Bootstrap  Optimism  Corrected
# --------------------------------------------------------------------------- 
# Cal Slope (Class 1)           1.0362     1.0298    -0.0064     1.0426
# Cal Intercept (Class 1)       0.4222     0.4135    -0.0087     0.4309
# Cal Slope (Class 2)           1.0738     1.0855     0.0117     1.0622
# Cal Intercept (Class 2)      -1.5766    -1.5690     0.0076    -1.5842
# Cal Slope (Class 3)           1.2374     1.2056    -0.0318     1.2692
# Cal Intercept (Class 3)       0.9505     0.9232    -0.0273     0.9778
# 
# === Bootstrap 95%% confidence intervals ===
#   Discrimination metrics:
# PDI                  [0.5452, 0.6122]
# C-stat (2vs1)        [0.6955, 0.8015]
# C-stat (3vs1)        [0.7557, 0.8189]
# Nagelkerke R2        [0.3020, 0.4130]
# Brier Score          [0.1759, 0.1900]
# 
#   Calibration Parameter:
# Cal Slope (Class 1)       [0.9635, 1.1053]
# Cal Intercept (Class 1)   [0.3353, 0.5003]
# Cal Slope (Class 2)       [0.9404, 1.2933]
# Cal Intercept (Class 2)   [-1.6633, -1.4408]
# Cal Slope (Class 3)       [1.0912, 1.3795]
# Cal Intercept (Class 3)   [0.8199, 1.0480]