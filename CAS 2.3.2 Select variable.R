# ==============================================================================
# Data imputation
# ==============================================================================
# Methodology: Two-stage variable selection approach
#   Step 1: Multiple imputation + Univariate screening (30 ¡ú 12-15 features)
#   Step 2: LASSO regression + Stepwise selection (15 ¡ú 6 features)
#
# Key features:
#   Robust handling of missing data through multiple imputation
#   Stability analysis across 5 imputed datasets
#   Combined statistical and machine learning approaches
#   Clinical relevance consideration through effect size thresholds
#
# Selection criteria:
#   Univariate: p < 0.05, significant in ¡Ý3/5 imputations
#   Multivariate: Top 6 variables based on combined selection frequency  
#   Final selection: Variables with frequency ¡Ý 0.5 considered clinically relevant
# 
# Actual selection results (based on stability analysis):
#   Age, Hypertension, CEA, Smoking: Perfect stability (1.0)
#   GFR: High stability (0.6)  
#   NEU: Moderate stability (0.5) - included for clinical relevance
# ==============================================================================
### Load data
CAS_clean <- CAS_complete_data[,c(6, 10:39)]

### Convert data types
preprocess_data <- function(data) {
  cat("\n=== Performing data preprocessing ===\n")
  for(i in 2:ncol(data)) {
    data[[i]][data[[i]] == "NA"] <- NA
    data[[i]] <- as.numeric(data[[i]])
    }
  
  ## Handle the outcome variable
  data$Outcome[data$Outcome == "NA"] <- NA
  data$Outcome <- as.factor(data$Outcome)
  
  ## Basic information
  cat("Data dimensions:", nrow(data), "¡Á", ncol(data), "\n")
  cat("Original distribution of the outcome variable:\n")
  print(table(data$Outcome, useNA = "always"))
  
  ## If you wish to modify the reference category, use the following code
  # cat("\nAdjusting the reference category...\n")
  # data$Outcome <- relevel(data$Outcome, ref = "1")
  # cat("Distribution of the outcome variable after adjustment (1 as the reference category):\n")
  # print(table(data$Outcome, useNA = "always"))
  
  ## Missing value situation
  missing_count <- sapply(data, function(x) sum(is.na(x)))
  cat("\nTop 10 variables with the most missing values:\n")
  print(sort(missing_count, decreasing = TRUE)[1:10])
  return(data)
  }
CAS_processed <- preprocess_data(CAS_clean)



# ==============================================================================
# Data driven variable selection
# ==============================================================================
library(mice)   # version 3.17.0
library(nnet)   # version 7.3-20
library(glmnet) # version 4.1-8
library(pROC)   # version 1.18.5

### Univariate analysis
# Purpose: Filter 30+ variables down to 12-15 using univariate significance
# Method: Multiple imputation (m=5) + likelihood ratio tests
# Criteria: Variables significant (p < 0.05) in ¡Ý 3/5 imputations
# Output: Reduced variable set for advanced multivariate selection
step1_robust_imputation_screening <- function(data, m, p_threshold) {
  cat("\n=== Step 1: Multiple imputation + Univariate filtering ===\n")
  cat("Aim: 30 ¡ú 12-15 features\n")
  start_time <- Sys.time()

  ## Handling the Outcome
  cat("Number of raw data samples:", nrow(data), "\n")
  cat("Missing number of Outcome:", sum(is.na(data$Outcome)), "\n")
  data_clean <- data[!is.na(data$Outcome), ]
  cat("Sample size after removing missing Outcome:", nrow(data_clean), "\n")
  cat("Distribution of Outcome:\n")
  print(table(data_clean$Outcome))
  if(!is.factor(data_clean$Outcome)) {
    data_clean$Outcome <- as.factor(data_clean$Outcome)
    }
  
  ##################################
  # Multiple imputation
  ##################################
  cat("\n=== Step 1.1: Perform multiple imputation ===\n")
  imp_data <- mice(data_clean, m, printFlag = F, seed = 123)
  cat("Completed\n")
  all_vars <- names(data_clean)[names(data_clean) != "Outcome"]
  cat("Number of variables to be filtered:", length(all_vars), "\n")
  
  ##################################
  # Univariate filtering
  ##################################
  cat("\n=== Step 1.2: Perform univariate analysis ===\n")
  ## Storage of univariate analysis results
  univar_results <- data.frame(
    Variable = all_vars,
    Mean_P_Value = NA, SD_P_Value = NA, Mean_AbsCoef = NA, 
    Significant_Count = 0, Model_Success_Count = 0,
    stringsAsFactors = FALSE
    )
  pb <- txtProgressBar(min = 0, max = length(all_vars), style = 3)
  for(j in 1:length(all_vars)) {
    var <- all_vars[j]
    p_values <- c()
    coefs <- c()
    success_count <- 0
    for(i in 1:m) {
      complete_data <- complete(imp_data, i)
      tryCatch({
        # Extract and analyze data
        analysis_data <- complete_data[, c("Outcome", var), drop = FALSE]
        # Ensure that there are no missing values in the data
        if(any(is.na(analysis_data))) {next}
        # Check sample size
        if(nrow(analysis_data) < 20) {next}
        # Check variable variability
        if(length(unique(analysis_data[[var]])) < 2) {next}
        # Ensure that Outcome is a factor
        analysis_data$Outcome <- as.factor(analysis_data$Outcome)
        # Remove possible abnormal levels
        analysis_data$Outcome <- droplevels(analysis_data$Outcome)
        # Check factor level
        if(length(levels(analysis_data$Outcome)) < 2) {next}
        # Create data box
        clean_data <- data.frame(
          Outcome = analysis_data$Outcome,
          var_value = as.numeric(analysis_data[[var]]),
          stringsAsFactors = FALSE
          )
        # Recheck data integrity
        if(any(is.na(clean_data)) || any(is.infinite(clean_data$var_value))) {next} 
        # Check extreme values
        if(any(abs(clean_data$var_value) > 1e10)) {next} 
        tryCatch({
          # Fitting a univariate model
          model_uni <- multinom(Outcome ~ var_value, data = clean_data, 
                                trace = FALSE, maxit = 1000, Hess = TRUE)
          # Fitting zero model
          null_model <- multinom(Outcome ~ 1, data = clean_data, 
                                 trace = FALSE, maxit = 1000, Hess = TRUE) 
          # Likelihood ratio test (LRT) of univariate model
          log_lik_full <- logLik(model_uni)   
          # LRT of zero model
          log_lik_null <- logLik(null_model)
          # Calculate LRT
          lr_stat <- -2 * (as.numeric(log_lik_null) - as.numeric(log_lik_full))
          df_diff <- attr(log_lik_full, "df") - attr(log_lik_null, "df")
          if(df_diff > 0 && lr_stat >= 0) {
            p_val <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
            if(!is.na(p_val) && is.finite(p_val) && p_val >= 0 && p_val <= 1) {
              p_values <- c(p_values, p_val)
              success_count <- success_count + 1
              coef_val <- coef(model_uni)
              if(is.matrix(coef_val)) {
                # Multi classification situation
                if(ncol(coef_val) >= 2) {
                  coef_abs <- abs(coef_val[1, 2])
                  } else {
                  coef_abs <- NA
                  }
                } else {
                # Binary classification situation
                if(length(coef_val) >= 2) {
                  coef_abs <- abs(coef_val[2])
                  } else {
                    coef_abs <- NA
                    }
                  }
              if(!is.na(coef_abs) && is.finite(coef_abs)) {
                coefs <- c(coefs, coef_abs)
                }
              }
            }
          }, error = function(e2) {
            # If multinom fails, try other methods
            if(length(levels(clean_data$Outcome)) == 2) {
              # Binary case, use glm
              tryCatch({
                binary_outcome <- as.numeric(clean_data$Outcome) - 1
                glm_model <- glm(binary_outcome ~ var_value, data = clean_data, family = binomial())
                # Use Wald test
                summary_glm <- summary(glm_model)
                # P value for variable value
                p_val <- summary_glm$coefficients[2, 4]
                if(!is.na(p_val) && is.finite(p_val) && p_val >= 0 && p_val <= 1) {
                  p_values <<- c(p_values, p_val)
                  success_count <<- success_count + 1
                  coef_abs <- abs(summary_glm$coefficients[2, 1])
                  if(!is.na(coef_abs) && is.finite(coef_abs)) {
                    coefs <<- c(coefs, coef_abs)
                    }
                  }
                }, error = function(e3) {
                  # GLM also failed, record but continue
                })
              }
            })
        }, error = function(e) {
          cat(sprintf("\nFeature %s, imputation %d processing failed: %s\n", var, i, e$message))
        })
      }
    
    ## Summary the results
    univar_results$Model_Success_Count[j] <- success_count
    if(length(p_values) > 0) {
      univar_results$Mean_P_Value[j] <- mean(p_values, na.rm = TRUE)
      univar_results$SD_P_Value[j] <- sd(p_values, na.rm = TRUE)
      univar_results$Significant_Count[j] <- sum(p_values < p_threshold, na.rm = TRUE)
      }
    if(length(coefs) > 0) {
      univar_results$Mean_AbsCoef[j] <- mean(coefs, na.rm = TRUE)
      }
    setTxtProgressBar(pb, j)
    }
  close(pb)
  
  ##################################
  # Filter valid results
  ##################################
  valid_results <- univar_results[!is.na(univar_results$Mean_P_Value), ]
  cat(sprintf("\nNumber of variables successfully analyzed: %d\n", nrow(valid_results)))
  if(nrow(valid_results) == 0) {
    cat("\nError! All variable model fittings failed\n")
    cat("\n== Diagnostic Information ==\n")
    test_data <- complete(imp_data, 1)
    cat("First 5 rows of the first imputed dataset:\n")
    print(head(test_data, 5))
    cat("\nOutcome variable information:\n")
    print(table(test_data$Outcome, useNA = "ifany"))
    cat("\nInformation of the first predictor variable:\n")
    first_var <- all_vars[1]
    cat("Variable name:", first_var, "\n")
    print(summary(test_data[[first_var]]))
    return(list(
      error = "All variable fittings failed",
      diagnostic_data = test_data[1:10, c("Outcome", first_var)]
      ))
    }
  
  ## Sort by p-value
  valid_results <- valid_results[order(valid_results$Mean_P_Value), ]
  # Variable selection strategy
  primary_selection <- valid_results$Variable[
    valid_results$Mean_P_Value < p_threshold & 
      valid_results$Significant_Count >= 3 &
      valid_results$Model_Success_Count >= 3
    ]
  # Adjust selection quantity
  if(length(primary_selection) > 15) {
    selected_vars <- valid_results$Variable[1:15]
    selection_method <- "Too many variables, select top 15 most significant variables"
    } else if(length(primary_selection) < 8) {
    selected_vars <- valid_results$Variable[1:min(12, nrow(valid_results))]
    selection_method <- "Few eligible variables, select top 12 variables"
    } else {
    selected_vars <- primary_selection
    selection_method <- "Select variables based on significance criteria"
    }
  cat(sprintf("\n%s\n", selection_method))
  cat("\n=== Univariate filtering results summary ===\n")
  result_table <- valid_results[valid_results$Variable %in% selected_vars, ]
  result_table$Selected <- "Yes"
  result_table[, c("Mean_P_Value", "Mean_AbsCoef", "SD_P_Value")] <- 
    round(result_table[, c("Mean_P_Value", "Mean_AbsCoef", "SD_P_Value")], 6)
  print(result_table[, c("Variable", "Mean_P_Value", "Significant_Count", 
                         "Mean_AbsCoef", "Model_Success_Count", "Selected")])
  
  ## Show details of significant variables
  sig_vars <- valid_results[valid_results$Mean_P_Value < p_threshold, ]
  if(nrow(sig_vars) > 0) {
    cat(sprintf("\nSignificant variables (p < %.2f): %d¸ö\n", p_threshold, nrow(sig_vars)))
    for(i in 1:min(10, nrow(sig_vars))) {
      cat(sprintf("  %s: p=%.6f, |coef|=%.4f, success count=%d/5\n", 
                  sig_vars$Variable[i], sig_vars$Mean_P_Value[i], 
                  sig_vars$Mean_AbsCoef[i], sig_vars$Model_Success_Count[i]))
      }
    } else {
      cat(sprintf("\nNo significant variables (p < %.2f)\n", p_threshold))
      }
  end_time <- Sys.time()
  
  cat(sprintf("\n=== Step 1 completed, selected %d variables, time taken: %.2f minutes ===\n", 
              length(selected_vars), as.numeric(difftime(end_time, start_time, units = "mins"))))
  return(list(
    data_clean = data_clean,
    selected_vars = selected_vars,
    results_table = valid_results,
    imputed_data = imp_data,
    selection_method = selection_method,
    significant_vars = if(nrow(sig_vars) > 0) sig_vars$Variable else character(0)
    ))
  }
step1_result <- step1_robust_imputation_screening(CAS_processed, m = 5, p_threshold = 0.05)


### LASSO and stepwise regression
# Purpose: Further reduce 12-15 variables to 6 final predictors
# Methods: 
#   LASSO regression with cross-validation (alpha=1)
#   Stepwise regression using AIC criterion
#   Stability analysis across 5 imputed datasets
# Selection: Top 6 variables based on combined selection frequency
# Note: Final selection includes NEU (frequency 0.5) for clinical relevance
step2_advanced_selection <- function(step1_result, max_vars, lasso_alpha, cv_folds) {
  cat("\n=== Step 2: LASSO regression + Stepwise regression ===\n")
  cat("Aim: 15 ¡ú 6 features\n")
  start_time <- Sys.time()
  
  ## Get data from step1 result
  selected_vars <- step1_result$selected_vars
  imp_data <- step1_result$imputed_data
  cat("Number of variables selected in step 1:", length(selected_vars), "\n")
  cat("Variable list:", paste(selected_vars, collapse = ", "), "\n")
  cat("\nVariable name diagnostics\n")
  test_data <- complete(imp_data, 1)
  cat("Column names in imputed dataset:\n")
  cat(paste(names(test_data), collapse = ", "), "\n")
  all_names <- names(test_data)
  for(var in selected_vars) {
    if(!var %in% all_names) {
      similar_names <- all_names[grepl(var, all_names, ignore.case = TRUE)]
      if(length(similar_names) > 0) {
        cat(sprintf("Variable '%s' not found, possible matches: %s\n", 
                    var, paste(similar_names, collapse = ", ")))
        } else {
          cat(sprintf("Variable '%s' not found and no similar matches\n", var))
          }
      }
    }
  
  # Initialize empty lists for results
  lasso_selections <- list()
  stepwise_selections <- list()
  
  cat("\n==================================\n")
  cat("\nStarting variable selection on 5 imputed datasets\n")
  for(i in 1:5) {
    cat(sprintf("\n--- Imputed dataset %d ---\n", i))
    complete_data <- complete(imp_data, i)
    # Check if variables exist in dataset
    available_vars <- intersect(selected_vars, names(complete_data))
    missing_vars <- setdiff(selected_vars, names(complete_data))
    if(length(missing_vars) > 0) {
      cat("Warning: The following variables are not in the dataset:", paste(missing_vars, collapse = ", "), "\n")
      cat("Will use available", length(available_vars), "variables\n")
      }
    if(length(available_vars) == 0) {
      cat("Error: No available variables, skipping this imputation\n")
      lasso_selections[[i]] <- character(0)
      stepwise_selections[[i]] <- character(0)
      next
      }
    # Prepare data for analysis
    analysis_data <- complete_data[, c("Outcome", available_vars), drop = FALSE]
    analysis_data$Outcome <- as.factor(analysis_data$Outcome)
    analysis_data$Outcome <- droplevels(analysis_data$Outcome)
    # Check data quality
    if(any(is.na(analysis_data))) {
      cat("Warning: Dataset", i, "contains missing values, skipping\n")
      lasso_selections[[i]] <- character(0)
      stepwise_selections[[i]] <- character(0)
      next
      }
    if(length(levels(analysis_data$Outcome)) < 2) {
      cat("Warning: Dataset", i, "has insufficient Outcome levels, skipping\n")
      lasso_selections[[i]] <- character(0)
      stepwise_selections[[i]] <- character(0)
      next
      }
    n_samples <- nrow(analysis_data)
    n_vars <- length(available_vars)
    cat("Samples:", n_samples, "Variables:", n_vars, "\n")

    ##################################
    # LASSO regression
    ##################################
    cat("\n=== Step 2.1: LASSO regression ===\n")
    tryCatch({
      # Prepare data for LASSO
      x_matrix <- as.matrix(analysis_data[, available_vars, drop = FALSE])
      y_vector <- analysis_data$Outcome
      if(any(is.na(x_matrix)) || any(is.infinite(x_matrix))) {
        stop("X matrix contains NA or infinite values")
        } # Data quality check
      if(length(levels(y_vector)) == 2) {
        # Binary classification
        family_type <- "binomial"
        y_numeric <- as.numeric(y_vector) - 1
        } else {
          # Multiclass classification
          family_type <- "multinomial"
          y_numeric <- y_vector
          }
      # Cross-validation for optimal lambda
      set.seed(123 + i)
      cv_lasso <- cv.glmnet(x_matrix, y_numeric, family = family_type, alpha = lasso_alpha,
                            nfolds = min(cv_folds, n_samples), type.measure = "class")
      optimal_lambda <- cv_lasso$lambda.min # Use lambda.min
      lasso_model <- glmnet(x_matrix, y_numeric, family = family_type,
                            alpha = lasso_alpha, lambda = optimal_lambda) # Fit final LASSO model
      # Extract non-zero coefficients
      if(family_type == "multinomial") {
        coefs <- coef(lasso_model)[[1]]
        } else {
          coefs <- coef(lasso_model)
          }
      nonzero_indices <- which(coefs[-1] != 0)
      lasso_selected <- available_vars[nonzero_indices]
      cat("LASSO selected variables:", length(lasso_selected), "\n")
      if(length(lasso_selected) > 0) {
        cat("LASSO selected variables:", paste(lasso_selected, collapse = ", "), "\n")
        lasso_selections[[i]] <- lasso_selected
        } else {
          cat("LASSO selected no variables\n")
          lasso_selections[[i]] <- character(0)
          }
      }, error = function(e) {
        cat("LASSO regression failed:", e$message, "\n")
        lasso_selections[[i]] <- character(0)
        })
    
    ##################################
    # Stepwise regression
    ##################################
    cat("\n=== Step 2.2: Stepwise regression ===\n")
    stepwise_selections[[i]] <- character(0)
    tryCatch({
      step_data <- analysis_data
      current_vars <- available_vars
      is_binary <- length(levels(step_data$Outcome)) == 2
      if(length(current_vars) <= 3) {
        cat("Few variables, keeping all\n")
        stepwise_selections[[i]] <- current_vars
        } else {
          cat("Sufficient variables, performing AIC backward stepwise selection...\n")
          full_formula <- as.formula(paste("Outcome ~", paste(current_vars, collapse = " + ")))
          if(is_binary) {
            full_model <- glm(full_formula, data = step_data, family = binomial())
            } else {
              full_model <- multinom(full_formula, data = step_data, trace = FALSE, maxit = 1000)
              }
          current_aic <- AIC(full_model)
          best_vars <- current_vars
          
          ## Stepwise variable removal
          improved <- TRUE
          iteration <- 0
          max_iterations <- length(current_vars) - 1  # Keep at least 1 variable
          while(improved && length(current_vars) > 1 && iteration < max_iterations) {
            improved <- FALSE
            iteration <- iteration + 1
            best_test_aic <- current_aic
            best_test_vars <- current_vars
            # Try removing each variable
            for(j in 1:length(current_vars)) {
              test_vars <- current_vars[-j]
              if(length(test_vars) == 0) next
              test_formula <- as.formula(paste("Outcome ~", paste(test_vars, collapse = " + ")))
              tryCatch({
                if(is_binary) {
                  test_model <- glm(test_formula, data = step_data, family = binomial())
                  } else {
                    test_model <- multinom(test_formula, data = step_data, 
                                           trace = FALSE, maxit = 1000)
                    }
                test_aic <- AIC(test_model)
                # If AIC improves, record this configuration
                if(test_aic < best_test_aic) {
                  best_test_aic <- test_aic
                  best_test_vars <- test_vars
                  improved <- TRUE
                  }
                }, error = function(e_inner) {
                  # Skip failed models
                  })
              }
            # Update current best configuration
            if(improved) {
              current_vars <- best_test_vars
              current_aic <- best_test_aic
              best_vars <- current_vars
              }
            }
          stepwise_selections[[i]] <- best_vars
          }
      cat("Stepwise selected variables:", length(stepwise_selections[[i]]), "\n")
      if(length(stepwise_selections[[i]]) > 0) {
        cat("Stepwise selected variables:", paste(stepwise_selections[[i]], collapse = ", "), "\n")
        } else {
          cat("Stepwise selected no variables\n")
          }
      }, error = function(e) {
        cat("Stepwise regression failed:", e$message, "\n")
      
      ## Fallback: use LASSO results or top variables
      if(length(lasso_selections[[i]]) > 0) {
        cat("Using LASSO results as fallback\n")
        stepwise_selections[[i]] <- lasso_selections[[i]]
        } else {
          cat("Using top 3 variables as fallback\n")
          stepwise_selections[[i]] <- available_vars[1:min(3, length(available_vars))]
          }
        })
    }
  
  ##################################
  # Summary and stability analysis
  ##################################
  cat("\n=== Step 2.3: Variable selection stability analysis ===\n")
  cat("Selection results across imputed datasets:\n")
  for(i in 1:5) {
    cat(sprintf("Imputation %d - LASSO: %s\n", i, ifelse(length(lasso_selections[[i]]) > 0, 
                                                         paste(lasso_selections[[i]], 
                                                               collapse = ", "), "None")))
    cat(sprintf("Imputation %d - Stepwise: %s\n", i, ifelse(length(stepwise_selections[[i]]) > 0, 
                                                            paste(stepwise_selections[[i]], 
                                                                  collapse = ", "), "None")))
    }
  
  ## Count selection frequency for each variable
  all_lasso_vars <- unique(unlist(lasso_selections[lengths(lasso_selections) > 0]))
  all_stepwise_vars <- unique(unlist(stepwise_selections[lengths(stepwise_selections) > 0]))
  all_possible_vars <- unique(c(all_lasso_vars, all_stepwise_vars))
  if(length(all_possible_vars) == 0) {
    cat("Warning: No variables selected!\n")
    return(list(
      final_vars = character(0), stability_table = data.frame(),
      lasso_selections = lasso_selections, stepwise_selections = stepwise_selections,
      selection_strategy = "No variables selected", performance_results = list(),
      imputed_data = imp_data
      ))
    }
  # Calculate selection frequencies
  lasso_counts <- sapply(all_possible_vars, function(var) {
    sum(sapply(lasso_selections, function(sel) length(sel) > 0 && var %in% sel))
    })
  stepwise_counts <- sapply(all_possible_vars, function(var) {
    sum(sapply(stepwise_selections, function(sel) length(sel) > 0 && var %in% sel))
    })
  combined_counts <- lasso_counts + stepwise_counts
  
  ## Create stability summary table
  stability_table <- data.frame(
    Variable = all_possible_vars, LASSO_Count = lasso_counts,
    Stepwise_Count = stepwise_counts, Combined_Count = combined_counts,
    LASSO_Frequency = round(lasso_counts / 5, 3), 
    Stepwise_Frequency = round(stepwise_counts / 5, 3),
    Combined_Frequency = round(combined_counts / 10, 3),
    stringsAsFactors = FALSE
    )
  # Sort by combined frequency
  stability_table <- stability_table[order(stability_table$Combined_Count, decreasing = TRUE), ]
  cat("\nVariable selection stability summary:\n")
  print(stability_table)
  
  ##################################
  # Final variable selection strategy
  ##################################
  cat("\n=== Step 2.4: Final variable selection ===\n")
  # Strategy 1: Combined frequency >= 0.6
  # Note: Select top 6 variables based on combined frequency
  # NEU (frequency 0.5) is included due to clinical relevance for inflammation assessment
  high_stability_vars <- stability_table$Variable[stability_table$Combined_Frequency >= 0.6]
  # Strategy 2: Supported by both methods with combined frequency >= 0.4
  both_methods_vars <- stability_table$Variable[
    stability_table$LASSO_Frequency >= 0.4 & 
      stability_table$Stepwise_Frequency >= 0.4 &
      stability_table$Combined_Frequency >= 0.4
    ]
  # Strategy 3: Highly supported by one method (frequency >= 0.8)
  highly_supported_vars <- stability_table$Variable[
    (stability_table$LASSO_Frequency >= 0.8 | stability_table$Stepwise_Frequency >= 0.8) &
      stability_table$Combined_Frequency >= 0.3
    ]
  # Strategy 4: If above strategies select too few variables, use combined frequency
  candidate_vars <- unique(c(high_stability_vars, both_methods_vars, highly_supported_vars))
  
  if(length(candidate_vars) == 0) {
    # If no variables meet strict criteria, relax criteria
    min_freq <- max(0.2, min(stability_table$Combined_Frequency[1:min(6, nrow(stability_table))]))
    candidate_vars <- stability_table$Variable[stability_table$Combined_Frequency >= min_freq]
    candidate_vars <- candidate_vars[1:min(max_vars, length(candidate_vars))]
    selection_strategy <- paste("Based on relaxed criteria (frequency >=", min_freq, ")")
    } else if(length(candidate_vars) > max_vars) {
      # If too many candidate variables, select highest frequency ones
      candidate_indices <- which(stability_table$Variable %in% candidate_vars)
      candidate_vars <- stability_table$Variable[candidate_indices[1:max_vars]]
      selection_strategy <- paste("Based on multiple criteria (top", max_vars, "variables)")
      } else {
        selection_strategy <- "Comprehensive selection based on multiple criteria"
        }
  final_vars <- candidate_vars
  cat("Selection strategy:", selection_strategy, "\n")
  cat("Final number of selected variables:", length(final_vars), "\n")
  cat("Final variable list:", paste(final_vars, collapse = ", "), "\n")
  
  ## Show stability details for final variables
  if(length(final_vars) > 0) {
    cat("\nStability details for final variables:\n")
    final_stability <- stability_table[stability_table$Variable %in% final_vars, ]
    print(final_stability[, c("Variable", "LASSO_Count", "Stepwise_Count", 
                              "Combined_Count", "Combined_Frequency")])
    }
  
  ##################################
  # Performance evaluation
  ##################################
  cat("\n=== Step 2.5: Model performance evaluation ===\n")
  performance_results <- list()
  if(length(final_vars) > 0) {
    for(i in 1:5) {
      complete_data <- complete(imp_data, i)
      eval_data <- complete_data[, c("Outcome", final_vars), drop = FALSE]
      eval_data$Outcome <- as.factor(eval_data$Outcome)
      eval_data$Outcome <- droplevels(eval_data$Outcome)
      tryCatch({
        if(length(levels(eval_data$Outcome)) == 2) {
          # Binary classification performance evaluation
          formula_eval <- as.formula(paste("Outcome ~", paste(final_vars, collapse = " + ")))
          eval_model <- glm(formula_eval, data = eval_data, family = binomial())
          # Calculate AUC
          pred_prob <- predict(eval_model, type = "response")
          roc_obj <- roc(eval_data$Outcome, pred_prob, quiet = TRUE)
          auc_value <- as.numeric(auc(roc_obj))
          performance_results[[i]] <- list(
            AIC = AIC(eval_model),
            AUC = auc_value,
            Variables = length(final_vars)
            )
          cat(sprintf("Imputation %d: AIC=%.2f, AUC=%.3f, Variables=%d\n", 
                      i, AIC(eval_model), auc_value, length(final_vars)))
          } else {
          # Multiclass performance evaluation
            formula_eval <- as.formula(paste("Outcome ~", paste(final_vars, collapse = " + ")))
            eval_model <- multinom(formula_eval, data = eval_data, trace = FALSE)
            performance_results[[i]] <- list(
              AIC = AIC(eval_model),
              Variables = length(final_vars)
              )
            cat(sprintf("Imputation %d: AIC=%.2f, Variables=%d\n", 
                        i, AIC(eval_model), length(final_vars)))
            }
        }, error = function(e) {
          cat(sprintf("Imputation %d performance evaluation failed: %s\n", i, e$message))
          })
      }
    
    ## Calculate average performance
    if(length(performance_results) > 0) {
      avg_aic <- mean(sapply(performance_results, function(x) x$AIC), na.rm = TRUE)
      cat(sprintf("\nAverage AIC: %.2f\n", avg_aic))
      if(length(levels(complete(imp_data, 1)$Outcome)) == 2) {
        auc_values <- sapply(performance_results, function(x) x$AUC)
        if(any(!is.na(auc_values))) {
          avg_auc <- mean(auc_values, na.rm = TRUE)
          cat(sprintf("Average AIC: %.3f\n", avg_auc))
          }
        }
      }
    }
  end_time <- Sys.time()
  cat(sprintf("\n=== Step 2 completed, selected %d variables, time taken: %.2f minutes ===\n", 
              length(final_vars), as.numeric(difftime(end_time, start_time, units = "mins"))))
  return(list(
    final_vars = final_vars,
    stability_table = stability_table,
    lasso_selections = lasso_selections,
    stepwise_selections = stepwise_selections,
    selection_strategy = selection_strategy,
    performance_results = performance_results,
    imputed_data = imp_data
    ))
  }
step2_result <- step2_advanced_selection(step1_result, max_vars = 6, lasso_alpha = 1, cv_folds = 10)
