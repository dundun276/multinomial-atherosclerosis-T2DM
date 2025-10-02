# ==============================================================================
# Get stacked imputed data
# ==============================================================================
# Base methodology framework:
#   Primary reference: Gehringer et al. J Clin Epidemiol. 2024;174:111481. 
#   Code from: Gehringer CK. multinomial_risk_prediction_RA GitHub repository:
#              https://github.com/celinagehringer/multinomial_risk_prediction_RA
#              (File: Section3.2_MLR_development.R)
#
# Independent innovations:
#   Class weighting for imbalanced data
#   Multinomial nomogram with dual-outcome visualization  
#   Custom result organization and prediction system
#
# Adapted by: Yongqi Zheng | Contact: dundun276@126.com
# Usage: For academic non-commercial research only. Citation required.
# ==============================================================================
library(mice)  # version 3.17.0
library(nnet)  # version 7.3-20
library(dplyr) # version 1.1.4

### Impute data set and stack datasets
# Methodology adapted from Gehringer et al. (2024)
multi_imp <- step1_result$imputed_data
stacked_data_og <- complete(multi_imp, action = "long")
outcome_counts <- table(stacked_data_og$Outcome)
multi_mo <- with(multi_imp,
                 multinom(Outcome ~ Age + Smoking + Hypertension + CEA + GFR + NEU, model = TRUE))


### Modify the multinomial regression fit and add weights
## Weighting is used to improve the predictive ability of the model for minority classes.
# Intuitively speaking, the intercept that affects minority classes is relatively low.
# A lower intercept will lower the initial predicted value, 
# making it more difficult for individuals to classify as a minority class.
cat("Original category distribution:\n", outcome_counts)
class_weights <- sum(outcome_counts) / (length(outcome_counts) * outcome_counts)
cat("\nBalanced weights:\n", class_weights)

multi_mo_weighted <- with(multi_imp, {
  current_data <- data.frame(
    Outcome = Outcome,
    Age = Age, 
    Smoking = Smoking, 
    Hypertension = Hypertension, 
    CEA = CEA, 
    GFR = GFR, 
    NEU = NEU
    )
  # Create weights
  weights <- rep(1, nrow(current_data))
  for(level in names(class_weights)) {
    weights[current_data$Outcome == level] <- class_weights[level]
    }
  # Fit the weighted model
  multinom(Outcome ~ Age + Smoking + Hypertension + CEA + GFR + NEU, 
           weights = weights,
           model = TRUE,
           trace = FALSE,
           maxit = 1000)
  })


### Summary of imputation dataset results
# Methodology adapted from Gehringer et al. (2024)
# Custom function developed to better organize multinomial model results
custom_pool_summary <- function(pool_obj) {
  summary_df <- summary(pool_obj)
  # Intercept and 6 variables 
  n_vars <- 7
  # There are only 2 comparisons: 2vs1 and 3vs1
  n_outcomes <- 2
  summary_df <- summary_df[1:(n_outcomes * n_vars), ]
  summary_df$outcome_level <- rep(1:n_outcomes, each = n_vars)
  summary_df$comparison <- rep(c("Outcome2_vs_Outcome1", "Outcome3_vs_Outcome1"), each = n_vars)
  return(summary_df)
  }
pooled_sum <- custom_pool_summary(multi_mo_weighted) %>%
  mutate(
    Odds.ratio = exp(estimate),
    Confint_lower = exp(estimate - 1.96 * std.error),
    Confint_upper = exp(estimate + 1.96 * std.error)
    )


### Display Model Summary
display_model_coefficients <- function(pooled_sum) {
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Multinomial Logistic Regression Model Coefficients\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  summary_table <- pooled_sum %>%
    mutate(
      Comparison = case_when(
        comparison == "Outcome2_vs_Outcome1" ~ "Outcome 2 vs 1",
        comparison == "Outcome3_vs_Outcome1" ~ "Outcome 3 vs 1",
        TRUE ~ comparison
        ),
      Coefficient = round(estimate, 4),
      SE = round(std.error, 4),
      P_value = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ as.character(round(p.value, 3))
        ),
      OR = round(Odds.ratio, 3),
      CI_95 = paste0("(", round(Confint_lower, 3), "-", round(Confint_upper, 3), ")")
      ) %>%
    select(Comparison, Variable = term, Coefficient, SE, P_value, OR, CI_95)
  for(comp in unique(summary_table$Comparison)) {
    cat("\n", comp, ":\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")
    comp_data <- summary_table %>% filter(Comparison == comp)
    print(comp_data %>% select(-Comparison))
    }
  cat("\nNote: OR = Odds Ratio, CI_95 = 95% Confidence Interval\n")
  }
display_model_coefficients(pooled_sum)



# ==============================================================================
# RMS-based multinomial logistic regression nomogram
# ==============================================================================
# Completely independent development:
# This nomogram implementation was developed from scratch based on:
#   6 variables selected through data-driven selection in previous analysis
#   Clinical requirements for carotid atherosclerosis prediction
#
# Key innovations:
#   Dual-line visualization for multinomial outcomes
#   Integrated linear predictor scales
#   Custom scoring system with clinical relevance
# ==============================================================================
library(rms)   # version 6.8-0
library(mice)  # version 3.17.0
library(nnet)  # version 7.3-20
library(ggplot2)   # version 3.5.1
library(gridExtra) # version 2.3
library(grid)  # version 4.2.3
library(dplyr) # version 1.1.4
library(RColorBrewer) # version 1.1-3

### Data Preparation Function
# Variables selected through data-driven selection in previous analysis
prepare_nomogram_data <- function(pooled_sum, imputed_data) {
  # Extract coefficients for Outcome 2 vs 1
  coefs_outcome_2vs1 <- pooled_sum %>%
    filter(comparison == "Outcome2_vs_Outcome1") %>%
    select(term, estimate, std.error, p.value, Odds.ratio, Confint_lower, Confint_upper)
  # Extract coefficients for Outcome 3 vs 1  
  coefs_outcome_3vs1 <- pooled_sum %>%
    filter(comparison == "Outcome3_vs_Outcome1") %>%
    select(term, estimate, std.error, p.value, Odds.ratio, Confint_lower, Confint_upper)
  # Get sample data for variable ranges
  sample_data <- complete(imputed_data, 1)
  
  ## Customize settings for continuous and categorical variables
  # Define variable ranges based on actual data
  var_ranges <- list(
    Age = seq(floor(min(sample_data$Age, na.rm = TRUE) / 5) * 5, 
              ceiling(max(sample_data$Age, na.rm = TRUE) / 5) * 5, by = 5),
    Smoking = c(0, 1),
    Hypertension = c(0, 1),
    CEA = seq(0, ceiling(max(sample_data$CEA, na.rm = TRUE) / 5) * 5, by = 5),
    GFR = seq(0, ceiling(max(sample_data$GFR, na.rm = TRUE)), by = 1),
    NEU = seq(0, ceiling(max(sample_data$NEU, na.rm = TRUE)), by = 1)
    )
  # Variable information
  var_info <- list(
    Age = list(label = "Age", unit = "years", type = "continuous"),
    Smoking = list(label = "Smoking status", unit = "", type = "binary",
                   labels = c("No" = 0, "Yes" = 1)),
    Hypertension = list(label = "Hypertension", unit = "", type = "binary", 
                        labels = c("No" = 0, "Yes" = 1)),
    CEA = list(label = "CEA", unit = "ng/mL", type = "continuous"),
    GFR = list(label = "GFR", unit = "mL/min", type = "continuous"),
    NEU = list(label = "NEU", unit = expression(paste("¡Á10"^"9", "/L")), type = "continuous")
    )
  return(list(
    coefs_outcome1 = coefs_outcome_2vs1,
    coefs_outcome2 = coefs_outcome_3vs1,
    var_ranges = var_ranges,
    var_info = var_info,
    sample_data = sample_data
    ))
  }


### Create Individual Binary Nomograms for Each Outcome
create_binary_nomogram_data <- function(coefs, var_ranges, var_info, outcome_name) {
  # Create a mock dataset for rms
  n_points <- 1000
  mock_data <- data.frame(
    Outcome = rbinom(n_points, 1, 0.5)  # Binary outcome for this comparison
    )
  # Add predictor variables
  for(var in names(var_ranges)) {
    if(var_info[[var]]$type == "continuous") {
      mock_data[[var]] <- runif(n_points, 
                                min(var_ranges[[var]]), 
                                max(var_ranges[[var]]))
      } else {
        mock_data[[var]] <- sample(var_ranges[[var]], n_points, replace = TRUE)
        }
    }
  # Set up rms environment
  dd <- datadist(mock_data)
  options(datadist = 'dd')
  # Create formula
  formula_str <- paste("Outcome ~", paste(names(var_ranges), collapse = " + "))
  # Fit binary logistic regression model
  model <- lrm(as.formula(formula_str), data = mock_data, x = TRUE, y = TRUE)
  # Replace coefficients with our pooled estimates
  # Find matching coefficients
  model_coef_names <- names(coef(model))
  for(i in 1:length(model_coef_names)) {
    coef_name <- model_coef_names[i]
    if(coef_name == "Intercept") {
      matching_coef <- coefs$estimate[coefs$term == "(Intercept)"]
      } else {
        matching_coef <- coefs$estimate[coefs$term == coef_name]
        }
    if(length(matching_coef) > 0) {
      model$coefficients[i] <- matching_coef
      }
    }
  return(list(model = model, dd = dd))
  }


### Function to determine optimal tick marks
get_optimal_ticks <- function(var_range, scores, var_info, max_ticks = 8) {
  if(var_info$type == "binary") {
    return(1:length(var_range))
    }
  score_range <- max(scores) - min(scores)
  if(score_range < 0.1) {
    return(1)
    } else if(score_range < 10) {
      tick_indices <- c(1, length(var_range))
      if(length(var_range) > 4) {
        mid_point <- ceiling(length(var_range)/2)
        tick_indices <- c(1, mid_point, length(var_range))
        }
      } else if(score_range < 30) {
        n_ticks <- min(5, length(var_range))
        tick_indices <- round(seq(1, length(var_range), length.out = n_ticks))
        } else {
          n_ticks <- min(max_ticks, length(var_range))
          tick_indices <- round(seq(1, length(var_range), length.out = n_ticks))
          }
  return(unique(tick_indices))
  }


### Create Enhanced Dual-Line Nomogram with Dual Total Score and LP Lines
create_dual_line_nomogram <- function(pooled_sum, imputed_data, colors, width, height) {
  ##################################
  # Data Preparation
  ##################################
  nom_data <- prepare_nomogram_data(pooled_sum, imputed_data)
  var_names <- names(nom_data$var_ranges)
  
  ## Create binary model data
  # Submodel 1: Normal group vs. Thickening group
  Submodel1_data <- create_binary_nomogram_data(
    nom_data$coefs_outcome1, nom_data$var_ranges, 
    nom_data$var_info, "Outcome 2 vs 1"
    )
  # Submodel 2: Normal group vs. Plaque group
  Submodel2_data <- create_binary_nomogram_data(
    nom_data$coefs_outcome2, nom_data$var_ranges, 
    nom_data$var_info, "Outcome 3 vs 1" 
    )
  
  ##################################
  # Calculation Functions
  ##################################
  ## Calculate maximum contribution for each outcome
  calc_max_contribution <- function(coefs, var_ranges) {
    max_contributions <- sapply(var_names, function(var) {
      var_range <- var_ranges[[var]]
      var_min <- min(var_range)
      var_max <- max(var_range)
      coef <- coefs$estimate[coefs$term == var]
      coef <- if(length(coef) == 0) 0 else as.numeric(coef)
      abs((var_max - var_min) * coef)
      })
    return(max(max_contributions, na.rm = TRUE))
    }
  
  ## Calculate separate total score ranges
  calculate_score_ranges <- function() {
    max_contrib_1 <- calc_max_contribution(nom_data$coefs_outcome1, nom_data$var_ranges)
    max_contrib_2 <- calc_max_contribution(nom_data$coefs_outcome2, nom_data$var_ranges)
    global_max_contrib <- max(max_contrib_1, max_contrib_2)
    calc_total_range <- function(coefs) {
      total_min <- total_max <- 0
      for(var in var_names) {
        var_range <- nom_data$var_ranges[[var]]
        var_min <- min(var_range)
        var_max <- max(var_range)
        coef <- coefs$estimate[coefs$term == var]
        coef <- if(length(coef) == 0) 0 else as.numeric(coef)
        scores <- (var_range - var_min) * coef
        min_score <- min(scores)
        if(min_score < 0) scores <- scores - min_score
        if(global_max_contrib > 0) {
          normalized_scores <- (scores / global_max_contrib) * 100
          } else {
            normalized_scores <- rep(0, length(scores))
            }
        total_min <- total_min + min(normalized_scores)
        total_max <- total_max + max(normalized_scores)
        }
      return(list(min_score = total_min, max_score = total_max))
      }
    return(list(
      outcome1 = calc_total_range(nom_data$coefs_outcome1),
      outcome2 = calc_total_range(nom_data$coefs_outcome2),
      max_contribution_outcome1 = max_contrib_1,
      max_contribution_outcome2 = max_contrib_2,
      global_max_contribution = global_max_contrib
      ))
    }
  
  ## Calculate linear predictor (LP) ranges
  calculate_lp_ranges <- function() {
    calc_lp_range <- function(coefs_data) {
      intercept <- coefs_data$estimate[coefs_data$term == "(Intercept)"]
      lp_min <- lp_max <- intercept
      for(var in var_names) {
        var_range <- nom_data$var_ranges[[var]]
        var_min <- min(var_range)
        var_max <- max(var_range)
        coef_raw <- coefs_data$estimate[coefs_data$term == var]
        coef <- if(length(coef_raw) == 0) 0 else as.numeric(coef_raw)
        if(coef >= 0) {
          lp_min <- lp_min + coef * var_min
          lp_max <- lp_max + coef * var_max
          } else {
            lp_min <- lp_min + coef * var_max
            lp_max <- lp_max + coef * var_min
            }
        }
      return(list(min = lp_min, max = lp_max))
      }
    return(list(
      outcome1 = calc_lp_range(nom_data$coefs_outcome1),
      outcome2 = calc_lp_range(nom_data$coefs_outcome2)
      ))
    }
  
  ##################################
  # Drawing Helper Functions
  ##################################
  ## Extract coefficient information
  extract_coeff_info <- function(var) {
    coef_1_raw <- nom_data$coefs_outcome1$estimate[nom_data$coefs_outcome1$term == var]
    coef_2_raw <- nom_data$coefs_outcome2$estimate[nom_data$coefs_outcome2$term == var]
    or_1 <- nom_data$coefs_outcome1$Odds.ratio[nom_data$coefs_outcome1$term == var]
    or_2 <- nom_data$coefs_outcome2$Odds.ratio[nom_data$coefs_outcome2$term == var]
    return(list(
      coef_1 = if(length(coef_1_raw) == 0) 0 else as.numeric(coef_1_raw),
      coef_2 = if(length(coef_2_raw) == 0) 0 else as.numeric(coef_2_raw),
      or_1 = if(length(or_1) == 0) 1 else as.numeric(or_1),
      or_2 = if(length(or_2) == 0) 1 else as.numeric(or_2)
      ))
    }
  
  ## Calculate variable score positions
  calculate_position <- function(var_range, coef, global_max_contrib) {
    var_min <- min(var_range)
    scores <- (var_range - var_min) * coef
    min_score <- min(scores)
    if(min_score < 0) scores <- scores - min_score
    if(global_max_contrib > 0) {
      x_positions <- (scores / global_max_contrib) * 100
      } else {
        x_positions <- rep(0, length(var_range))
        }
    return(pmax(x_positions, 0))
    }
  
  ## Get variable tick marks
  get_variable_ticks <- function(var, var_range) {
    # When you want to customize the display
    if (var == "Age") {
      return(which(var_range %in% c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)))
      } else if (var == "CEA") {
        return(which(var_range %in% c(0, 5, 10, 15, 20)))
        } else if (var == "GFR") {
          return(which(var_range %in% c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 155)))
          } else if (var == "NEU") {
            # If not, only use the following 3 lines of code
            var_length <- length(var_range)
            middle_idx <- round(var_length / 2)
            return(c(1, middle_idx, var_length))
            } else {
              # For other variables, return evenly distributed tick marks
              return(seq(1, length(var_range), length.out = min(6, length(var_range))))
              }
    }
  
  ## Draw variable lines
  draw_variable_line <- function(x_positions, y_pos, color, var_range, var_info, 
                                 all_ticks, or_value, var_name, cex_scale, line_offset = 0) {
    y_line <- y_pos + line_offset
    is_or_one <- abs(or_value - 1) < 0.001
    if(is_or_one) {
      points(0, y_line, pch = 19, col = color, cex = cex_scale)
      # Regarding the almost no statistical effect treatment (OR = 1.000) of GFR in this model
      if(var_name == "GFR" && line_offset < 0) {
        text(5, y_line, "OR = 1.000", cex = 0.9, col = color)
        text(0, y_line - 0.25, as.character(var_range[1]), cex = cex_scale, col = color)
        } else {
          text_y <- if(line_offset > 0) y_line + 0.25 else y_line - 0.25
          text(0, text_y, "OR = 1.000", cex = 0.9, col = color)
          opposite_y <- if(line_offset > 0) y_line - 0.25 else y_line + 0.25
          text(0, opposite_y, as.character(var_range[1]), cex = cex_scale, col = color)
          }
      } else {
        line_start <- max(0, min(x_positions))
        line_end <- max(x_positions)
        if(line_end > line_start) {
          segments(line_start, y_line, line_end, y_line, col = color, lwd = 3)
          } else {
            segments(line_start - 2, y_line, line_start + 2, y_line, col = color, lwd = 3)
            }
        for(j in all_ticks) {
          if(j <= length(var_range)) {
            x_pos <- max(0, x_positions[j])
            if(x_pos >= 0 && x_pos <= 100) {
              segments(x_pos, y_line - 0.05, x_pos, y_line + 0.05, col = color, lwd = 2)
              value_label <- if(var_info$type == "binary" && !is.null(var_info$labels)) {
                names(var_info$labels)[var_info$labels == var_range[j]]
                } else {
                  as.character(var_range[j])
                  }
              text_y <- if(line_offset > 0) y_line + 0.25 else y_line - 0.25
              # Special handling for the CP line of the GFR: only show the first and last labels
              # Because it has too many labels and the effect is too narrow
              show_label <- TRUE
              if(var_name == "GFR" && line_offset < 0) {
                show_label <- (j == all_ticks[1] || j == all_ticks[length(all_ticks)])
                }
              if(show_label) {
                text(x_pos, text_y, value_label, cex = cex_scale, col = color)
                }
              }
            }
          }
        }
    }
  
  ## Generate tick marks
  generate_ticks <- function(min_score, max_score) {
    range_size <- max_score - min_score
    step <- if(range_size <= 10) 2 
    else if(range_size <= 50) 5 
    else if(range_size <= 100) 10 
    else if(range_size <= 200) 20 
    else 25
    ticks <- seq(0, ceiling(max_score), by = step)
    return(sort(unique(c(0, ticks, max_score))))
    }
  
  ## Draw scale tick marks
  draw_scale_ticks <- function(y_pos, color, ticks, max_value, cex_scale, label_offset = 0.25) {
    for(score in ticks) {
      x_pos <- (score / max_value) * 100
      if(x_pos >= 0 && x_pos <= 100) {
        segments(x_pos, y_pos - 0.05, x_pos, y_pos + 0.05, col = color, lwd = 2)
        text(x_pos, y_pos + label_offset, sprintf("%.0f", score), cex = cex_scale, col = color)
        }
      }
    }
  
  ##################################
  # Main Drawing Functions
  ##################################
  ## Setup plot environment
  setup_plot <- function(y_positions, n_vars, cex_main) {
    par(mfrow = c(1, 1), mar = c(0, 0, 2.5, 1))
    plot(0, 0, type = "n", xlim = c(-30, 108), ylim = c(0, n_vars + 3.65),
         xlab = "", ylab = "", axes = FALSE,
         main = "Multinomial Logistic Regression Nomogram",
         cex.main = cex_main, font.main = 2)
    }
  
  ## Draw points scale
  draw_points_scale <- function(y_points, cex_scale, cex_label) {
    points_seq <- seq(0, 100, 5)
    for(i in seq_along(points_seq)) {
      x_pos <- points_seq[i]
      segments(x_pos, y_points - 0.1, x_pos, y_points + 0.1, lwd = 2)
      if(points_seq[i] %% 10 == 0) {
        text(x_pos, y_points - 0.3, points_seq[i], cex = cex_scale + 0.1)
        }
      }
    segments(0, y_points, 100, y_points, lwd = 3)
    text(-25, y_points, "Points", adj = 0, cex = cex_label, font = 2)
    }
  
  ## Draw all variables
  draw_all_variables <- function(y_positions, separate_total_ranges, cex_scale, cex_label) {
    for(i in seq_along(var_names)) {
      var <- var_names[i]
      y_var <- y_positions[i + 1]
      var_info <- nom_data$var_info[[var]]
      var_range <- nom_data$var_ranges[[var]]
      # Get coefficient information
      coeff_info <- extract_coeff_info(var)
      # Format variable label
      if(var == "NEU" && is.expression(var_info$unit)) {
        var_label <- paste0(var_info$label, " (¡Á10^9/L)")
        } else {
          var_label <- paste0(var_info$label,
                              if(var_info$unit != "" && !is.expression(var_info$unit)) 
                                paste0(" (", var_info$unit, ")") else "")
          }
      text(-25, y_var, var_label, adj = 0, cex = cex_label, font = 2)
      # Calculate positions
      x_pos_1 <- calculate_position(var_range, coeff_info$coef_1, 
                                    separate_total_ranges$global_max_contribution)
      x_pos_2 <- calculate_position(var_range, coeff_info$coef_2, 
                                    separate_total_ranges$global_max_contribution)
      # Get tick marks
      all_ticks <- get_variable_ticks(var, var_range)
      # Draw lines
      draw_variable_line(x_pos_1, y_var, colors[1], var_range, var_info, 
                         all_ticks, coeff_info$or_1, var, cex_scale, 0.1)
      draw_variable_line(x_pos_2, y_var, colors[2], var_range, var_info, 
                         all_ticks, coeff_info$or_2, var, cex_scale, -0.1)
      }
    }
  
  ## Draw total score scales
  draw_total_scores <- function(y_total, separate_total_ranges, cex_scale, cex_label) {
    y_total_1 <- y_total + 0.08
    y_total_2 <- y_total - 0.08
    # Draw lines
    segments(0, y_total_1, 100, y_total_1, col = colors[1], lwd = 3)
    segments(0, y_total_2, 100, y_total_2, col = colors[2], lwd = 3)
    # Generate and draw tick marks
    total_ticks_1 <- generate_ticks(separate_total_ranges$outcome1$min_score, 
                                    separate_total_ranges$outcome1$max_score)
    total_ticks_2 <- generate_ticks(separate_total_ranges$outcome2$min_score, 
                                    separate_total_ranges$outcome2$max_score)
    draw_scale_ticks(y_total_1, colors[1], total_ticks_1, 
                     separate_total_ranges$outcome1$max_score, cex_scale, 0.25)
    draw_scale_ticks(y_total_2, colors[2], total_ticks_2, 
                     separate_total_ranges$outcome2$max_score, cex_scale, -0.25)
    # Labels
    text(-25, y_total, "Total Score", adj = 0, cex = cex_label, font = 2, col = "black")
    }
  
  ## Draw LP scales
  draw_lp_scales <- function(y_lp, lp_ranges, cex_scale, cex_label) {
    y_lp_1 <- y_lp + 0.08
    y_lp_2 <- y_lp - 0.08
    # Draw lines
    segments(0, y_lp_1, 100, y_lp_1, col = colors[1], lwd = 3)
    segments(0, y_lp_2, 100, y_lp_2, col = colors[2], lwd = 3)
    # Special function for generating LP tick marks
    generate_lp_ticks <- function(min_val, max_val, target_count = 5) {
      range_size <- max_val - min_val
      if(target_count <= 2) return(c(min_val, max_val))
      inner_count <- target_count - 2
      if(inner_count > 0) {
        step_size <- range_size / (target_count - 1)
        uniform_positions <- seq(min_val + step_size, max_val - step_size, length.out = inner_count)
        adjusted_inner <- sapply(uniform_positions, function(x) {
          candidates <- c(floor(x), ceiling(x))
          if(range_size > 8) {
            candidates <- c(candidates, floor(x * 2) / 2, ceiling(x * 2) / 2)
            }
          candidates <- candidates[abs(candidates - min_val) > 0.5 & abs(candidates - max_val) > 0.5]
          if(length(candidates) > 0) {
            return(candidates[which.min(abs(candidates - x))])
            } else {
              return(round(x))
              }
          })
        # Special handling
        if(min_val < -6 && max_val > 4) {
          suggested_inner <- c(-4, -2, 0, 2, 4, 6)
          adjusted_inner <- suggested_inner[suggested_inner > min_val & suggested_inner < max_val]
          } else if(min_val < -3 && max_val > 6) {
            suggested_inner <- c(-2, 0, 2, 4, 6)
            adjusted_inner <- suggested_inner[suggested_inner > min_val & suggested_inner < max_val]
            }
        all_ticks <- c(min_val, adjusted_inner, max_val)
        return(sort(unique(all_ticks)))
        } else {
          return(c(min_val, max_val))
          }
      }
    # Generate tick marks
    lp_ticks_1 <- generate_lp_ticks(lp_ranges$outcome1$min, lp_ranges$outcome1$max, 5)
    lp_ticks_2 <- generate_lp_ticks(lp_ranges$outcome2$min, lp_ranges$outcome2$max, 5)
    # Draw tick marks
    lp_range_1 <- lp_ranges$outcome1$max - lp_ranges$outcome1$min
    lp_range_2 <- lp_ranges$outcome2$max - lp_ranges$outcome2$min
    for(lp_val in lp_ticks_1) {
      x_pos <- ((lp_val - lp_ranges$outcome1$min) / lp_range_1) * 100
      if(x_pos >= -5 && x_pos <= 105) {
        x_pos_clipped <- pmax(0, pmin(100, x_pos))
        segments(x_pos_clipped, y_lp_1 - 0.05, x_pos_clipped, y_lp_1 + 0.05, col = colors[1], lwd = 2)
        text(x_pos_clipped, y_lp_1 + 0.25, sprintf("%.2f", lp_val), cex = cex_scale, col = colors[1])
        }
      }
    for(lp_val in lp_ticks_2) {
      x_pos <- ((lp_val - lp_ranges$outcome2$min) / lp_range_2) * 100
      if(x_pos >= -5 && x_pos <= 105) {
        x_pos_clipped <- pmax(0, pmin(100, x_pos))
        segments(x_pos_clipped, y_lp_2 - 0.05, x_pos_clipped, y_lp_2 + 0.05, col = colors[2], lwd = 2)
        text(x_pos_clipped, y_lp_2 - 0.25, sprintf("%.2f", lp_val), cex = cex_scale, col = colors[2])
        }
      }
    # Labels
    text(-25, y_lp, "Linear Predictor", adj = 0, cex = cex_label, font = 2, col = "black")
    }
  
  ## Add legend
  add_legend <- function() {
    legend("topright", legend = c("CIMT thickening", "CP presence"),
           col = colors, lwd = 3, cex = 1.15, bty = "n")
    }
  
  ##################################
  # Main Drawing Function
  ##################################
  draw_custom_nomogram <- function() {
    ## Calculate data
    separate_total_ranges <- calculate_score_ranges()
    lp_ranges <- calculate_lp_ranges()
    
    ## Set parameters
    cex_main <- 1.9
    cex_label <- 1.3
    cex_scale <- 1.1
    n_vars <- length(nom_data$var_ranges)
    
    ## Calculate positions
    base_positions <- seq(n_vars + 3, 1, by = -1.06)
    y_positions <- base_positions
    # NEU position adjustment
    neu_index <- which(names(nom_data$var_ranges) == "NEU")
    if(length(neu_index) > 0) {
      neu_pos_index <- neu_index + 1
      for(i in 1:(length(base_positions) - neu_pos_index)) {
        y_positions[neu_pos_index + i] <- y_positions[neu_pos_index + i] - 0.5
        }
      }
    
    ## Execute drawing
    setup_plot(y_positions, n_vars, cex_main)
    draw_points_scale(y_positions[1], cex_scale, cex_label)
    draw_all_variables(y_positions, separate_total_ranges, cex_scale, cex_label)
    add_legend()
    draw_total_scores(1.6, separate_total_ranges, cex_scale, cex_label)
    draw_lp_scales(0.6, lp_ranges, cex_scale, cex_label)
    
    ## Print information
    cat("Outcome 2 Total Score Range: ", separate_total_ranges$outcome1$min_score, 
        " to ", separate_total_ranges$outcome1$max_score, "\n")
    cat("Outcome 3 Total Score Range: ", separate_total_ranges$outcome2$min_score, 
        " to ", separate_total_ranges$outcome2$max_score, "\n")
    cat("Max contribution Outcome 2: ", separate_total_ranges$max_contribution_outcome1, "\n")
    cat("Max contribution Outcome 3: ", separate_total_ranges$max_contribution_outcome2, "\n")
    return(separate_total_ranges)
    }
  
  ##################################
  # Execute Drawing
  ##################################
  png("CAS_nomogram.png", width = 14, height = 7, units = "in", res = 600, pointsize = 12)
  separate_total_ranges <- draw_custom_nomogram()
  dev.off()
  return(list(
    nomogram_data = nom_data,
    Submodel1 = Submodel1_data$model,
    Submodel2 = Submodel1_data$model,
    separate_total_ranges = separate_total_ranges
    ))
  }


### Enhanced Prediction Function
predict_from_nomogram <- function(nomogram_result, patient_data) {
  coefs_2 <- nomogram_result$nomogram_data$coefs_outcome1
  coefs_3 <- nomogram_result$nomogram_data$coefs_outcome2
  # Calculate linear predictors
  lp_2 <- coefs_2$estimate[coefs_2$term == "(Intercept)"]
  lp_3 <- coefs_3$estimate[coefs_3$term == "(Intercept)"]
  # Initialize variables
  individual_scores_2 <- c()
  individual_scores_3 <- c()
  total_score_2 <- 0
  total_score_3 <- 0
  # Get normalization factor from nomogram result
  max_contribution_overall <- nomogram_result$separate_total_ranges$global_max_contribution
  for(var in names(patient_data)) {
    coef_2 <- coefs_2$estimate[coefs_2$term == var]
    coef_3 <- coefs_3$estimate[coefs_3$term == var]
    if(length(coef_2) > 0 && !is.na(coef_2)) {
      lp_2 <- lp_2 + coef_2 * patient_data[[var]]
      # Calculate nomogram score for this variable (Outcome 2)
      var_ranges <- nomogram_result$nomogram_data$var_ranges
      if(var %in% names(var_ranges)) {
        var_range <- var_ranges[[var]]
        var_min <- min(var_range)
        raw_score_2 <- (patient_data[[var]] - var_min) * coef_2
        # Handle negative coefficients
        if(coef_2 < 0) {
          var_max <- max(var_range)
          max_raw_score <- (var_max - var_min) * coef_2
          raw_score_2 <- raw_score_2 - max_raw_score
          }
        nomogram_score_2 <- (abs(raw_score_2) / max_contribution_overall) * 100
        individual_scores_2 <- c(individual_scores_2, nomogram_score_2)
        total_score_2 <- total_score_2 + nomogram_score_2
        }
      }
    if(length(coef_3) > 0 && !is.na(coef_3)) {
      lp_3 <- lp_3 + coef_3 * patient_data[[var]]
      # Calculate nomogram score for this variable (Outcome 3)
      if(var %in% names(var_ranges)) {
        var_range <- var_ranges[[var]]
        var_min <- min(var_range)
        raw_score_3 <- (patient_data[[var]] - var_min) * coef_3
        # Handle negative coefficients
        if(coef_3 < 0) {
          var_max <- max(var_range)
          max_raw_score <- (var_max - var_min) * coef_3
          raw_score_3 <- raw_score_3 - max_raw_score
          }
        nomogram_score_3 <- (abs(raw_score_3) / max_contribution_overall) * 100
        individual_scores_3 <- c(individual_scores_3, nomogram_score_3)
        total_score_3 <- total_score_3 + nomogram_score_3
        }
      }
    }
  # Calculate probabilities using softmax
  exp_lp2 <- exp(lp_2)
  exp_lp3 <- exp(lp_3)
  prob_1 <- 1 / (1 + exp_lp2 + exp_lp3)
  prob_2 <- exp_lp2 / (1 + exp_lp2 + exp_lp3)
  prob_3 <- exp_lp3 / (1 + exp_lp2 + exp_lp3)
  # Create main results data frame
  results <- data.frame(
    Outcome_1_Probability = round(prob_1, 4),
    Outcome_2_Probability = round(prob_2, 4),
    Outcome_3_Probability = round(prob_3, 4),
    Linear_Predictor_2vs1 = round(lp_2, 4),
    Linear_Predictor_3vs1 = round(lp_3, 4),
    Nomogram_Score_Outcome2 = round(total_score_2, 2),
    Nomogram_Score_Outcome3 = round(total_score_3, 2)
    )
  # Add individual variable contributions only if they exist
  if(length(individual_scores_2) > 0) {
    names(individual_scores_2) <- paste0(names(patient_data), "_Score_Outcome2")
    individual_scores_2_df <- data.frame(t(individual_scores_2))
    results <- cbind(results, individual_scores_2_df)
    }
  if(length(individual_scores_3) > 0) {
    names(individual_scores_3) <- paste0(names(patient_data), "_Score_Outcome3")
    individual_scores_3_df <- data.frame(t(individual_scores_3))
    results <- cbind(results, individual_scores_3_df)
    }
  return(results)
  }


### Main Execution
cat("Creating RMS-based Multinomial Nomogram...\n")
nomogram_result <- create_dual_line_nomogram(
  pooled_sum, 
  multi_imp,
  colors = c("#1f77b4","#FF6666"),
  # The following are the colors available for selection
  # colors = c("#1f77b4", "#ff7f0e"),
  # colors = c("#1f77b4", "#AF7AC5"),
  # colors = c("#1f77b4", "#5F9E93"),
  # colors = c("#5F9E93","#FF6666"),
  width = 15, height = 12
  )


### Example prediction
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Prediction Example\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
example_patient <- list(
  Age = 69,
  Smoking = 1,
  Hypertension = 1,
  CEA = 4.71,
  GFR = 80.17,
  NEU = 3.89
  )
cat("\nPrediction Results:\n")
predicted_results <- predict_from_nomogram(nomogram_result, example_patient)
print(predicted_results)
# Outcome_1_Probability     Outcome_2_Probability     Outcome_3_Probability  Linear_Predictor_2vs1
#         0.0278                   0.3278                      0.6444                2.4664
# Linear_Predictor_3vs1    Nomogram_Score_Outcome2   Nomogram_Score_Outcome3   Age_Score_Outcome2
#         3.1424                   136.69                      116.56              66.68135
# Smoking_Score_Outcome2 Hypertension_Score_Outcome2    CEA_Score_Outcome2     GFR_Score_Outcome2
#       11.03292                 17.09091                    6.782323               30.5214
# NEU_Score_Outcome2         Age_Score_Outcome3       Smoking_Score_Outcome3
#       4.580774                   70.000                    15.92051
# Hypertension_Score_Outcome3 CEA_Score_Outcome3         GFR_Score_Outcome3    NEU_Score_Outcome3
#       13.43346                 8.974955                    3.068011              5.164907

example_patient <- list(
  Age = 60,
  Smoking = 0,
  Hypertension = 1,
  CEA = 1.83,
  GFR = 103.42,
  NEU = 3.85
  )
# Outcome_1_Probability     Outcome_2_Probability     Outcome_3_Probability  Linear_Predictor_2vs1
#         0.1463                   0.5139                      0.3398                1.2567
# Linear_Predictor_3vs1    Nomogram_Score_Outcome2   Nomogram_Score_Outcome3   Age_Score_Outcome2
#         0.8429                   118.13                       81.29              54.43375
# Smoking_Score_Outcome2 Hypertension_Score_Outcome2    CEA_Score_Outcome2     GFR_Score_Outcome2
#              0                 17.09091                     2.63517              39.37288
# NEU_Score_Outcome2         Age_Score_Outcome3       Smoking_Score_Outcome3
#       4.600887                 57.14286                           0
# Hypertension_Score_Outcome3 CEA_Score_Outcome3         GFR_Score_Outcome3    NEU_Score_Outcome3
#       13.43346                 3.487084                    2.114767              5.111798

example_patient <- list(
  Age = 55,
  Smoking = 0,
  Hypertension = 0,
  CEA = 0.57,
  GFR = 78.82,
  NEU = 2.69
  )
# Outcome_1_Probability     Outcome_2_Probability     Outcome_3_Probability  Linear_Predictor_2vs1
#          0.534                    0.198                      0.2679               -0.9919
# Linear_Predictor_3vs1    Nomogram_Score_Outcome2   Nomogram_Score_Outcome3   Age_Score_Outcome2
#        -0.6897                    83.64                       57.78              47.62953
# Smoking_Score_Outcome2 Hypertension_Score_Outcome2    CEA_Score_Outcome2     GFR_Score_Outcome2
#              0                        0                   0.8207907              30.00745
# NEU_Score_Outcome2         Age_Score_Outcome3       Smoking_Score_Outcome3
#        5.184169                  50.000                           0
# Hypertension_Score_Outcome3 CEA_Score_Outcome3         GFR_Score_Outcome3    NEU_Score_Outcome3
#              0                 1.086141                     3.12336               3.57162