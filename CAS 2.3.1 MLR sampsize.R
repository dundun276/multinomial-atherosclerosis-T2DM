# ==============================================================================
# Sample size calculation for the development of a multinomial prediction model
# ==============================================================================
# Methodology framework:
#   Primary reference: Pate et al. Stat Methods Med Res. 2023;32:555-571. doi:10.1177/09622802231151220
#   Code from: Pate A. MRC-multi-outcome GitHub repository: 
#              https://github.com/alexpate30/MRC-multi-outcome
#              (File: p3.1 worked example.R)
#
# Key methodological adaption:
#   Replaced the need for prior c-statistics (original Pate approach) with conservative R2 approximation for each pairwise submodel
#   Gehringer et al. J Clin Epidemiol. 2024;174:111481. doi:10.1016/j.jclinepi.2024.111481
#   GitHub: https://github.com/celinagehringer/multinomial_risk_prediction_RA
#           (File: Section3.2.2_MLR_sampsize.R)
#
# Adapted by: Yongqi Zheng | Contact: dundun276@126.com
# Usage: For academic non-commercial research only. Citation required.
#
# Study-specific parameters:
#   Outcome categories: Normal (n = 467), Thickening (n = 102), Plaque (n = 554)
#   Predictors: Q = 6
#   Conservative assumption: Overall Cox-Snell R2 = 0.15 (due to lack of prior c-statistics)
# ==============================================================================
### STEP'S 1 and 2: Identify values for Q, p_k, p_k_r, max(R2_CS), R2_CS_adj and 
## R2_CS_adj_k_r
# Define the number of events in each category
EV1 <- 467 # Normal group
EV2 <- 102 # Thickening group
EV3 <- 554 # Plaque group

# Define K - number of outcome categories
K <- 3

# Define Q - number of predictor parameters
Q <- 6

## Define p_k and p_k_r
# Define total number of events
n.total <- EV1 + EV2 + EV3 

# Calculate p_k
p.1 <- EV1/n.total
p.2 <- EV2/n.total
p.3 <- EV3/n.total
p.1 # 0.4158504
p.2 # 0.09082814
p.3 # 0.4933215

# Calculate p_k_r
p.1.2 <- (EV1 + EV2)/n.total
p.1.3 <- (EV1 + EV3)/n.total
p.2.3 <- (EV2 + EV3)/n.total

## Calculate max(R2_CS)
max_R2_CS <- 1 - ((p.1^p.1) * (p.2^p.2) * (p.3^p.3))^2
max_R2_CS # 0.8447424

## Calculate R2_CS_adj
# Calculate an estimte of R2_CS_app, based off R2_NAGEL = 0.15
R2_CS_adj <- 0.15 * max_R2_CS
R2_CS_adj # 0.1267114

## Calculate R2_CS_adj_k_r
## Since no prior information on c-statistics is available,  we will use the conservative approach of R2 = 0.15
# Calculate pairwise outcome proportions (phi), of category k relative to category i
phi.1.2 <- EV2/(EV1 + EV2)
phi.1.3 <- EV3/(EV1 + EV3)
phi.2.3 <- EV3/(EV2 + EV3)
phi.1.2 # 0.1792619
phi.1.3 # 0.5426053
phi.2.3 # 0.8445122

R2_app_kr.1.2 <- 0.15 * (1 - ((phi.1.2^(phi.1.2)) * ((1-phi.1.2)^(1-phi.1.2)))^2)
R2_app_kr.1.3 <- 0.15 * (1 - ((phi.1.3^(phi.1.3)) * ((1-phi.1.3)^(1-phi.1.3)))^2)
R2_app_kr.2.3 <- 0.15 * (1 - ((phi.2.3^(phi.2.3)) * ((1-phi.2.3)^(1-phi.2.3)))^2)
R2_app_kr.1.2 # 0.09143773
R2_app_kr.1.3 # 0.1122264
R2_app_kr.2.3 # 0.08679315



### STEP 3: Criterion (i), ensure minimal overfitting (shrinkage factor S = 0.9)
## Let S be the level of shrinkage we are targeting
S <- 0.9

## Calculate m_k_r
m.1.2 <- Q/((S - 1) * log(1 - R2_app_kr.1.2/S))
m.1.3 <- Q/((S - 1) * log(1 - R2_app_kr.1.3/S))
m.2.3 <- Q/((S - 1) * log(1 - R2_app_kr.2.3/S))

## Calculate n_k_r for criterion (i) for each submodel
N_C1.1.2 <- m.1.2/p.1.2
N_C1.1.3 <- m.1.3/p.1.3
N_C1.2.3 <- m.2.3/p.2.3
N_C1.1.2 # 1105.297
N_C1.1.3 # 495.5108
N_C1.2.3 # 1012.86

## Take the ceiling of the maximum of these as the sample size for criteiron (i)
N_C1 <- ceiling(max(N_C1.1.2, N_C1.1.3, N_C1.2.3))
N_C1 # 1106

## Now calculate number of each event we expect to see in a datast of this size
N_C1 * p.1 # 459.9305
N_C1 * p.2 # 100.4559
N_C1 * p.3 # 545.6135



### STEP 4: Criterion (ii), ensure precise estimate of overall model performance (within 5% of max R2)
N_C2 <- 4 * Q/((R2_CS_adj/(R2_CS_adj + 0.05 * max_R2_CS) - 1) * log(1 - R2_CS_adj - 0.05 * max_R2_CS))
N_C2 <- ceiling(N_C2)
N_C2 # 519

# Now calculate number of each event we expect to see in a datast of this size
N_C2 * p.1 # 215.8264
N_C2 * p.2 # 47.1398
N_C2 * p.3 # 256.0338



### STEP 5: Criterion (iii), ensure precise estimate of outcome-specific risks (margin of error < 0.05)
N_C3.1 <- qchisq(1 - 0.05/K, 1) * p.1 * (1 - p.1)/0.05^2 
N_C3.2 <- qchisq(1 - 0.05/K, 1) * p.2 * (1 - p.2)/0.05^2 
N_C3.3 <- qchisq(1 - 0.05/K, 1) * p.3 * (1 - p.3)/0.05^2 
N_C3.1 # 556.8807
N_C3.2 # 189.3073
N_C3.3 # 573.0117

N_C3 <- ceiling(max(N_C3.1, N_C3.2, N_C3.3))
N_C3 # 574

## Now calculate number of each event we expect to see in a dataset of this size
N_C3 * p.1 # 238.6981
N_C3 * p.2 # 52.13535
N_C3 * p.3 # 283.1665



### STEP 6: Take the maximum sample size across all three criteria
N_C1  # Criterion (i): Minimize overfitting (shrinkage approach)
N_C2  # Criterion (ii): Precise overall model performance (R2 precision)
N_C3  # Criterion (iii): Precise outcome-specific risk estimates (margin of error)
N <- max(N_C1, N_C2, N_C3)
N # 1106

### Minimum sample size: n = 1106