source("Scripts/00_packages_functions.R")

TxSmry <- function(x) {
  factor(x,
         levels = c("Treatment combination",
                    "None", 
                    "Surg only", 
                    "SACT only",
                    "RTx only")
  )}

cohort6 <- readRDS("Scratch_data/cohort6.Rds")
cohort6 <- cohort6 %>% 
  mutate(tx_smry = TxSmry(tx_smry),
         any_tx = if_else(tx_smry == "None", 0, 1))
#cohort6 <- cohort6 %>% filter(sex_name == "Male") # when checking to see if outputs were similar in males and females
#cohort6 <- cohort6 %>% filter(sex_name == "Female") # when checking to see if outputs were similar in males and females

## Quick method to test point estimates only - 1 year survival ---- 

# Fit the mediation model
m_mod <- glm(
  any_tx ~ egfr + age + sex + cancer_type + mm_count + n_cancers + m_stage,
  data = cohort6,
  family = binomial
)

# Fit the outcome model
y_mod <- glm(
  surv1y ~ ns(egfr, df = 5) * any_tx + age + sex + cancer_type + mm_count + n_cancers + m_stage,
  data = cohort6, 
  family = binomial
)

# Set exposure values
a0 <- 20 # Low value
a1 <- 70 # High value - selected as this is the median and close to mean value and consistent with lowest adjusted likelihood of "no treatment" across cancer types

# Compute counterfactual mediator values
cohort6$m_a0 <- predict(m_mod, 
                        newdata = transform(cohort6, egfr = a0), 
                        type = "response")
cohort6$m_a1 <- predict(m_mod, 
                        newdata = transform(cohort6, egfr = a1), 
                        type = "response")

# Helper function: predic probability for any A, M

pred_y <- function(A,M) {
  nd <- cohort6 
  nd$egfr <- A
  nd$any_tx <- M
  predict(y_mod, newdata = nd, type = "response")
}

# Four required potential outcomes
Y_a1_m_a1 <- pred_y(a1, cohort6$m_a1)
Y_a1_m_a0 <- pred_y(a1, cohort6$m_a0)
Y_a0_m_a1 <- pred_y(a0, cohort6$m_a1)
Y_a0_m_a0 <- pred_y(a0, cohort6$m_a0)
mean(Y_a1_m_a0 - Y_a0_m_a0)

# Controlled direct effect (CDE)
CDE <- mean(Y_a1_m_a0 - Y_a0_m_a0, na.rm = TRUE)
# Pure indirect effect (PIE)
PIE <- mean(Y_a0_m_a1 - Y_a0_m_a0, na.rm = TRUE)
# Reference interaction (INTref)
INTref <- mean(Y_a1_m_a0 - Y_a0_m_a0, na.rm = TRUE) - CDE - PIE
# Mediated interaction (INTmed)
INTmed <- mean(Y_a1_m_a1 - Y_a1_m_a0, na.rm = TRUE) - PIE
# Total effect
TE <- mean(Y_a1_m_a1 - Y_a0_m_a0, na.rm = TRUE)

## Output
decomp <- list(
  CDE = CDE, 
  PIE = PIE,
  INTref = INTref,
  INTmed = INTmed,
  Total_Effect = TE
)
decomp


## Estimates and confidence intervals for 4-way decomposition - 1 YEAR SURVIVAL ----

decomp_fun1 <- function(data, indices, a0 = 20, a1 = 70, df_spline = 5) {
  
  d <- data[indices, ]
  
  # Fit the mediation model
  m_mod <- glm(
    any_tx ~ egfr + age + sex + cancer_type + mm_count + n_cancers + m_stage,
    data = d,
    family = binomial
  )
  
  # Fit the outcome model with splines
  y_mod <- glm(
    surv1y ~ ns(egfr, df = df_spline) * any_tx + age + sex + cancer_type + mm_count + n_cancers + m_stage,
    data = d, 
    family = binomial
  )
  
  # Mediator probabilities under A = a0 and A = a1
  
  p_a0 <- predict(m_mod, newdata = transform(d, egfr = a0), type = "response")
  p_a1 <- predict(m_mod, newdata = transform(d, egfr = a1), type = "response")
  
  pred_y_bin <- function(A,M) {
    nd <- d
    nd$egfr <- A
    nd$any_tx <- M
    predict(y_mod, newdata = nd, type = "response")
  }
  
  Ey_a1_M1 <- pred_y_bin(a1, 1)
  Ey_a1_M0 <- pred_y_bin(a1, 0)
  Ey_a0_M1 <- pred_y_bin(a0, 1)
  Ey_a0_M0 <- pred_y_bin(a0, 0)
  
  # Integrate over M:
  Y_a1_m_a1 <- p_a1 * Ey_a1_M1 + (1 - p_a1) * Ey_a1_M0
  Y_a1_m_a0 <- p_a0 * Ey_a1_M1 + (1 - p_a0) * Ey_a1_M0
  Y_a0_m_a1 <- p_a1 * Ey_a0_M1 + (1 - p_a1) * Ey_a0_M0
  Y_a0_m_a0 <- p_a0 * Ey_a0_M1 + (1 - p_a0) * Ey_a0_M0
  
  
  # Components
  CDE <- mean(Y_a1_m_a0 - Y_a0_m_a0, na.rm = TRUE)
  PIE <- mean(Y_a0_m_a1 - Y_a0_m_a0, na.rm = TRUE)
  INTref <- mean(Y_a1_m_a0 - Y_a0_m_a0, na.rm = TRUE) - CDE - PIE
  INTmed <- mean(Y_a1_m_a1 - Y_a1_m_a0, na.rm = TRUE) - PIE
  TE <- mean(Y_a1_m_a1 - Y_a0_m_a0, na.rm = TRUE)
  
  return(c(CDE = CDE, PIE = PIE, INTref = INTref, INTmed = INTmed, TE = TE))
  
}

set.seed(123)
boot_res1 <- boot(
  data = cohort6,
  statistic = decomp_fun1, 
  R=1000
)

effect_names <- colnames(boot_res1$t)
if (is.null(effect_names) || any(effect_names == "")) {
  effect_names <- c("CDE", "PIE", "INTref", "INTmed", "TE")
}

ci_table_anytx_1y <- data.frame(
  effect = effect_names, 
  estimate = as.numeric(boot_res1$t0),
  ci_lower = apply(boot_res1$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  ci_upper = apply(boot_res1$t, 2, quantile, probs = 0.975, na.rm = TRUE)
)

ci_table_anytx_1y
ci_table_anytx_1y %>% write.csv("Outputs/decomposition_estimates_ci_anytx_1y.csv")

## Estimates and confidence intervals for 4-way decomposition - 5 YEAR SURVIVAL ----

decomp_fun2 <- function(data, indices, a0 = 20, a1 = 70, df_spline = 5) {
  
  d <- data[indices, ]
  
  # Fit the mediation model
  m_mod <- glm(
    any_tx ~ egfr + age + sex + cancer_type + mm_count + n_cancers + m_stage,
    data = d,
    family = binomial
  )
  
  # Fit the outcome model with splines
  y_mod <- glm(
    surv5y ~ ns(egfr, df = df_spline) * any_tx + age + sex + cancer_type + mm_count + n_cancers + m_stage,
    data = d, 
    family = binomial
  )
  
  # Mediator probabilities under A = a0 and A = a1
  
  p_a0 <- predict(m_mod, newdata = transform(d, egfr = a0), type = "response")
  p_a1 <- predict(m_mod, newdata = transform(d, egfr = a1), type = "response")
  
  pred_y_bin <- function(A,M) {
    nd <- d
    nd$egfr <- A
    nd$any_tx <- M
    predict(y_mod, newdata = nd, type = "response")
  }
  
  Ey_a1_M1 <- pred_y_bin(a1, 1)
  Ey_a1_M0 <- pred_y_bin(a1, 0)
  Ey_a0_M1 <- pred_y_bin(a0, 1)
  Ey_a0_M0 <- pred_y_bin(a0, 0)
  
  # Integrate over M:
  Y_a1_m_a1 <- p_a1 * Ey_a1_M1 + (1 - p_a1) * Ey_a1_M0
  Y_a1_m_a0 <- p_a0 * Ey_a1_M1 + (1 - p_a0) * Ey_a1_M0
  Y_a0_m_a1 <- p_a1 * Ey_a0_M1 + (1 - p_a1) * Ey_a0_M0
  Y_a0_m_a0 <- p_a0 * Ey_a0_M1 + (1 - p_a0) * Ey_a0_M0
  
  
  # Components
  CDE <- mean(Y_a1_m_a0 - Y_a0_m_a0, na.rm = TRUE)
  PIE <- mean(Y_a0_m_a1 - Y_a0_m_a0, na.rm = TRUE)
  INTref <- mean(Y_a1_m_a0 - Y_a0_m_a0, na.rm = TRUE) - CDE - PIE
  INTmed <- mean(Y_a1_m_a1 - Y_a1_m_a0, na.rm = TRUE) - PIE
  TE <- mean(Y_a1_m_a1 - Y_a0_m_a0, na.rm = TRUE)
  
  return(c(CDE = CDE, PIE = PIE, INTref = INTref, INTmed = INTmed, TE = TE))
  
}

set.seed(123)
boot_res2 <- boot(
  data = cohort6,
  statistic = decomp_fun2, 
  R=1000
)

effect_names <- colnames(boot_res2$t)
if (is.null(effect_names) || any(effect_names == "")) {
  effect_names <- c("CDE", "PIE", "INTref", "INTmed", "TE")
}

ci_table_anytx_5y <- data.frame(
  effect = effect_names, 
  estimate = as.numeric(boot_res2$t0),
  ci_lower = apply(boot_res2$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  ci_upper = apply(boot_res2$t, 2, quantile, probs = 0.975, na.rm = TRUE)
)

ci_table_anytx_5y
ci_table_anytx_5y %>% write.csv("Outputs/decomposition_estimates_ci_anytx_5y.csv")
