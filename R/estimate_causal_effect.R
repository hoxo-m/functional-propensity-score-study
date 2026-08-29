box::use(weight = ./estimate_SFPS_weights)
box::use(stats[lm, coef])

estimate_causal_effect <- function(
    outcome_Y, functional_treatment_X, time_grid, covariates_C, 
    pve_threshold, pve_threshold_star, use_weight = TRUE) {
  
  if (use_weight) {
    pve_threshold <- as.double(pve_threshold)
    weights <- weight$estimate_SFPS_weights(functional_treatment_X, time_grid, 
                                            covariates_C, pve_threshold)
  } else {
    weights <- NULL
  }

  pve_threshold_star <- as.double(pve_threshold_star)
  fpca <- weight$functional_PCA(
    functional_treatment_X, time_grid, pve_threshold_star)
  fpc_scores_A <- fpca$xiEst
  
  regression_data <- data.frame(outcome = outcome_Y, fpc_scores_A)
  fit <- lm(outcome ~ ., data = regression_data, weights = weights)

  effect_coefficients <- coef(fit)[-1]
  
  drop(fpca$phi %*% effect_coefficients)
}
