box::use(stats[rnorm])

sqrt2 <- sqrt(2)
sqrt3 <- sqrt(3)

evaluate_eigenfunctions <- function(t) {
  phi1 <- sqrt2 * sin(2 * pi * 1 * t)
  phi2 <- sqrt2 * cos(2 * pi * 1 * t)
  phi3 <- sqrt2 * sin(2 * pi * 2 * t)
  phi4 <- sqrt2 * cos(2 * pi * 2 * t)
  phi5 <- sqrt2 * sin(2 * pi * 3 * t)
  phi6 <- sqrt2 * cos(2 * pi * 3 * t)
  
  rbind(phi1, phi2, phi3, phi4, phi5, phi6)
}

generate_covariates_C <- function(latent_variables_Z, treatment_covariate_relation) {
  n <- nrow(latent_variables_Z)
  
  covariates_C <- matrix(nrow = n, ncol = 3L)
  if (treatment_covariate_relation == "Linear") {
    covariates_C[, 1] <- latent_variables_Z[, 1]
  } else {
    covariates_C[, 1] <- (latent_variables_Z[, 1] + 0.5)^2
  }
  covariates_C[, 2] <- 0.2 * latent_variables_Z[, 2]
  covariates_C[, 3] <- 0.2 * latent_variables_Z[, 3]
  
  covariate_errors_W <- rnorm(n = n * 3L, mean = 0, sd = c(1, 0.5, 0.5))
  covariate_errors_W <- matrix(covariate_errors_W, nrow = n, ncol = 3L, byrow = TRUE)
  
  covariates_C + covariate_errors_W
}

#' @export
true_effect_coefficients <- c(2, 1, 0.5, 0.5, 0, 0)

generate_outcome_Y <- function(fpc_scores_A, covariate_effect_g_C) {
  n <- nrow(fpc_scores_A)
  
  integral_term <- fpc_scores_A %*% true_effect_coefficients # (n,1)=(n,6)x(6,1)
  integral_term <- drop(integral_term) # to n-vector
  
  # outcome_errors <- rnorm(n = 1L, mean = 0, sd = 25)
  outcome_errors <- rnorm(n, mean = 0, sd = 5)
  
  1 + integral_term + covariate_effect_g_C + outcome_errors
}

#' @export
simulate_data <- function(
    n = 200L, 
    treatment_covariate_relation = c("Linear", "Non-linear"),
    covariate_outcome_relation = c("Linear", "Non-linear")) {
  
  treatment_covariate_relation <- match.arg(treatment_covariate_relation)
  covariate_outcome_relation <- match.arg(covariate_outcome_relation)
  
  time_grid <- seq(0, 1, length.out = 51)
  # n_t <- length(time_grid)

  latent_variables_Z <- rnorm(n * 6L, mean = 0, sd = 1)
  latent_variables_Z <- matrix(latent_variables_Z, nrow = n, ncol = 6L) # (n,6)

  fpc_score_scales <- c(4, 2 * sqrt3, 2 * sqrt2, 2, 1, 1 / sqrt2)
  fpc_score_scales <- matrix(fpc_score_scales, nrow = n, ncol = 6L, byrow = TRUE)
  fpc_scores_A <- latent_variables_Z * fpc_score_scales # (n,6)

  eigenfunction_phi_values <- evaluate_eigenfunctions(time_grid) # (6,n_t)
  
  functional_treatment_X <- 
    fpc_scores_A %*% eigenfunction_phi_values # (n,n_t) = (n,6) x (6,n_t)

  covariates_C <- 
    generate_covariates_C(latent_variables_Z, treatment_covariate_relation) # (n,3)
  
  covariate_effect_g_C <- 
    if (covariate_outcome_relation == "Linear") {
      2 * covariates_C[, 1]
    } else {
      2 * covariates_C[, 1] + covariates_C[, 2]^2
    }

  outcome_Y <- generate_outcome_Y(fpc_scores_A, covariate_effect_g_C)
  
  list(time_grid = time_grid, functional_treatment_X = functional_treatment_X, 
       covariates_C = covariates_C, outcome_Y = outcome_Y, 
       true_fpc_scores_A = fpc_scores_A, 
       true_eigenfunction_phi_values = eigenfunction_phi_values)
}
