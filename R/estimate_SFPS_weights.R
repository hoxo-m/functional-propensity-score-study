box::use(stats[optim])
box::use(fdapace[FPCA])

#' @export
functional_PCA <- function(functional_treatment_X, time_grid, pve_threshold) {
  n <- nrow(functional_treatment_X)
  
  Ly <- list()
  Lt <- list()
  for (i in seq_len(n)) {
    Ly[[i]] <- functional_treatment_X[i, ]
    Lt[[i]] <- time_grid
  }
  
  FPCA(Ly, Lt, optns = list(FVEthreshold = pve_threshold))
}

#' @export
estimate_SFPS_weights <- function(
    functional_treatment_X, time_grid, covariates_C, pve_threshold = 0.95) {
  
  n <- nrow(functional_treatment_X)

  fpca <- functional_PCA(functional_treatment_X, time_grid, pve_threshold)
  
  truncate_level <- fpca$selectK
  fpc_scores_A <- fpca$xiEst

  standardized_fpc_scores_A <- scale(fpc_scores_A) # (n,truncate_level)
  standardized_covariates_C <- scale(covariates_C) # (n,3)

  L <- truncate_level
  p <- ncol(standardized_covariates_C)
  
  interaction_moments_vec_AC <- matrix(nrow = n, ncol = L * p)
  for (i in seq_len(n)) {
    tAi <- t(standardized_fpc_scores_A[i, , drop = FALSE])
    Ci <- standardized_covariates_C[i, , drop = FALSE]
    temp <- tAi %*% Ci
    interaction_moments_vec_AC[i, ] <- as.vector(temp)
  }
  
  stacked_balancing_moments_g <- cbind(
    standardized_fpc_scores_A, 
    standardized_covariates_C, 
    interaction_moments_vec_AC
  )
  
  dual_objective <- function(theta) {
    linear_predictor <- drop(-stacked_balancing_moments_g %*% theta)
    m <- max(linear_predictor)
    log(sum(exp(linear_predictor - m))) + m
  }
  
  theta <- double(ncol(stacked_balancing_moments_g))
  opt <- optim(theta, dual_objective, method = "BFGS")
  theta_opt <- opt$par
  
  if (opt$convergence != 0) {
    warning("SFPS weight optimization did not converge.")
  }
  
  linear_predictor <- drop(-stacked_balancing_moments_g %*% theta_opt)
  linear_predictor <- linear_predictor - max(linear_predictor)
  
  weights <- exp(linear_predictor)
  weights / sum(weights)
}
