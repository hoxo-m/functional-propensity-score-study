
``` r
library(tidyverse)

df_setting <- tribble(
  ~ `Treat-Cov`, ~ `Cov-Outcome`, ~ `PVE`, ~ `PVE*`, ~ `Weight`, 
  "Linear",     "Linear",     "",     "0.95", "Unweighted", 
  "Linear",     "Linear",     "",     "0.99", "Unweighted", 
  "Linear",     "Linear",     "0.95", "0.95", "Weighted",
  "Linear",     "Linear",     "0.95", "0.99", "Weighted",
  "Linear",     "Linear",     "0.99", "0.95", "Weighted",
  "Linear",     "Linear",     "0.99", "0.99", "Weighted",
  "Non-linear", "Linear",     "",     "0.95", "Unweighted", 
  "Non-linear", "Linear",     "",     "0.99", "Unweighted", 
  "Non-linear", "Linear",     "0.95", "0.95", "Weighted",
  "Non-linear", "Linear",     "0.95", "0.99", "Weighted",
  "Non-linear", "Linear",     "0.99", "0.95", "Weighted",
  "Non-linear", "Linear",     "0.99", "0.99", "Weighted",
  "Linear",     "Non-linear", "",     "0.95", "Unweighted", 
  "Linear",     "Non-linear", "",     "0.99", "Unweighted", 
  "Linear",     "Non-linear", "0.95", "0.95", "Weighted",
  "Linear",     "Non-linear", "0.95", "0.99", "Weighted",
  "Linear",     "Non-linear", "0.99", "0.95", "Weighted",
  "Linear",     "Non-linear", "0.99", "0.99", "Weighted",
  "Non-linear", "Non-linear", "",     "0.95", "Unweighted", 
  "Non-linear", "Non-linear", "",     "0.99", "Unweighted", 
  "Non-linear", "Non-linear", "0.95", "0.95", "Weighted",
  "Non-linear", "Non-linear", "0.95", "0.99", "Weighted",
  "Non-linear", "Non-linear", "0.99", "0.95", "Weighted",
  "Non-linear", "Non-linear", "0.99", "0.99", "Weighted",
)
```

``` r
box::use(./R/data_generation)
box::use(estimate = ./R/estimate_causal_effect)

n_iter <- 200 # Number of iteration
n_size <- 200 # Data size

dummy_data <- data_generation$simulate_data()
true_effect_coefficients <- c(2, 1, 0.5, 0.5, 0, 0)
mu_true <- drop(t(dummy_data$true_eigenfunction_phi_values) %*% 
                  true_effect_coefficients)

n_time = length(dummy_data$time_grid)

dx <- diff(dummy_data$time_grid)
my_integrate <- function(x) {
  sum(dx * (x[-1] + x[-length(x)]) / 2)
}

execute_experiment <- function(settings) {
  mu_hat_matrix <- matrix(nrow = n_iter, ncol = n_time)
  ise_vec <- double(n_iter)
  for (j in seq_len(n_iter)) {
    data <- data_generation$simulate_data(
      n = n_size, 
      treatment_covariate_relation = settings$`Treat-Cov`, 
      covariate_outcome_relation = settings$`Cov-Outcome`)
    
    pve <- settings$PVE
    pve_star <- settings$`PVE*`
    use_weight <- settings$Weight == "Weighted"
    
    mu_hat <- estimate$estimate_causal_effect(
      data$outcome_Y, data$functional_treatment_X, data$time_grid,
      data$covariates_C, pve, pve_star, use_weight)
    
    mu_hat_matrix[j, ] <- mu_hat
    
    squared_error <- (mu_hat - mu_true)^2
    ise_vec[j] <- my_integrate(squared_error)
  }
  
  isb <- my_integrate((colMeans(mu_hat_matrix) - mu_true)^2)
  data.frame(MISE = median(ise_vec), AISE = mean(ise_vec), ISB = isb)
}
```

``` r
library(future)
library(furrr)

workers <- min(availableCores(), nrow(df_setting) / 2)

plan(multisession, workers = workers)

df_setting1 <- df_setting |> filter(`PVE*` == 0.95)

list_settings <- df_setting1 |> split(seq_len(nrow(df_setting1)))
result <- future_map_dfr(list_settings, execute_experiment,
                         .options = furrr_options(seed = 314))

df1 <- cbind(df_setting1, round(result, digits = 4))
```

``` r
knitr::kable(df1)
```

| Treat-Cov  | Cov-Outcome | PVE  | PVE\* | Weight     |   MISE |   AISE |    ISB |
|:-----------|:------------|:-----|:------|:-----------|-------:|-------:|-------:|
| Linear     | Linear      |      | 0.95  | Unweighted | 0.3124 | 0.3258 | 0.2428 |
| Linear     | Linear      | 0.95 | 0.95  | Weighted   | 0.1679 | 0.2008 | 0.0006 |
| Linear     | Linear      | 0.99 | 0.95  | Weighted   | 0.1897 | 0.2566 | 0.0040 |
| Non-linear | Linear      |      | 0.95  | Unweighted | 0.3497 | 0.3759 | 0.2627 |
| Non-linear | Linear      | 0.95 | 0.95  | Weighted   | 0.1100 | 0.1382 | 0.0001 |
| Non-linear | Linear      | 0.99 | 0.95  | Weighted   | 0.1197 | 0.1537 | 0.0003 |
| Linear     | Non-linear  |      | 0.95  | Unweighted | 0.3066 | 0.3270 | 0.2474 |
| Linear     | Non-linear  | 0.95 | 0.95  | Weighted   | 0.1816 | 0.2234 | 0.0011 |
| Linear     | Non-linear  | 0.99 | 0.95  | Weighted   | 0.1980 | 0.2494 | 0.0010 |
| Non-linear | Non-linear  |      | 0.95  | Unweighted | 0.3422 | 0.3604 | 0.2454 |
| Non-linear | Non-linear  | 0.95 | 0.95  | Weighted   | 0.1029 | 0.1345 | 0.0004 |
| Non-linear | Non-linear  | 0.99 | 0.95  | Weighted   | 0.1190 | 0.1527 | 0.0005 |

``` r
df_setting2 <- df_setting |> filter(`PVE*` == 0.99)

list_settings <- df_setting2 |> split(seq_len(nrow(df_setting2)))
result <- future_map_dfr(list_settings, execute_experiment,
                         .options = furrr_options(seed = 314))

df2 <- cbind(df_setting2, round(result, digits = 4))
```

``` r
knitr::kable(df2)
```

| Treat-Cov  | Cov-Outcome | PVE  | PVE\* | Weight     |   MISE |   AISE |    ISB |
|:-----------|:------------|:-----|:------|:-----------|-------:|-------:|-------:|
| Linear     | Linear      |      | 0.99  | Unweighted | 0.6247 | 0.8058 | 0.2454 |
| Linear     | Linear      | 0.95 | 0.99  | Weighted   | 0.8945 | 1.3330 | 0.0079 |
| Linear     | Linear      | 0.99 | 0.99  | Weighted   | 0.9519 | 1.3348 | 0.0173 |
| Non-linear | Linear      |      | 0.99  | Unweighted | 0.7626 | 0.8904 | 0.2644 |
| Non-linear | Linear      | 0.95 | 0.99  | Weighted   | 0.6879 | 0.9840 | 0.0012 |
| Non-linear | Linear      | 0.99 | 0.99  | Weighted   | 0.6179 | 0.8772 | 0.0024 |
| Linear     | Non-linear  |      | 0.99  | Unweighted | 0.6382 | 0.7723 | 0.2500 |
| Linear     | Non-linear  | 0.95 | 0.99  | Weighted   | 1.2044 | 1.4140 | 0.0087 |
| Linear     | Non-linear  | 0.99 | 0.99  | Weighted   | 0.8308 | 1.2589 | 0.0042 |
| Non-linear | Non-linear  |      | 0.99  | Unweighted | 0.7292 | 0.9061 | 0.2461 |
| Non-linear | Non-linear  | 0.95 | 0.99  | Weighted   | 0.7991 | 1.1496 | 0.0026 |
| Non-linear | Non-linear  | 0.99 | 0.99  | Weighted   | 0.6261 | 0.8800 | 0.0083 |
