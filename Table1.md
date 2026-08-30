
``` r
library(tidyverse)

simulation_settings <- tribble(
  ~ `Treat-Cov`, ~ `Cov-Outcome`, ~ `PVE`, ~ `PVE*`, ~ `Weight`, 
    "Linear",      "Linear",        "",      "0.95",   "Unweighted", 
    "Linear",      "Linear",        "",      "0.99",   "Unweighted", 
    "Linear",      "Linear",        "0.95",  "0.95",   "Weighted",
    "Linear",      "Linear",        "0.95",  "0.99",   "Weighted",
    "Linear",      "Linear",        "0.99",  "0.95",   "Weighted",
    "Linear",      "Linear",        "0.99",  "0.99",   "Weighted",
    "Non-linear",  "Linear",        "",      "0.95",   "Unweighted", 
    "Non-linear",  "Linear",        "",      "0.99",   "Unweighted", 
    "Non-linear",  "Linear",        "0.95",  "0.95",   "Weighted",
    "Non-linear",  "Linear",        "0.95",  "0.99",   "Weighted",
    "Non-linear",  "Linear",        "0.99",  "0.95",   "Weighted",
    "Non-linear",  "Linear",        "0.99",  "0.99",   "Weighted",
    "Linear",      "Non-linear",    "",      "0.95",   "Unweighted", 
    "Linear",      "Non-linear",    "",      "0.99",   "Unweighted", 
    "Linear",      "Non-linear",    "0.95",  "0.95",   "Weighted",
    "Linear",      "Non-linear",    "0.95",  "0.99",   "Weighted",
    "Linear",      "Non-linear",    "0.99",  "0.95",   "Weighted",
    "Linear",      "Non-linear",    "0.99",  "0.99",   "Weighted",
    "Non-linear",  "Non-linear",    "",      "0.95",   "Unweighted", 
    "Non-linear",  "Non-linear",    "",      "0.99",   "Unweighted", 
    "Non-linear",  "Non-linear",    "0.95",  "0.95",   "Weighted",
    "Non-linear",  "Non-linear",    "0.95",  "0.99",   "Weighted",
    "Non-linear",  "Non-linear",    "0.99",  "0.95",   "Weighted",
    "Non-linear",  "Non-linear",    "0.99",  "0.99",   "Weighted",
)
```

``` r
box::use(./R/data_generation)
box::use(estimate = ./R/estimate_causal_effect)

n_replications <- 200
sample_size <- 200

get_data_impl <- function(treatment_covariate_relation, 
                          covariate_outcome_relation) {
  map(
    seq_len(n_replications), 
    \(dummy) data_generation$simulate_data(
      sample_size, treatment_covariate_relation, covariate_outcome_relation)
  )
}
get_data <- memoise::memoise(get_data_impl)

set.seed(314)
list_settings <- simulation_settings |> split(seq_len(nrow(simulation_settings)))
list_data <- list_settings |> 
  map(\(settings) get_data(settings$`Treat-Cov`, settings$`Cov-Outcome`))

dummy_data <- list_data[[1]][[1]]
true_effects_mu <- drop(t(dummy_data$true_eigenfunction_phi_values) %*% 
                          data_generation$true_effect_coefficients)

n_time_points <- length(dummy_data$time_grid)

dx <- diff(dummy_data$time_grid)
my_integrate <- function(x) {
  sum(dx * (x[-1] + x[-length(x)]) / 2)
}

run_simulation <- function(data_list, settings) {
  estimated_effects <- matrix(nrow = n_replications, ncol = n_time_points)
  ise <- double(n_replications)
  for (j in seq_len(n_replications)) {
    data <- data_list[[j]]
    
    pve <- settings$PVE
    pve_star <- settings$`PVE*`
    use_weight <- settings$Weight == "Weighted"
    
    mu_hat <- estimate$estimate_causal_effect(
      data$outcome_Y, data$functional_treatment_X, data$time_grid,
      data$covariates_C, pve, pve_star, use_weight)
    
    estimated_effects[j, ] <- mu_hat
    
    squared_error <- (mu_hat - true_effects_mu)^2
    ise[j] <- my_integrate(squared_error)
  }
  
  isb <- my_integrate((colMeans(estimated_effects) - true_effects_mu)^2)
  tibble(MISE = median(ise), AISE = mean(ise), ISB = isb)
}
```

``` r
library(future)
library(furrr)

n_workers <- min(availableCores(), nrow(simulation_settings))
plan(multisession, workers = n_workers)

result <- future_map2_dfr(list_data, list_settings, run_simulation)

df <- cbind(simulation_settings, result)
```

``` r
knitr::kable(df |> filter(`PVE*` == "0.95"), digits = 4)
```

| Treat-Cov  | Cov-Outcome | PVE  | PVE\* | Weight     |   MISE |   AISE |    ISB |
|:-----------|:------------|:-----|:------|:-----------|-------:|-------:|-------:|
| Linear     | Linear      |      | 0.95  | Unweighted | 0.3327 | 0.3440 | 0.2521 |
| Linear     | Linear      | 0.95 | 0.95  | Weighted   | 0.1709 | 0.2091 | 0.0006 |
| Linear     | Linear      | 0.99 | 0.95  | Weighted   | 0.1881 | 0.2300 | 0.0009 |
| Non-linear | Linear      |      | 0.95  | Unweighted | 0.3426 | 0.3732 | 0.2595 |
| Non-linear | Linear      | 0.95 | 0.95  | Weighted   | 0.1056 | 0.1411 | 0.0002 |
| Non-linear | Linear      | 0.99 | 0.95  | Weighted   | 0.1077 | 0.1420 | 0.0003 |
| Linear     | Non-linear  |      | 0.95  | Unweighted | 0.3103 | 0.3287 | 0.2496 |
| Linear     | Non-linear  | 0.95 | 0.95  | Weighted   | 0.1837 | 0.2142 | 0.0013 |
| Linear     | Non-linear  | 0.99 | 0.95  | Weighted   | 0.2144 | 0.2405 | 0.0020 |
| Non-linear | Non-linear  |      | 0.95  | Unweighted | 0.3427 | 0.3669 | 0.2421 |
| Non-linear | Non-linear  | 0.95 | 0.95  | Weighted   | 0.1241 | 0.1482 | 0.0011 |
| Non-linear | Non-linear  | 0.99 | 0.95  | Weighted   | 0.1189 | 0.1500 | 0.0011 |

``` r
knitr::kable(df |> filter(`PVE*` == "0.99"), digits = 4)
```

| Treat-Cov  | Cov-Outcome | PVE  | PVE\* | Weight     |   MISE |   AISE |    ISB |
|:-----------|:------------|:-----|:------|:-----------|-------:|-------:|-------:|
| Linear     | Linear      |      | 0.99  | Unweighted | 0.6283 | 0.7387 | 0.2533 |
| Linear     | Linear      | 0.95 | 0.99  | Weighted   | 0.7530 | 1.0895 | 0.0014 |
| Linear     | Linear      | 0.99 | 0.99  | Weighted   | 0.7520 | 1.0877 | 0.0024 |
| Non-linear | Linear      |      | 0.99  | Unweighted | 0.7641 | 0.9073 | 0.2619 |
| Non-linear | Linear      | 0.95 | 0.99  | Weighted   | 0.7420 | 0.9954 | 0.0012 |
| Non-linear | Linear      | 0.99 | 0.99  | Weighted   | 0.5586 | 0.7222 | 0.0025 |
| Linear     | Non-linear  |      | 0.99  | Unweighted | 0.5776 | 0.6886 | 0.2511 |
| Linear     | Non-linear  | 0.95 | 0.99  | Weighted   | 0.9300 | 1.1843 | 0.0026 |
| Linear     | Non-linear  | 0.99 | 0.99  | Weighted   | 0.9576 | 1.1908 | 0.0065 |
| Non-linear | Non-linear  |      | 0.99  | Unweighted | 0.7192 | 0.8713 | 0.2483 |
| Non-linear | Non-linear  | 0.95 | 0.99  | Weighted   | 0.8674 | 1.1178 | 0.0130 |
| Non-linear | Non-linear  | 0.99 | 0.99  | Weighted   | 0.6020 | 0.8279 | 0.0080 |
