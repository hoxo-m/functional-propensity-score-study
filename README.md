Functional Propensity Score Study
================
Koji Makiyama

<!-- README.md is generated from README.Rmd Please edit that file -->

This repository contains an independent study and reproduction of the
methods proposed in:

> Ciardulli, S., Fontana, N., Vantini, S., & Ieva, F. (2026).
> *Generalized propensity score weighting for functional causal
> inference framework*. arXiv:2608.03200.

An implementation by the authors is available as the
[`FPScausal`](https://CRAN.R-project.org/package=FPScausal) R package on
CRAN.

The code in this repository is an independent reimplementation created
for methodological study and reproducibility.

The main goals of this repository are to:

- understand the proposed functional propensity score weighting method;
- independently implement the simulation design and SFPS estimator;
- reproduce the numerical results reported in the paper; and
- investigate the behavior and assumptions of the proposed method.

This is a study and reproduction repository, not an official
implementation of the method.

## Current status

The simulation study reported in Table 1 is currently implemented.

- [Reproduction of Table 1](Table1.md)

The implementation includes:

- generation of the functional treatment, covariates, and outcome under
  the four simulation settings;
- functional principal component analysis (FPCA);
- estimation of the proposed SFPS weights;
- weighted estimation of the causal effect function; and
- evaluation using MISE, AISE, and ISB.

## Implementation

The code in this repository is implemented independently of the authors’
`FPScausal` package. This allows the individual steps of the method,
including FPCA, construction of the balancing moments, SFPS weight
estimation, and causal effect estimation, to be examined separately.

The main R code is contained in the `R/` directory:

- `R/data_generation.R` implements the simulation data-generating
  process;
- `R/estimate_SFPS_weights.R` implements FPCA and estimation of the SFPS
  weights; and
- `R/estimate_causal_effect.R` estimates the causal effect function
  using weighted functional principal component scores.

The implementation is written independently from the official
`FPScausal` package so that the individual steps of the proposed method
can be examined and validated separately.

### Implementation notes

There is an ambiguity in the parameterization of the normal
distributions used in the paper.

For the outcome error, the paper specifies

$$\epsilon_i \sim N(0, 25).$$

We interpret 25 as the variance, so the R implementation uses
`rnorm(..., sd = 5)`.

In Supplementary Material C, some covariate errors are specified as
$N(0, 0.5)$. For these terms, this reproduction uses `0.5` directly as
the standard deviation, i.e. `rnorm(..., sd = 0.5)`.

This was because reproducing the values in Table 1 was impossible when
aligning with either notational convention.

## Reproduction of Table 1

`Table1.Rmd` reproduces the simulation experiment in Table 1 of the
paper.

For each of the four combinations of treatment–covariate and
covariate–outcome relationships, 200 data sets with a sample size of 200
are generated. The proposed SFPS-weighted estimator is compared with the
unweighted estimator.

The truncation levels used for weighting and outcome estimation are
selected separately using the PVE thresholds denoted by `PVE` and
`PVE*`, respectively.

See [Table1.md](Table1.md) for the full simulation setup, code, and
results.

## Requirements

The analyses are written in R. The main packages currently used include:

- `fdapace`
- `box`
- `tidyverse`
- `future`
- `furrr`
- `memoise`

## Notes

This repository is intended for methodological study and
reproducibility. Results may differ slightly from those reported in the
original paper because the implementation is independent and may use
different numerical routines.
