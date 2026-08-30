Functional Propensity Score Study
================
Koji Makiyama

<!-- README.md is generated from README.Rmd Please edit that file -->

This repository contains an independent study and reproduction of the
methods proposed in:

> Ciardulli, S., Fontana, N., Vantini, S., & Ieva, F. (2026).
> *Generalized propensity score weighting for functional causal
> inference framework*. arXiv:2608.03200.

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

## Further investigations

In addition to reproducing the reported simulation results, this
repository may be used to investigate several aspects of the method,
including:

- the role of FPCA truncation in SFPS estimation;
- the balancing conditions used to estimate the SFPS weights;
- the dual optimization problem used in the proposed method;
- comparison with the method of Zhang, Xue, and Wang (2021); and
- comparison with multivariate generalized propensity score methods.

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
