Notes on the Derivation of the SFPS Weights
================

<!-- SFPSDerivationNotes.md is generated from SFPSDerivationNotes.Rmd.
Please edit SFPSDerivationNotes.Rmd instead. -->

This document records questions that arose while independently deriving
the SFPS weighting method proposed in:

> Ciardulli, S., Fontana, N., Vantini, S., & Ieva, F. (2026).
> *Generalized propensity score weighting for functional causal
> inference framework*. arXiv:2608.03200.

These notes are intended to document the derivation process and should
not be interpreted as definitive claims of errors in the paper. Some of
the points below may result from omitted intermediate steps, differences
in parameterization, or misunderstandings of the derivation.

## 1. Inner optimization over the weights

The paper formulates the weight estimation problem using a Lagrangian
and then considers a dual problem of the form

$$\max_{\gamma, \lambda} \inf_w \mathcal{L}(w, \gamma, \lambda).$$

For the terms involving an individual weight $w_i$, the displayed
Lagrangian leads to an expression of the form

$$\ell_i(w_i) = \log w_i + w_i(-n+\gamma^\top g_i),$$

up to terms that do not depend on $w_i$.

The first-order condition is

$$\frac{\partial \ell_i}{\partial w_i} = \frac{1}{w_i} - n + \gamma^\top g_i = 0,$$

which gives

$$w_i = \frac{1}{n-\gamma^\top g_i}.$$

However,

$$\frac{\partial^2 \ell_i}{\partial w_i^2} = - \frac{1}{w_i^2} < 0.$$

Therefore, the stationary point is a maximum of $\ell_i(w_i)$ rather
than a minimum.

Moreover, if $w_i > 0$ and values arbitrarily close to zero are allowed,

$$\log w_i \rightarrow -\infty \qquad \text{as} \qquad w_i \rightarrow 0^+.$$

Consequently, the inner infimum appears to be unbounded below.

This raises the following question:

> How does the stationary solution above correspond to the inner
> $\inf_w$ in the stated dual problem?

A sign change in the primal objective, a different domain for the
weights, or an additional step in the dual derivation might resolve this
issue.

## 2. Does standardization imply $\bar{g}=0$?

The balancing vector contains the standardized FPC scores, standardized
covariates, and their interaction terms. Schematically,

$$g_i = 
\begin{bmatrix}
  A_i \\
  C_i \\
  \mathrm{vec}(A_i C_i^\top)
\end{bmatrix}.$$

Because $A$ and $C$ are centered during standardization,

$$\frac{1}{n}\sum_{i=1}^n A_i=0,
\qquad
\frac{1}{n}\sum_{i=1}^n C_i=0.$$

These conditions, however, do not generally imply

$$\frac{1}{n} \sum_{i=1}^n A_i C_i^\top = 0.$$

The latter quantity is essentially the sample cross-covariance between
the FPC scores and the covariates and is generally nonzero when the
treatment and covariates are associated.

Therefore,

$$\overline{g} = \frac{1}{n}\sum_{i=1}^n g_i$$

does not appear to be zero merely because $A$ and $C$ have been
standardized.

This point matters if $\overline{g} = 0$ is used to simplify terms in
the dual objective.

The following numerical example illustrates this point:

``` r
library(FPScausal)
library(fdapace)

n <- 200

data <- simulate_fps_data(n, include_functional_cov = FALSE, seed = 314)
time_grid <- seq(0, 1, length.out = 51)

Ly <- data$X |> split(seq_len(n))
Lt <- rep(list(time_grid), times = n)
fpca <- FPCA(Ly, Lt)

A <- t(scale(fpca$xiEst)) # (6, n)
C <- t(scale(data$C))     # (3, n)

vec_AC <- matrix(nrow = 6*3, ncol = n)
for (i in seq_len(n)) {
  vec_AC[, i] <- as.vector(A[, i, drop = FALSE] %*% t(C[, i]))
}

g <- rbind(A, C, vec_AC)
g_bar <- rowMeans(g)

round(g_bar, digits = 16)
```

    ##  [1]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
    ##  [6]  0.000000000  0.000000000  0.000000000  0.000000000 -0.723664508
    ## [11]  0.158523694  0.022831518 -0.040351459 -0.049843781 -0.098803154
    ## [16] -0.150461667 -0.336578227  0.009993705 -0.089925354  0.108218632
    ## [21]  0.201455570 -0.034099848  0.171293932 -0.450929043  0.038681773
    ## [26]  0.086003859 -0.113051780

## Status

These issues are currently under investigation.
