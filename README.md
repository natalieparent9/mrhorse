
# mrhorse

Provides the R code to implement the methods described in:

Grant, AJ and Burgess, S (2024). [A Bayesian approach to Mendelian
randomization using summary statistics in the univariable and
multivariable settings with correlated
pleiotropy](https://doi.org/10.1016/j.ajhg.2023.12.002). *American
Journal of Human Genetics*. 111(1):165-180. doi:
<https://doi.org/10.1016/j.ajhg.2023.12.002>

## Installation

You can install the development version of mrhorse from
[GitHub](https://github.com/) with:

``` r
library(devtools)
source_URL("https://github.com/aj-grant/mr_horse.R")
```

## Implementation

The function `mr_horse()` implements the univariable MR-Horse method.
The required input is a data frame with column headings:

- `betaY`: estimate of genetic association with the outcome
- `betaYse`: standard error of the estimate of genetic association with
  the outcome
- `betaX`: estimate of genetic association with the exposure
- `betaXse`: standard error of the estimate of genetic association with
  the exposure

Each row in the data frame represents a genetic variant.

Optional arguments are:

- `n.chains`: number of chains to run (default = 3)
- `variable.names`: vector of parameters to save in the MCMC output. The
  causal effect is “theta”, and will always be saved. Other relevant
  parameters include “alpha” (pleiotropic effects for each genetic
  variant) and “rho” (correlations between pleiotropic effects and
  genetic variant-exposure effects for each variant)
- `n.iter`: number of iterations (in addition to burn-in, default =
  10000)
- `n.burnin`: number of iterations for burn-in (default = 10000)
- `stan`: fit the model using rstan (default = FALSE, uses JAGS)
- `n.cores`: number of cores to use in parallel when running multiple
  chains (default = parallelly::availableCores())

Output from the `mr_horse()` function include:

- `$MR_Estimate`: a data frame with the causal effect estimate (which is
  the posterior mean), standard deviation (i.e., the posterior standard
  deviation), upper and lower bounds of the 95% credible interval, and
  the R-hat value
- `$MR_Coda`: full MCMC samples for all parameters in `variable.names`

JAGS plotting tools can be used with the MCMC output, for example using
`traceplot()` and `densplot()`.

The function `mvmr_horse` implements the multivariable MVMR-Horse
method. The required inputs are as above, but where `betaX` is replaced
by the K columns `betaX1`, `betaX2`, …, and `betaXse` is replaced by the
K columns `betaX1se`, `betaX2se`, … .

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mrhorse)
```

``` r
# Fit model with JAGS
MREx = mr_horse(data_ex)
```

``` r
# Check estimates
MREx$MR_Estimate
#>   Estimate    SD 2.5% quantile 97.5% quantile  Rhat
#> 1    0.097 0.018         0.064          0.132 1.001
```