
# mrhorse

Provides the R code to implement the methods described in:

Grant, AJ and Burgess, S (2024). [A Bayesian approach to Mendelian
randomization using summary statistics in the univariable and
multivariable settings with correlated
pleiotropy](https://doi.org/10.1016/j.ajhg.2023.12.002). *American
Journal of Human Genetics*. 111(1):165-180. doi:
<https://doi.org/10.1016/j.ajhg.2023.12.002>

The paper introduces a method for performing both univariable and
multivariable Mendelian randomization to examine the causal effect of
one or more exposures on an outcome using genetic association data. The
methods take as input estimates of the associations between genetic
instruments and the exposure(s) and outcome from GWAS summary
statistics. Causal effect estimation is performed in a Bayesian
framework using a horseshoe shrinkage prior to account for both
correlated and uncorrelated pleiotropic effects.

## Installation

In order to run the code,
[JAGS](https://sourceforge.net/projects/mcmc-jags/) first needs to be
installed.

You can install the development version of mrhorse from
[GitHub](https://github.com/) with:

``` r
library(devtools)
devtools::install_github("natalieparent9/mrhorse")
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/6p/3_l0pbvj53lbsc774mttyysr0000gn/T/RtmpjwlRUF/remotes11152523bd813/natalieparent9-mrhorse-63191e6/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/6p/3_l0pbvj53lbsc774mttyysr0000gn/T/RtmpjwlRUF/remotes11152523bd813/natalieparent9-mrhorse-63191e6/DESCRIPTION’
#>   ─  preparing ‘mrhorse’:
#>      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  cleaning src
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>        NB: this package now depends on R (>= 3.5.0)
#>        WARNING: Added dependency on R >= 3.5.0 because serialized objects in
#>      serialize/load version 3 cannot be read in older versions of R.
#>      File(s) containing such objects:
#>        ‘mrhorse/README_cache/gfm/unnamed-chunk-3_fc0074dc90533e48417417587eb05da2.RData’
#>        ‘mrhorse/README_cache/gfm/unnamed-chunk-3_fc0074dc90533e48417417587eb05da2.rdx’
#>        ‘mrhorse/README_cache/gfm/unnamed-chunk-4_06618c0614294ee3282d94231388fd59.RData’
#>        ‘mrhorse/README_cache/gfm/unnamed-chunk-4_06618c0614294ee3282d94231388fd59.rdx’
#>        ‘mrhorse/README_cache/gfm/unnamed-chunk-5_cec2b9a5f79319fa17e2a46445c4d61e.RData’
#>        ‘mrhorse/README_cache/gfm/unnamed-chunk-5_cec2b9a5f79319fa17e2a46445c4d61e.rdx’
#>   ─  building ‘mrhorse_0.0.0.9000.tar.gz’
#>      
#> 
library(mrhorse)
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
- `return_fit`: return the fitted JAGS or Stan model object, default is
  FALSE
- `fixed_tau`: optional fixed value for tau, default is to estimate tau
- `omega`: correlation parameter for betaX and betaY, default is 0
  assuming independence (no sample overlap)

Output from the `mr_horse()` function include:

- `$MR_Estimate`: a data frame with the causal effect estimate (which is
  the posterior mean), standard deviation (i.e., the posterior standard
  deviation), upper and lower bounds of the 95% credible interval, and
  the R-hat value
- `$MR_Coda`: full MCMC samples for all parameters in `variable.names`,
  note that burnin iterations are discared
- `$Fit`: the fitted model object (if reuturn_fit was set to TRUE)

JAGS plotting tools can be used with the MCMC output, for example using
`traceplot()` and `densplot()`.

The function `mvmr_horse` implements the multivariable MVMR-Horse
method. The required inputs are as above, but where `betaX` is replaced
by the K columns `betaX1`, `betaX2`, …, and `betaXse` is replaced by the
K columns `betaX1se`, `betaX2se`, … .

## Examples

### Univariable MR

The built in data data_ex contains simulated genetic association
estimates between 100 genetic instruments and an exposure and outcome,
as well as their corresponding standard errors. This is taken from the
first replication of the first simulation study (that is, where 20% of
variants are pleiotropic, and pleiotropy is balanced). The dataset can
be analysed as follows.

``` r
set.seed(20230531)
MREx = mr_horse(data_ex, n.iter=1000, n.burnin=1000)
MREx$MR_Estimate
#>   Estimate    SD 2.5% quantile 97.5% quantile Rhat
#> 1    0.098 0.018         0.064          0.135    1
```

View diagnostic plots

``` r
coda::traceplot(MREx$MR_Coda[, "theta"])
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="50%" />

``` r
coda::densplot(MREx$MR_Coda[, "theta"])
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="50%" />

Note the R2jags and rstan package also have a traceplot and densplot
function

### Multivariable MR

The csv file dat_mv_ex.csv contains a dataframe containing simulated
genetic association estimates between 100 genetic instruments and two
exposures and an outcome, as well as their corresponding standard
errors. This is taken from the first replication of the multivariable
simulation study (that is, where 20% of variants are pleiotropic, and
pleiotropy is balanced). The dataset can be analysed as follows.

``` r
set.seed(100)
MVMREx = mvmr_horse(data_mv_ex, n.iter=1000, n.burnin=1000)
#> Fitting model with 2 exposures and 100 variants
MVMREx$MR_Estimate
#>   Parameter Estimate    SD 2.5% quantile 97.5% quantile  Rhat
#> 1  theta[1]    0.100 0.019         0.063          0.138 1.024
#> 2  theta[2]    0.106 0.018         0.070          0.140 1.011
```
