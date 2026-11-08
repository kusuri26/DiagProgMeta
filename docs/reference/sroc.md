# Bivariate random-effects model for the SROC curve without publication bias

This function fits the bivariate random-effects model of the
logit-transformed sensitivity and specificity (Reitsma et al., 2005) by
the maximum likelihood method, and calculates the SAUC and its
confidence interval. Publication bias is not considered.

## Usage

``` r
sroc(
  TP,
  FN,
  FP,
  TN,
  data = NULL,
  cc.value = 0.5,
  cc.type = c("single", "all"),
  ci.level = 0.95,
  parset = list()
)
```

## Arguments

- TP:

  Specify the column of true positives (TP).

- FN:

  Specify the column of false negatives (FN).

- FP:

  Specify the column of false positives (FP).

- TN:

  Specify the column of true negatives (TN).

- data:

  A data frame containing the columns of `TP`, `FN`, `FP`, and `TN`. One
  row corresponds to one study. If `data = NULL` (default), `TP`, `FN`,
  `FP`, and `TN` are evaluated in the calling environment (e.g., numeric
  vectors).

- cc.value:

  Value added to the cells for the “continuity correction”, within (0,
  1\]. The default is 0.5.

- cc.type:

  Type of the “continuity correction” (`c("single", "all")`). The
  default is `"single"`. `"single"`: `cc.value` is added to all the
  cells of the studies with at least one zero cell. `"all"`: `cc.value`
  is added to all the cells of all the studies if at least one study has
  a zero cell.

- ci.level:

  The confidence level of the confidence intervals, within (0, 1). The
  default is 0.95.

- parset:

  A list of the settings for the estimation, which updates the default
  settings. Bad initial values and bounds may cause non-convergence.

  - `inital.values`: The initial values of
    `c(mu1, mu2, tau1, tau2, rho)`. The default is
    `c(0.1, 0.1, 0.1, 0.1, -0.1)`. If it is `NULL` or not a vector of
    length 5, the estimates by
    [`mixmeta`](https://rdrr.io/pkg/mixmeta/man/mixmeta.html) are used.

  - `mu.bound`: \\\mu_1\\ and \\\mu_2\\ are searched within
    `[-mu.bound, mu.bound]`. The default is 5.

  - `tau.bound`: \\\tau_1\\ and \\\tau_2\\ are searched within
    `[eps, tau.bound]`. The default is 5.

  - `eps`: The lower bound of \\\tau_1\\ and \\\tau_2\\. The default is
    `1e-5`.

## Value

The output of [`nlminb`](https://rdrr.io/r/stats/nlminb.html) (`par` is
`(mu1, mu2, tau1, tau2, rho)`; `tau1` and `tau2` are standard
deviations) with the following additional elements. If
[`nlminb`](https://rdrr.io/r/stats/nlminb.html) failed, an object of
class `"try-error"`.

- `var.ml`:

  The inverse of the numerical Hessian matrix of the negative
  log-likelihood.

- `sauc.ci`:

  SAUC and its confidence interval.

- `mu1.ci`, `mu2.ci`:

  The estimates and confidence intervals of \\\mu_1\\ and \\\mu_2\\, and
  of the sensitivity and specificity.

- `par.all`:

  All the estimates: `mu1`, `mu2`, `tau1^2`, `tau2^2`, `tau12`, `sauc`,
  `sens`, `spec`.

## Details

- Bivariate model:

  See details in [`srocsa.t`](reference/srocsa.t.md).

- SAUC:

  See [`sauc`](reference/sauc.md).

- Continuity correction:

  <https://en.wikipedia.org/wiki/Continuity_correction>

## References

Reitsma JB, Glas AS, Rutjes AW, Scholten RJ, Bossuyt PM, Zwinderman AH
(2005). Bivariate analysis of sensitivity and specificity produces
informative summary measures in diagnostic reviews. *Journal of Clinical
Epidemiology*, 58(10), 982–990.
[doi:10.1016/j.jclinepi.2005.02.022](https://doi.org/10.1016/j.jclinepi.2005.02.022)

## See also

[`srocsa.t`](reference/srocsa.t.md), [`sauc`](reference/sauc.md)

## Examples

``` r
## load data
data(IVD)

fit <- sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
fit$par.all
#>        mu1        mu2     tau1^2     tau2^2      tau12       sauc       sens 
#>  1.4839917  1.8236493  0.3569361  0.6825891 -0.2094967  0.8735863  0.8151747 
#>       spec 
#>  0.8610034 
fit$sauc.ci
#>      sauc   sauc.lb   sauc.ub 
#> 0.8735863 0.8057723 0.9200719 
```
