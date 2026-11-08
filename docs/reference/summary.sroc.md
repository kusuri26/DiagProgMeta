# Summary of 'sroc' objects

Summarize the results of [`reitsma.sroc`](reference/reitsma.sroc.md) or
[`hattori.tsroc`](reference/hattori.tsroc.md) (for each time point):
convergence of the optimizer, the log-likelihood, the estimates of
\\(\mu_1, \mu_2, \tau_1, \tau_2, \rho)\\ with standard errors (SE) and
confidence intervals, the summary sensitivity, specificity, and SAUC,
and all the estimates (`par.all`).

## Usage

``` r
# S3 method for class 'sroc'
summary(object, ...)
```

## Arguments

- object:

  A fitted object returned by the function
  [`reitsma.sroc`](reference/reitsma.sroc.md) or
  [`hattori.tsroc`](reference/hattori.tsroc.md).

- ...:

  Not used.

## Value

An object of class `"summary.sroc"`: a list with

- `call`:

  The call of the fitting function (if recorded).

- `time.dependent`:

  `TRUE` for the results of
  [`hattori.tsroc`](reference/hattori.tsroc.md).

- `ci.level`:

  The confidence level.

- `fits`:

  A list with one element for each fit (time point), each a list of
  `ok`, `K` (number of studies), `conv`, `message`, `logLik`, `coef`
  (matrix of the estimates, SEs, and confidence limits), `sauc`, `sens`,
  `spec` (estimate and confidence limits), and `par.all`.

It is printed by
[`print.summary.sroc`](reference/print.summary.sroc.md).

## Details

The SEs are the square roots of the diagonal elements of `var.ml` (the
inverse of the numerical Hessian matrix of the negative log-likelihood).
The confidence intervals are the Wald intervals of \\\mu_1\\ and
\\\mu_2\\, of \\\log\tau_1\\ and \\\log\tau_2\\, and of
\\\mathrm{atanh}(\rho)\\ (Fisher's z), back-transformed; so they stay
within the parameter space. The confidence level is `ci.level` used in
the fitting function (0.95 by default).

## See also

[`print.sroc`](reference/print.sroc.md),
[`print.summary.sroc`](reference/print.summary.sroc.md),
[`plot.sroc`](reference/plot.sroc.md)

## Examples

``` r
data(IVD)
fit <- reitsma.sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
summary(fit)
#> SROC curve: bivariate random-effects model (Reitsma et al., 2005)
#> Publication bias is not considered
#> 
#> Call:
#> reitsma.sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
#> 
#> ============================================================
#> Results
#> ============================================================
#> Number of studies: 33
#> Convergence:       yes (relative convergence (4))
#> Log-likelihood:    -31.076
#> 
#> Parameters (tau1, tau2: between-study SDs; rho: correlation)
#> ------------------------------------------------------------
#>      Estimate    SE          95% CI
#> mu1     1.484 0.177  [1.136, 1.832]
#> mu2     1.824 0.154  [1.522, 2.125]
#> tau1    0.597 0.175  [0.336, 1.062]
#> tau2    0.826 0.117  [0.626, 1.090]
#> rho    -0.424 0.236 [-0.769, 0.111]
#> 
#> Summary measures
#> ------------------------------------------------------------
#>             Estimate         95% CI
#> Sensitivity    0.815 [0.757, 0.862]
#> Specificity    0.861 [0.821, 0.893]
#> SAUC           0.874 [0.806, 0.920]
#> 
#> All estimates (par.all)
#> ------------------------------------------------------------
#>    mu1    mu2 tau1^2 tau2^2  tau12   sauc   sens   spec 
#>  1.484  1.824  0.357  0.683 -0.209  0.874  0.815  0.861 
#> 

## the coefficient table as a matrix
s <- summary(fit)
s$fits[[1]]$coef
#>        Estimate        SE      lower     upper
#> mu1   1.4839917 0.1774776  1.1361421 1.8318413
#> mu2   1.8236493 0.1538733  1.5220632 2.1252355
#> tau1  0.5974413 0.1753902  0.3360539 1.0621394
#> tau2  0.8261895 0.1167857  0.6262653 1.0899360
#> rho  -0.4244263 0.2363535 -0.7690938 0.1114831

## for hattori.tsroc(), one element of s$fits for each time point
data(Ki67)
Ki67_HR <- Ki67$HR
Ki67_HR$lnHR     <- log(Ki67_HR$HR)
Ki67_HR$lnHR.var <- ((log(Ki67_HR$ci.up) - log(Ki67_HR$ci.low)) / (2 * qnorm(0.975)))^2
fit2 <- hattori.tsroc(
  time.points = c(24, 36, 60),
  study.os = study,
  time.os = t, n1.os = n1, n0.os = n0,
  s1.os = s1, s0.os = s0,
  data.os = Ki67$OS,
  study.lnhr = study, lnhr = lnHR, var.lnHR = lnHR.var,
  data.lnhr = Ki67_HR,
  time.med = mct_mo, n1.med = n1, n0.med = n0,
  s1.med = s1_mct, s0.med = s0_mct,
  data.mt = Ki67$MT
)   
s2 <- summary(fit2)
s2$fits[["t = 36"]]$coef
#>        Estimate         SE       lower       upper
#> mu1   0.6731155 0.18530620  0.30992207 1.036309012
#> mu2   0.2837232 0.11469250  0.05893001 0.508516337
#> tau1  0.6556385 0.20886375  0.35115678 1.224130751
#> tau2  0.5121485 0.09521945  0.35574533 0.737314255
#> rho  -0.8772330 0.16089307 -0.99156170 0.004680657
```
