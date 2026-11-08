# Plot the SROC curves

Plot method for objects of class `"sroc"` returned by
[`reitsma.sroc`](reference/reitsma.sroc.md) (one SROC curve) and
[`hattori.tsroc`](reference/hattori.tsroc.md) (one time-dependent SROC
curve for each time point). For the output of
[`hattori.tsroc`](reference/hattori.tsroc.md), the curves of the time
points are drawn one per page by default; in an interactive session the
user is prompted ("Hit \<Return\> to see next plot") before each new
page. Alternatively, all the time points can be overlaid in one panel
(`overlay = TRUE`).

## Usage

``` r
# S3 method for class 'sroc'
plot(
  x,
  time = NULL,
  time.index = FALSE,
  overlay = FALSE,
  ask = NULL,
  add.spoint = TRUE,
  spoint.region = FALSE,
  plot.ci = FALSE,
  ci.level = 0.95,
  add.study = FALSE,
  study.size = c("fixed", "precision"),
  study.pch = 1,
  study.cex = 1,
  study.col = NULL,
  add.sauc = TRUE,
  legend.pos = "bottomright",
  sroc.cols = NULL,
  sroc.lty = 1,
  sroc.lwd = 2,
  ci.lty = 2,
  ci.lwd = 1,
  sp.pch = 19,
  sp.cex = 1.2,
  xlab = "FPR (1 - specificity)",
  ylab = "TPR (sensitivity)",
  main = NULL,
  addon = FALSE,
  n.grid = 300,
  ...
)
```

## Arguments

- x:

  An object of class `"sroc"` from
  [`reitsma.sroc`](reference/reitsma.sroc.md) or
  [`hattori.tsroc`](reference/hattori.tsroc.md).

- time:

  Time points to plot (only for
  [`hattori.tsroc`](reference/hattori.tsroc.md) results). Either numeric
  values of the time points (e.g., `c(24, 60)`) or integer indices when
  `time.index = TRUE`. Default `NULL` plots all the time points.

- time.index:

  Logical; if `TRUE`, `time` is interpreted as the positions of the time
  points.

- overlay:

  Logical; if `TRUE`, all the selected time points are drawn in one
  panel.

- ask:

  Logical; whether to wait for \<Return\> before each new page. Default
  `NULL` means `TRUE` when more than one page is drawn on an interactive
  device.

- add.spoint:

  Logical; whether to add the summary operating point.

- spoint.region:

  Logical; whether to add the confidence region of the summary point (an
  ellipse on the logit scale, back-transformed).

- plot.ci:

  Logical; whether to add pointwise confidence bands of the SROC curve
  (delta method using `var.ml`).

- ci.level:

  Confidence level of the bands and the region. Default 0.95.

- add.study:

  Logical; whether to add the observed studies (\\(1 - \mathrm{spec}\_i,
  \mathrm{sens}\_i)\\, after continuity correction) stored in the `data`
  element of the fit.

- study.size:

  Size of the study points: `"fixed"` (default) or `"precision"`,
  proportional to \\1/\sqrt{v_1 + v_2}\\ relative to its mean.

- study.pch, study.cex:

  Symbol and (base) size of the study points.

- study.col:

  Color of the study points. Default `NULL` uses the curve color,
  semi-transparent.

- add.sauc:

  Logical; whether to show the SAUC (and its CI) in a legend.

- legend.pos:

  Position of the legend (see
  [`legend`](https://rdrr.io/r/graphics/legend.html)).

- sroc.cols:

  Color(s) of the SROC curves. Recycled over time points when
  `overlay = TRUE`.

- sroc.lty, sroc.lwd:

  Line type and width of the SROC curves.

- ci.lty, ci.lwd:

  Line type and width of the confidence bands.

- sp.pch, sp.cex:

  Symbol and size of the summary points.

- xlab, ylab:

  Axis labels.

- main:

  Main title. Default `NULL` gives an automatic title (with the time
  point).

- addon:

  Logical; if `TRUE`, add to the current plot instead of starting a new
  one (only when one panel is drawn).

- n.grid:

  Number of grid points to draw each curve.

- ...:

  Other graphical arguments passed to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, a list (one element per drawn fit) of data frames with the
curve coordinates (`fpr`, `tpr`, and `lb`, `ub` when `plot.ci = TRUE`).

## Examples

``` r
data(IVD)
fit1 <- reitsma.sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
plot(fit1, plot.ci = TRUE, spoint.region = TRUE)

plot(fit1, add.study = TRUE)                             # add the studies

plot(fit1, add.study = TRUE, study.size = "precision")   # point size by precision


## time-dependent SROC curves
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
plot(fit2, ask = FALSE)             # one page per time point, no prompt



plot(fit2, time = 36)               # a single time point

plot(fit2, overlay = TRUE)          # all time points in one panel
```
