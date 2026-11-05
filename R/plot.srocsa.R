#' @title Plot the sensitivity analysis of the SROC curves
#'
#' @description
#' Plot method for 'srocsa' objects.
#' This function plots the SROC curves and the summary points for each marginal selection probability \eqn{p},
#' from the sensitivity analysis by \code{\link{srocsa.t}} or \code{\link{tsrocsa.lrt}}.
#' For \code{\link{srocsa.t}}, one plot is created.
#' For \code{\link{tsrocsa.lrt}}, one plot is created at each time point, with the title \code{time = xx};
#' in an interactive session the user is prompted ("Hit <Return> to see next plot")
#' before each new time point.
#' The same \eqn{p} has the same color and line type at all the time points.
#' Alternatively, all the time points can be overlaid in one panel (\code{overlay = TRUE}):
#' the color indicates the time point and the line type indicates \eqn{p}.
#'
#' @param x
#' A fitted object returned by the function \code{\link{srocsa.t}} or \code{\link{tsrocsa.lrt}}.
#'
#' @param time
#' Time points to plot (only for the results from \code{\link{tsrocsa.lrt}}).
#' If \code{time = NULL} (default), all the time points are plotted.
#'
#' @param p
#' Values of the marginal selection probability \eqn{p} to plot.
#' If \code{p = NULL} (default), all the values are plotted.
#'
#' @param overlay Logical; if \code{TRUE}, all the selected time points of \code{\link{tsrocsa.lrt}}
#' are drawn in one panel (ignored for \code{\link{srocsa.t}}, which has one panel).
#'
#' @param ask Logical; whether to wait for <Return> before each new time point.
#' Default \code{NULL} means \code{TRUE} when more than one plot is drawn on an interactive device.
#'
#' @param sroc.cols Colors of the SROC curves: one for each \eqn{p},
#' or one for each time point when \code{overlay = TRUE}.
#' Default uses grey colors (black for the first \eqn{p}), and
#' \code{\link[grDevices]{hcl.colors}} (palette \code{"Dark 3"}) for the time points when \code{overlay = TRUE}.
#' @param sroc.lty Line types of the SROC curves, recycled over \eqn{p}.
#' When \code{overlay = TRUE} and it is not specified, \code{1, 2, 3, ...} are used for the values of \eqn{p}.
#' @param sroc.lwd Line width of the SROC curves.
#'
#' @param add.spoint Whether to add the summary points on the SROC curves.
#' @param sp.pch Type of the summary points.
#' @param sp.cex Size of the summary points.
#'
#' @param plot.ci Whether to add the pointwise confidence bands of the SROC curves
#' (delta method using \code{var.ml}).
#' @param spoint.region Whether to add the confidence regions of the summary points.
#' @param ci.level Confidence level of the bands and regions. Default 0.95.
#' @param ci.lty,ci.lwd Line type and width of the bands and regions.
#'
#' @param add.study Whether to add the observed (published) studies
#' (\eqn{(1 - \mathrm{spec}_i, \mathrm{sens}_i)}{(1 - spec_i, sens_i)}, after continuity correction),
#' stored in the attribute \code{"data"} of \code{x}.
#' @param study.size Size of the study points: \code{"fixed"} (default) or \code{"precision"},
#' proportional to \eqn{1/\sqrt{v_1 + v_2}}{1/sqrt(v1 + v2)} relative to its mean.
#' @param study.pch,study.cex Symbol and (base) size of the study points.
#' @param study.col Color of the study points. Default \code{NULL} uses grey,
#' or the (semi-transparent) color of each time point when \code{overlay = TRUE}.
#'
#' @param add.sauc Whether to show the SAUC (and its CI) for each \eqn{p} in the legend.
#' @param legend.pos Position of the legend. See \code{\link[graphics]{legend}}.
#' \code{NULL} suppresses the legend.
#' @param legend.cex Size of the legend text.
#'
#' @param main Title of the plot.
#' If \code{main = NULL} (default), \code{time = xx} is used for the time-dependent results
#' (\code{"Time-dependent SROC curves"} when \code{overlay = TRUE}),
#' and no title for the results from \code{srocsa.t}.
#' @param xlab Label of x-axis.
#' @param ylab Label of y-axis.
#' @param n.grid Number of grid points to draw each curve.
#'
#' @param ... Other arguments in function \code{\link[graphics]{plot.default}}.
#'
#' @details
#' For \eqn{p < 1} in \code{\link{tsrocsa.lrt}}, the SROC curve is calculated from
#' \eqn{(\mu_1, \mu_2, \tau_1, \tau_2, \rho_1)}{(mu1, mu2, tau1, tau2, rho1)} of the trivariate model,
#' i.e., the parameters of the logit-transformed sensitivity and specificity.
#'
#' @importFrom grDevices gray.colors dev.interactive devAskNewPage hcl.colors adjustcolor
#' @importFrom graphics legend
#'
#' @examples
#' ## load data
#' data(IVD)
#'
#' fit1 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = c(1, 0.8, 0.6))
#' plot(fit1)
#'
#' ## set the colors and add the studies
#' plot(fit1, sroc.cols = c("black", "blue", "red"), add.study = TRUE)
#'
#' ## only p = 1 and p = 0.6, with confidence bands
#' plot(fit1, p = c(1, 0.6), plot.ci = TRUE, sroc.cols = c("black", "red"))
#'
#' \donttest{
#' ## time-dependent results
#' data(Ki67)
#' Ki67_HR <- Ki67$HR
#' Ki67_HR$lnHR     <- log(Ki67_HR$HR)
#' Ki67_HR$lnHR.var <- ((log(Ki67_HR$ci.up) - log(Ki67_HR$ci.low)) / (2 * qnorm(0.975)))^2
#' fit2 <- tsrocsa.lrt(
#'   time.points = c(36, 60),
#'   study.os = study,
#'   time.os = t, n1.os = n1, n0.os = n0,
#'   s1.os = s1, s0.os = s0,
#'   data.os = Ki67$OS,
#'   study.lnhr = study, lnhr = lnHR, var.lnHR = lnHR.var,
#'   data.lnhr = Ki67_HR,
#'   time.med = mct_mo, n1.med = n1, n0.med = n0,
#'   s1.med = s1_mct, s0.med = s0_mct,
#'   data.mt = Ki67$MT,
#'   senp.p = c(1, 0.8)
#' )
#' plot(fit2, ask = FALSE)              # one plot per time point, no prompt
#' plot(fit2, time = 60)                # one time point
#' plot(fit2, overlay = TRUE)           # all time points in one panel (color: time, line type: p)
#' plot(fit2, overlay = TRUE, p = 0.8)  # all time points for one p
#' }
#'
#' @seealso
#' \code{\link{srocsa.t}},
#' \code{\link{tsrocsa.lrt}},
#' \code{\link{plot.sroc}}
#'
#' @return \code{x}, invisibly. Called for the side effect of plotting.
#'
#' @method plot srocsa
#' @export
plot.srocsa <- function(
  x,
  time = NULL,
  p = NULL,
  overlay = FALSE,
  ask = NULL,
  sroc.cols = NULL,
  sroc.lty = 1,
  sroc.lwd = 2,
  add.spoint = TRUE,
  sp.pch = 19,
  sp.cex = 1.2,
  plot.ci = FALSE,
  spoint.region = FALSE,
  ci.level = 0.95,
  ci.lty = 2,
  ci.lwd = 1,
  add.study = FALSE,
  study.size = c("fixed", "precision"),
  study.pch = 1,
  study.cex = 1,
  study.col = NULL,
  add.sauc = TRUE,
  legend.pos = "bottomright",
  legend.cex = 0.85,
  main = NULL,
  xlab = "FPR (1 - specificity)",
  ylab = "TPR (sensitivity)",
  n.grid = 300,
  ...
){

  study.size <- match.arg(study.size)
  panels <- .srocsa_panels(x)
  is.td  <- !is.null(panels[[1]]$time)

  ## select time points ----
  if (!is.null(time)) {
    if (!is.td) stop("`time` can be used only for the results from tsrocsa.lrt.", call. = FALSE)
    times <- vapply(panels, function(b) b$time, character(1))
    keep  <- times %in% as.character(time) |
      suppressWarnings(as.numeric(times)) %in% suppressWarnings(as.numeric(time))
    if (!any(keep)) {
      stop("No results at time = ", paste(time, collapse = ", "),
           ". Available: ", paste(times, collapse = ", "), call. = FALSE)
    }
    panels <- panels[keep]
  }

  ## all p values (to fix one color / line type for each p across panels) ----
  p.all <- sort(unique(unlist(lapply(panels, `[[`, "p"))), decreasing = TRUE)
  if (!is.null(p)) {
    p.sel <- p.all[vapply(p.all, function(v) any(abs(v - p) < 1e-8), logical(1))]
    if (length(p.sel) == 0) {
      stop("No results at p = ", paste(p, collapse = ", "),
           ". Available: ", paste(p.all, collapse = ", "), call. = FALSE)
    }
  } else {
    p.sel <- p.all
  }

  k    <- length(p.all)
  cols <- if (is.null(sroc.cols)) gray.colors(k, gamma = 1, start = 0, end = 0.6) else rep_len(sroc.cols, k)
  ltys <- rep_len(sroc.lty, k)
  ## when the user gives colors for the selected p only, match them to p.sel
  if (!is.null(sroc.cols) && !is.null(p) && length(sroc.cols) == length(p.sel)) {
    cols[match(p.sel, p.all)] <- sroc.cols
  }
  if (!is.null(p) && length(sroc.lty) == length(p.sel)) {
    ltys[match(p.sel, p.all)] <- sroc.lty
  }

  ## overlay only makes sense for several time points ----
  overlay <- isTRUE(overlay) && is.td

  ## <Return> prompt between time points ----
  n.pages <- if (overlay) 1L else length(panels)
  if (!is.null(main)) main <- rep_len(main, n.pages)
  if (is.null(ask)) ask <- n.pages > 1 && dev.interactive()
  if (ask) {
    old.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(old.ask), add = TRUE)
  }

  ## curves to draw in a panel: selected p and successful fits ----
  use.of <- function(b) {
    use <- which(vapply(b$p, function(v) any(abs(v - p.sel) < 1e-8), logical(1)))
    ok  <- vapply(b$fits[use], .sroc_fit_ok, logical(1))
    if (any(!ok)) {
      warning("Estimation failed for p = ", paste(b$p[use][!ok], collapse = ", "),
              if (!is.null(b$time)) paste0(" at time = ", b$time),
              "; these curves are skipped.", call. = FALSE)
    }
    use[ok]
  }

  draw.one <- function(fit, col, lty) {
    .sroc_draw(fit, col = col, lty = lty, lwd = sroc.lwd,
               plot.ci = plot.ci, ci.lty = ci.lty, ci.lwd = ci.lwd,
               add.spoint = add.spoint, spoint.region = spoint.region,
               sp.pch = sp.pch, sp.cex = sp.cex,
               ci.level = ci.level, n.grid = n.grid)
  }

  sauc.lab <- function(fits) vapply(fits, .sroc_sauc_label, character(1))

  if (overlay) {

    ## ---------- all time points in one panel: color = time, line type = p ----------
    kt    <- length(panels)
    tcols <- if (is.null(sroc.cols)) hcl.colors(kt, palette = "Dark 3") else rep_len(sroc.cols, kt)
    ## line type for each p: 1, 2, ... over the selected p (solid when only one p),
    ## or the user's sroc.lty (same rule as the panel mode)
    plty  <- ltys
    if (missing(sroc.lty)) plty[match(p.sel, p.all)] <- seq_along(p.sel)

    .sroc_frame(xlab, ylab, if (!is.null(main)) main[1] else "Time-dependent SROC curves", ...)

    leg.lab <- character(0); leg.col <- character(0); leg.lty <- numeric(0)

    for (i in seq_along(panels)) {
      b   <- panels[[i]]
      use <- use.of(b)
      if (length(use) == 0) {
        warning("No SROC curve to plot at time = ", b$time, call. = FALSE)
        next
      }
      if (add.study) {
        .sroc_draw_studies(b$data,
                           col = if (is.null(study.col)) adjustcolor(tcols[i], alpha.f = 0.5) else rep_len(study.col, kt)[i],
                           study.size = study.size, study.pch = study.pch, study.cex = study.cex)
      }
      j.p <- match(round(b$p[use], 10), round(p.all, 10))
      for (m in seq_along(use)) draw.one(b$fits[[use[m]]], col = tcols[i], lty = plty[j.p[m]])

      lab <- paste0("time = ", b$time, if (length(p.sel) > 1 || length(use) > 1) paste0(", p = ", b$p[use]))
      if (add.sauc) lab <- paste0(lab, ": ", sauc.lab(b$fits[use]))
      leg.lab <- c(leg.lab, lab)
      leg.col <- c(leg.col, rep(tcols[i], length(use)))
      leg.lty <- c(leg.lty, plty[j.p])
    }

    if (!is.null(legend.pos) && length(leg.lab)) {
      legend(legend.pos, legend = leg.lab, col = leg.col, lty = leg.lty, lwd = sroc.lwd,
             pch = if (add.spoint) sp.pch else NA, bty = "n", cex = legend.cex)
    }

  } else {

    ## ---------- one panel per time point: color / line type = p ----------
    for (i in seq_along(panels)) {

      b   <- panels[[i]]
      use <- use.of(b)

      main.i <- if (!is.null(main)) main[i] else if (!is.null(b$time)) paste0("time = ", b$time) else NULL
      .sroc_frame(xlab, ylab, main.i, ...)

      if (length(use) == 0) {
        warning("No SROC curve to plot", if (!is.null(b$time)) paste0(" at time = ", b$time),
                call. = FALSE)
        next
      }

      if (add.study) {
        .sroc_draw_studies(b$data, col = if (is.null(study.col)) "grey55" else study.col,
                           study.size = study.size, study.pch = study.pch, study.cex = study.cex)
      }

      j.col <- match(round(b$p[use], 10), round(p.all, 10))
      for (m in seq_along(use)) draw.one(b$fits[[use[m]]], col = cols[j.col[m]], lty = ltys[j.col[m]])

      if (!is.null(legend.pos)) {
        lab <- paste0("p = ", b$p[use])
        if (add.sauc) lab <- paste0(lab, ": ", sauc.lab(b$fits[use]))
        legend(legend.pos, legend = lab,
               col = cols[j.col], lty = ltys[j.col], lwd = sroc.lwd,
               pch = if (add.spoint) sp.pch else NA,
               bty = "n", cex = legend.cex)
      }
    }
  }

  invisible(x)
}


## =====================================================================
## Internal: split an 'srocsa' object into panels
## each panel: list(time = NULL or "36", p = numeric, fits = list, data = data.frame or NULL)
## =====================================================================
.srocsa_panels <- function(x) {

  type <- attr(x, "type")
  xx   <- unclass(x)

  ## structure check when the attribute "type" is missing
  if (is.null(type)) {
    elems <- if (is.null(names(xx))) xx else xx[names(xx) != "call"]
    nested <- length(elems) > 0 && all(vapply(elems, function(e)
      is.list(e) && !inherits(e, "try-error") && is.null(e$par), logical(1)))
    type <- if (nested) "tsrocsa.lrt" else "srocsa.t"
  }

  if (identical(type, "tsrocsa.lrt")) {

    times <- names(xx)
    if (is.null(times)) times <- as.character(seq_along(xx))

    lapply(seq_along(xx), function(i) {
      e  <- xx[[i]]
      pv <- suppressWarnings(as.numeric(names(e)))
      if (is.null(names(e)) || anyNA(pv)) pv <- rep(NA_real_, length(e))
      list(time = times[i], p = pv, fits = unname(as.list(e)), data = attr(e, "data"))
    })

  } else {

    nm   <- names(xx)
    fits <- if (is.null(nm)) xx else xx[nm != "call"]
    pv   <- attr(x, "sa.p")
    if (is.null(pv) || length(pv) != length(fits)) {
      warning("The values of p are not found; the curves are numbered instead.", call. = FALSE)
      pv <- seq_along(fits)
    }
    list(list(time = NULL, p = as.numeric(pv), fits = unname(fits), data = attr(x, "data")))
  }
}
