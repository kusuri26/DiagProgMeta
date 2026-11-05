#' @title Plot the SROC curves
#'
#' @description
#' Plot method for objects of class \code{"sroc"} returned by
#' \code{\link{reitsma.sroc}} (one SROC curve) and \code{\link{hattori.tsroc}}
#' (one time-dependent SROC curve for each time point).
#' For the output of \code{\link{hattori.tsroc}}, the curves of the time points are drawn
#' one per page by default; in an interactive session the user is prompted
#' ("Hit <Return> to see next plot") before each new page.
#' Alternatively, all the time points can be overlaid in one panel (\code{overlay = TRUE}).
#'
#' @param x An object of class \code{"sroc"} from \code{\link{reitsma.sroc}} or \code{\link{hattori.tsroc}}.
#' @param time Time points to plot (only for \code{\link{hattori.tsroc}} results).
#' Either numeric values of the time points (e.g., \code{c(24, 60)}) or integer indices
#' when \code{time.index = TRUE}. Default \code{NULL} plots all the time points.
#' @param time.index Logical; if \code{TRUE}, \code{time} is interpreted as the positions of the time points.
#' @param overlay Logical; if \code{TRUE}, all the selected time points are drawn in one panel.
#' @param ask Logical; whether to wait for <Return> before each new page.
#' Default \code{NULL} means \code{TRUE} when more than one page is drawn on an interactive device.
#' @param add.spoint Logical; whether to add the summary operating point.
#' @param spoint.region Logical; whether to add the confidence region of the summary point
#' (an ellipse on the logit scale, back-transformed).
#' @param plot.ci Logical; whether to add pointwise confidence bands of the SROC curve
#' (delta method using \code{var.ml}).
#' @param ci.level Confidence level of the bands and the region. Default 0.95.
#' @param add.study Logical; whether to add the observed studies
#' (\eqn{(1 - \mathrm{spec}_i, \mathrm{sens}_i)}{(1 - spec_i, sens_i)}, after continuity correction)
#' stored in the \code{data} element of the fit.
#' @param study.size Size of the study points: \code{"fixed"} (default) or \code{"precision"},
#' proportional to \eqn{1/\sqrt{v_1 + v_2}}{1/sqrt(v1 + v2)} relative to its mean.
#' @param study.pch,study.cex Symbol and (base) size of the study points.
#' @param study.col Color of the study points. Default \code{NULL} uses the curve color, semi-transparent.
#' @param add.sauc Logical; whether to show the SAUC (and its CI) in a legend.
#' @param legend.pos Position of the legend (see \code{\link[graphics]{legend}}).
#' @param sroc.cols Color(s) of the SROC curves. Recycled over time points when \code{overlay = TRUE}.
#' @param sroc.lty,sroc.lwd Line type and width of the SROC curves.
#' @param ci.lty,ci.lwd Line type and width of the confidence bands.
#' @param sp.pch,sp.cex Symbol and size of the summary points.
#' @param xlab,ylab Axis labels.
#' @param main Main title. Default \code{NULL} gives an automatic title (with the time point).
#' @param addon Logical; if \code{TRUE}, add to the current plot instead of starting a new one
#' (only when one panel is drawn).
#' @param n.grid Number of grid points to draw each curve.
#' @param ... Other graphical arguments passed to \code{\link[graphics]{plot.default}}.
#'
#' @return Invisibly, a list (one element per drawn fit) of data frames with the
#' curve coordinates (\code{fpr}, \code{tpr}, and \code{lb}, \code{ub} when \code{plot.ci = TRUE}).
#'
#' @examples
#' data(IVD)
#' fit1 <- reitsma.sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
#' plot(fit1, plot.ci = TRUE, spoint.region = TRUE)
#' plot(fit1, add.study = TRUE)                             # add the studies
#' plot(fit1, add.study = TRUE, study.size = "precision")   # point size by precision
#'
#' ## time-dependent SROC curves
#' data(Ki67)
#' Ki67_HR <- Ki67$HR
#' Ki67_HR$lnHR     <- log(Ki67_HR$HR)
#' Ki67_HR$lnHR.var <- ((log(Ki67_HR$ci.up) - log(Ki67_HR$ci.low)) / (2 * qnorm(0.975)))^2
#' fit2 <- hattori.tsroc(
#'   time.points = c(24, 36, 60),
#'   study.os = study,
#'   time.os = t, n1.os = n1, n0.os = n0,
#'   s1.os = s1, s0.os = s0,
#'   data.os = Ki67$OS,
#'   study.lnhr = study, lnhr = lnHR, var.lnHR = lnHR.var,
#'   data.lnhr = Ki67_HR,
#'   time.med = mct_mo, n1.med = n1, n0.med = n0,
#'   s1.med = s1_mct, s0.med = s0_mct,
#'   data.mt = Ki67$MT
#' )
#' plot(fit2, ask = FALSE)             # one page per time point, no prompt
#' plot(fit2, time = 36)               # a single time point
#' plot(fit2, overlay = TRUE)          # all time points in one panel
#' @importFrom grDevices dev.interactive devAskNewPage hcl.colors adjustcolor
#' @importFrom graphics plot.default lines points polygon legend abline
#' @importFrom stats plogis qlogis qnorm qchisq
#'
#' @export
plot.sroc <- function(
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
) {

  fits <- .sroc_fits(x)
  is.td <- attr(fits, "td")

  ## ---- select time points (hattori.tsroc) ----
  if (!is.null(time)) {
    if (!is.td) {
      warning("'time' is ignored for a single (non time-dependent) SROC fit.")
    } else if (time.index) {
      if (any(time < 1 | time > length(fits))) stop("'time' index out of range.")
      fits <- fits[time]
    } else {
      tt <- .sroc_time_values(names(fits))
      idx <- match(time, tt)
      if (anyNA(idx)) {
        stop("Time point(s) not found: ", paste(time[is.na(idx)], collapse = ", "),
             ". Available: ", paste(tt, collapse = ", "))
      }
      fits <- fits[idx]
    }
  }

  ## ---- drop failed fits ----
  ok <- vapply(fits, .sroc_fit_ok, logical(1))
  if (any(!ok)) {
    warning("Estimation failed for: ", paste(names(fits)[!ok], collapse = ", "),
            "; these are skipped.")
    fits <- fits[ok]
  }
  if (length(fits) == 0) stop("No successful fit to plot.")

  k <- length(fits)
  if (is.null(sroc.cols)) {
    sroc.cols <- if (overlay && k > 1) hcl.colors(k, palette = "Dark 3") else "black"
  }
  sroc.cols <- rep_len(sroc.cols, k)
  if (!is.null(study.col)) study.col <- rep_len(study.col, k)
  study.size <- match.arg(study.size)

  ## ---- <Return> prompt between pages ----
  n.pages <- if (overlay) 1L else k
  if (is.null(ask)) ask <- n.pages > 1 && dev.interactive()
  if (ask && !addon) {
    old.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(old.ask), add = TRUE)
  }
  if (addon && n.pages > 1) {
    warning("'addon = TRUE' with several pages: only the first page is added to the current plot.")
  }

  out <- vector("list", k)
  names(out) <- names(fits)

  if (overlay) {
    ## ---------- all curves in one panel ----------
    if (is.null(main)) main <- if (is.td) "Time-dependent SROC curves" else "SROC curve"
    if (!addon) .sroc_frame(xlab, ylab, main, ...)
    leg <- character(k)
    for (i in seq_len(k)) {
      out[[i]] <- .sroc_draw(fits[[i]], col = sroc.cols[i],
                             lty = sroc.lty, lwd = sroc.lwd,
                             plot.ci = plot.ci, ci.lty = ci.lty, ci.lwd = ci.lwd,
                             add.spoint = add.spoint, spoint.region = spoint.region,
                             sp.pch = sp.pch, sp.cex = sp.cex,
                             ci.level = ci.level, n.grid = n.grid,
                             add.study = add.study, study.size = study.size,
                             study.pch = study.pch, study.cex = study.cex,
                             study.col = if (is.null(study.col)) NULL else study.col[i])
      lab <- if (is.td) names(fits)[i] else "SROC"
      leg[i] <- if (add.sauc) paste0(lab, ": ", .sroc_sauc_label(fits[[i]])) else lab
    }
    if (k > 1 || add.sauc) {
      legend(legend.pos, legend = leg, col = sroc.cols, lty = sroc.lty, lwd = sroc.lwd,
             pch = if (add.spoint) sp.pch else NA, bty = "n", cex = 0.85)
    }

  } else {
    ## ---------- one panel per fit ----------
    for (i in seq_len(k)) {
      main.i <- if (!is.null(main)) {
        rep_len(main, k)[i]
      } else if (is.td) {
        paste0("Time-dependent SROC curve (", names(fits)[i], ")")
      } else {
        "SROC curve"
      }
      if (!(addon && i == 1)) .sroc_frame(xlab, ylab, main.i, ...)
      out[[i]] <- .sroc_draw(fits[[i]], col = sroc.cols[i],
                             lty = sroc.lty, lwd = sroc.lwd,
                             plot.ci = plot.ci, ci.lty = ci.lty, ci.lwd = ci.lwd,
                             add.spoint = add.spoint, spoint.region = spoint.region,
                             sp.pch = sp.pch, sp.cex = sp.cex,
                             ci.level = ci.level, n.grid = n.grid,
                             add.study = add.study, study.size = study.size,
                             study.pch = study.pch, study.cex = study.cex,
                             study.col = if (is.null(study.col)) NULL else study.col[i])
      if (add.sauc) {
        legend(legend.pos, legend = .sroc_sauc_label(fits[[i]]), bty = "n", cex = 0.9)
      }
    }
  }

  invisible(out)
}


## =====================================================================
## Internal helpers
## =====================================================================

## Turn the object into a named list of single fits; attr "td" = time-dependent?
.sroc_fits <- function(x) {
  if (inherits(x, "try-error")) stop("The estimation failed (object is a 'try-error').")
  x <- unclass(x)
  if (is.list(x) && is.numeric(x$par)) {
    ## reitsma.sroc(): a single nlminb output
    fits <- list(SROC = x)
    attr(fits, "td") <- FALSE
  } else if (is.list(x) && length(x) > 0) {
    ## hattori.tsroc(): a list of nlminb outputs, one per time point
    fits <- x
    if (is.null(names(fits))) names(fits) <- paste0("t", seq_along(fits))
    attr(fits, "td") <- TRUE
  } else {
    stop("Unrecognized 'sroc' object.")
  }
  fits
}

## ---------------------------------------------------------------------
## SROC parameters (mu1, mu2, tau1, tau2, rho) and their variance matrix
## ---------------------------------------------------------------------
## Point estimates: from `par.all` when it reports mu1, mu2, tau1^2, tau2^2, tau12
## (reitsma.sroc, hattori.tsroc, srocsa.t, tsrocsa.lrt), otherwise from `par`.
##
## Variance (for CI bands / regions): `var.ml` is on the scale of `par`, so the
## 5 elements of `par` for logit-sens / logit-spec are located (by name when
## available) and converted to (mu1, mu2, tau1, tau2, rho) by one of
##   "rho"    : (mu1, mu2, tau1, tau2, rho or rho1)        tau = SD, correlation
##   "cov.sd" : (mu1, mu2, tau1, tau2, tau12)             tau = SD, covariance
##   "cov.var": (mu1, mu2, tau1^2, tau2^2, tau12)         variances, covariance
## with the delta method (numerical Jacobian). The mapping is accepted only when
## it reproduces `par.all`; otherwise the variance is NULL (band skipped).

.sroc_est_from_parall <- function(pa) {
  if (is.null(pa) || is.null(names(pa))) return(NULL)
  need <- c("mu1", "mu2", "tau1^2", "tau2^2", "tau12")
  if (!all(need %in% names(pa))) return(NULL)
  v <- as.numeric(pa[need])
  if (any(!is.finite(v)) || v[3] <= 0 || v[4] <= 0) return(NULL)
  c(mu1 = v[1], mu2 = v[2], tau1 = sqrt(v[3]), tau2 = sqrt(v[4]),
    rho = v[5] / sqrt(v[3] * v[4]))
}

.sroc_maps <- list(
  rho     = function(t) c(t[1], t[2], t[3], t[4], t[5]),
  cov.sd  = function(t) c(t[1], t[2], t[3], t[4], t[5] / (t[3] * t[4])),
  cov.var = function(t) c(t[1], t[2], sqrt(t[3]), sqrt(t[4]), t[5] / sqrt(t[3] * t[4]))
)

## candidate (index in par, mapping) pairs, most specific first
.sroc_candidates <- function(par) {
  nm <- names(par)
  out <- list()
  if (!is.null(nm) && all(c("mu1", "mu2") %in% nm)) {
    i12 <- match(c("mu1", "mu2"), nm)
    r <- intersect(c("rho", "rho1"), nm)
    if (all(c("tau1", "tau2") %in% nm) && length(r)) {
      out[[length(out) + 1]] <- list(idx = c(i12, match(c("tau1", "tau2", r[1]), nm)), map = "rho")
    }
    if ("tau12" %in% nm) {
      if (all(c("tau1", "tau2") %in% nm)) {
        idx <- c(i12, match(c("tau1", "tau2", "tau12"), nm))
        out[[length(out) + 1]] <- list(idx = idx, map = "cov.sd")
        out[[length(out) + 1]] <- list(idx = idx, map = "cov.var")
      }
      if (all(c("tau1^2", "tau2^2") %in% nm)) {
        out[[length(out) + 1]] <- list(idx = c(i12, match(c("tau1^2", "tau2^2", "tau12"), nm)),
                                       map = "cov.var")
      }
    }
  }
  if (length(out) == 0) {
    ## unnamed: bivariate (mu1, mu2, tau1, tau2, rho, ...) or
    ## trivariate (mu1, mu2, mu3, tau1, tau2, tau3, x12, x23, x13, ...)
    idx <- if (length(par) >= 9) c(1L, 2L, 4L, 5L, 7L) else 1:5
    for (m in names(.sroc_maps)) out[[length(out) + 1]] <- list(idx = idx, map = m)
  }
  out
}

.sroc_jacobian <- function(f, t) {
  J <- matrix(0, 5, 5)
  for (j in 1:5) {
    h <- 1e-6 * max(1, abs(t[j]))
    tp <- t; tm <- t; tp[j] <- t[j] + h; tm[j] <- t[j] - h
    J[, j] <- (f(tp) - f(tm)) / (2 * h)
  }
  J
}

## returns list(est = c(mu1, mu2, tau1, tau2, rho), V = 5x5 or NULL, note = reason if V is NULL)
.sroc_extract <- function(fit) {
  par <- fit$par
  est.all <- .sroc_est_from_parall(fit$par.all)
  cands <- if (is.numeric(par) && length(par) >= 5) .sroc_candidates(par) else list()

  ## choose the mapping of `par`
  chosen <- NULL
  for (cd in cands) {
    if (anyNA(cd$idx) || max(cd$idx) > length(par)) next
    e <- suppressWarnings(.sroc_maps[[cd$map]](as.numeric(par[cd$idx])))
    if (any(!is.finite(e))) next
    if (is.null(est.all) || all(abs(e - est.all) <= 1e-4 * pmax(1, abs(est.all)))) {
      chosen <- c(cd, list(est = e)); break
    }
  }

  est <- if (!is.null(est.all)) est.all else if (!is.null(chosen)) chosen$est else rep(NA_real_, 5)
  est <- stats::setNames(as.numeric(est), c("mu1", "mu2", "tau1", "tau2", "rho"))

  ## variance on the (mu1, mu2, tau1, tau2, rho) scale
  V <- NULL; note <- NULL
  if (is.null(chosen)) {
    note <- "the parameters in 'par' could not be matched to 'par.all'"
  } else if (is.null(fit$var.ml)) {
    note <- "'var.ml' is not available"
  } else {
    Vf <- as.matrix(fit$var.ml)
    if (nrow(Vf) < max(chosen$idx) || ncol(Vf) < max(chosen$idx)) {
      note <- "'var.ml' does not cover the SROC parameters"
    } else {
      Vs <- Vf[chosen$idx, chosen$idx, drop = FALSE]
      if (any(!is.finite(Vs))) {
        note <- "'var.ml' is not finite"
      } else {
        J <- .sroc_jacobian(.sroc_maps[[chosen$map]], as.numeric(par[chosen$idx]))
        V <- J %*% Vs %*% t(J)
      }
    }
  }
  list(est = est, V = V, note = note)
}

.sroc_par <- function(fit) .sroc_extract(fit)$est

.sroc_fit_ok <- function(fit) {
  !is.null(fit) && !inherits(fit, "try-error") && is.list(fit) &&
    all(is.finite(.sroc_par(fit)))
}

## "t = 24" -> 24
.sroc_time_values <- function(nm) {
  suppressWarnings(as.numeric(sub("^\\s*t\\s*=\\s*", "", nm)))
}

.sroc_frame <- function(xlab, ylab, main, ...) {
  plot.default(NULL, xlim = c(0, 1), ylim = c(0, 1),
               xlab = xlab, ylab = ylab, main = main, ...)
  abline(0, 1, col = "grey80", lty = 3)
}

## SROC curve (Reitsma / Hattori parameterization) with delta-method band
## logit(TPR) = mu1 - (tau1 * rho / tau2) * (logit(FPR) + mu2)
.sroc_curve <- function(par, V = NULL, ci.level = 0.95, n.grid = 300) {
  mu1 <- par[1]; mu2 <- par[2]; t1 <- par[3]; t2 <- par[4]; r <- par[5]
  ## grid: uniform on the logit scale (to reach the corners) + uniform on [0, 1]
  ## (to keep the middle smooth); logit range +-30 is visually 0 and 1
  fpr <- sort(unique(c(
    plogis(seq(-30, 30, length.out = n.grid)),
    seq(1 / (n.grid + 1), n.grid / (n.grid + 1), length.out = n.grid)
  )))
  q   <- qlogis(fpr) + mu2
  g   <- mu1 - (t1 * r / t2) * q
  res <- data.frame(fpr = fpr, tpr = plogis(g))

  if (!is.null(V)) {
    ## gradient of g w.r.t. (mu1, mu2, tau1, tau2, rho)
    G <- cbind(1,
               -t1 * r / t2,
               -r / t2 * q,
               t1 * r / t2^2 * q,
               -t1 / t2 * q)
    se <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    z  <- qnorm(1 - (1 - ci.level) / 2)
    ## the band widens without limit in the tails (se grows with |logit(FPR)|);
    ## draw it only for |logit(FPR)| <= 5, i.e. FPR in about [0.007, 0.993]
    inside <- abs(qlogis(fpr)) <= 5
    res$lb <- ifelse(inside, plogis(g - z * se), NA_real_)
    res$ub <- ifelse(inside, plogis(g + z * se), NA_real_)
  }
  res
}

## Study points (1 - spec_i, sens_i) from a data frame with y1 = logit(sens), y2 = logit(spec)
.sroc_draw_studies <- function(d, col, study.size = "fixed", study.pch = 1, study.cex = 1) {
  if (is.null(d) || is.null(d$y1) || is.null(d$y2)) {
    warning("No study-level data in the fit; studies are not added. ",
            "Refit with the current version of the fitting function.", call. = FALSE)
    return(invisible(NULL))
  }
  ok  <- is.finite(d$y1) & is.finite(d$y2)
  cex <- rep(study.cex, sum(ok))
  if (study.size == "precision" && !is.null(d$v1) && !is.null(d$v2)) {
    w   <- 1 / sqrt(d$v1[ok] + d$v2[ok])
    w[!is.finite(w)] <- NA
    cex <- study.cex * sqrt(w / mean(w, na.rm = TRUE))
    cex[is.na(cex)] <- study.cex
  }
  points(plogis(-d$y2[ok]), plogis(d$y1[ok]), pch = study.pch, cex = cex, col = col)
  invisible(NULL)
}

.sroc_draw <- function(fit, col, lty, lwd, plot.ci, ci.lty, ci.lwd,
                       add.spoint, spoint.region, sp.pch, sp.cex,
                       ci.level, n.grid,
                       add.study = FALSE, study.size = "fixed",
                       study.pch = 1, study.cex = 1, study.col = NULL) {
  ex  <- .sroc_extract(fit)
  par <- ex$est

  ## observed studies (drawn first, so the curve is on top)
  if (add.study) {
    .sroc_draw_studies(fit$data,
                       col = if (is.null(study.col)) adjustcolor(col, alpha.f = 0.6) else study.col,
                       study.size = study.size, study.pch = study.pch, study.cex = study.cex)
  }

  V <- if (plot.ci || spoint.region) ex$V else NULL
  if ((plot.ci || spoint.region) && is.null(V)) {
    warning("Confidence band/region skipped: ", ex$note, ".", call. = FALSE)
  }

  cv <- .sroc_curve(par, V = if (plot.ci) V else NULL,
                    ci.level = ci.level, n.grid = n.grid)
  lines(cv$fpr, cv$tpr, col = col, lty = lty, lwd = lwd)
  if (plot.ci && !is.null(cv$lb)) {
    lines(cv$fpr, cv$lb, col = col, lty = ci.lty, lwd = ci.lwd)
    lines(cv$fpr, cv$ub, col = col, lty = ci.lty, lwd = ci.lwd)
  }

  ## confidence region of the summary point: ellipse for (mu1, mu2) on logit scale
  if (spoint.region && !is.null(V)) {
    S  <- V[1:2, 1:2]
    ev <- eigen(S, symmetric = TRUE)
    if (all(ev$values > 0)) {
      a  <- seq(0, 2 * pi, length.out = 200)
      rad <- sqrt(qchisq(ci.level, df = 2))
      E  <- rad * ev$vectors %*% diag(sqrt(ev$values)) %*% rbind(cos(a), sin(a))
      m1 <- par[1] + E[1, ]
      m2 <- par[2] + E[2, ]
      polygon(1 - plogis(m2), plogis(m1), border = col, lty = ci.lty, lwd = ci.lwd)
    }
  }

  if (add.spoint) {
    points(1 - plogis(par[2]), plogis(par[1]), pch = sp.pch, cex = sp.cex, col = col)
  }
  cv
}

## "SAUC = 0.853 (0.801, 0.894)" from fit$sauc.ci (or par.all as fallback)
.sroc_sauc_label <- function(fit, digits = 3) {
  s <- fit$sauc.ci
  v <- NULL
  if (!is.null(s)) v <- suppressWarnings(as.numeric(unlist(s)))
  v <- v[is.finite(v)]
  fmt <- function(z) formatC(z, format = "f", digits = digits)
  if (length(v) >= 3) {
    return(paste0("SAUC = ", fmt(v[1]), " (", fmt(v[2]), ", ", fmt(v[3]), ")"))
  }
  if (length(v) >= 1) return(paste0("SAUC = ", fmt(v[1])))
  pa <- fit$par.all
  if (!is.null(pa) && "sauc" %in% names(pa)) return(paste0("SAUC = ", fmt(pa[["sauc"]])))
  "SAUC = NA"
}
