#' @title Summary of 'sroc' objects
#'
#' @description Summarize the results of \code{\link{reitsma.sroc}} or \code{\link{hattori.tsroc}}
#' (for each time point): convergence of the optimizer, the log-likelihood,
#' the estimates of \eqn{(\mu_1, \mu_2, \tau_1, \tau_2, \rho)}{(mu1, mu2, tau1, tau2, rho)}
#' with standard errors (SE) and confidence intervals,
#' the summary sensitivity, specificity, and SAUC, and all the estimates (\code{par.all}).
#'
#' @param object
#' A fitted object returned by the function \code{\link{reitsma.sroc}} or \code{\link{hattori.tsroc}}.
#' @param ... Not used.
#'
#' @details
#' The SEs are the square roots of the diagonal elements of \code{var.ml}
#' (the inverse of the numerical Hessian matrix of the negative log-likelihood).
#' The confidence intervals are the Wald intervals of \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2},
#' of \eqn{\log\tau_1}{log(tau1)} and \eqn{\log\tau_2}{log(tau2)}, and of
#' \eqn{\mathrm{atanh}(\rho)}{atanh(rho)} (Fisher's z), back-transformed;
#' so they stay within the parameter space.
#' The confidence level is \code{ci.level} used in the fitting function (0.95 by default).
#'
#' @return An object of class \code{"summary.sroc"}: a list with
#' \describe{
#' \item{\code{call}}{The call of the fitting function (if recorded).}
#' \item{\code{time.dependent}}{\code{TRUE} for the results of \code{\link{hattori.tsroc}}.}
#' \item{\code{ci.level}}{The confidence level.}
#' \item{\code{fits}}{A list with one element for each fit (time point), each a list of
#' \code{ok}, \code{K} (number of studies), \code{conv}, \code{message}, \code{logLik},
#' \code{coef} (matrix of the estimates, SEs, and confidence limits),
#' \code{sauc}, \code{sens}, \code{spec} (estimate and confidence limits), and \code{par.all}.}
#' }
#' It is printed by \code{\link{print.summary.sroc}}.
#'
#' @examples
#' data(IVD)
#' fit <- reitsma.sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
#' summary(fit)
#'
#' ## the coefficient table as a matrix
#' s <- summary(fit)
#' s$fits[[1]]$coef
#'
#' ## for hattori.tsroc(), one element of s$fits for each time point
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
#')   
#' s2 <- summary(fit2)
#' s2$fits[["t = 36"]]$coef
#'
#' @seealso
#' \code{\link{print.sroc}},
#' \code{\link{print.summary.sroc}},
#' \code{\link{plot.sroc}}
#'
#' @method summary sroc
#' @export
summary.sroc <- function(object, ...) {

  fits  <- .sroc_fits_safe(object)
  lev   <- .sroc_ci_level(object)
  st    <- lapply(seq_along(fits), function(i)
    .sroc_fit_stats(fits[[i]], lev, n.study = .sroc_n_study(object, fits, i)))
  names(st) <- names(fits)

  structure(
    list(
      call           = attr(object, "call"),
      time.dependent = isTRUE(attr(fits, "td")),
      ci.level       = lev,
      fits           = st
    ),
    class = "summary.sroc"
  )
}
