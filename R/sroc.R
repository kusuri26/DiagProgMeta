#' @title Bivariate random-effects model for the SROC curve without publication bias
#'
#' @description This function fits the bivariate random-effects model of the logit-transformed
#' sensitivity and specificity (Reitsma et al., 2005) by the maximum likelihood method,
#' and calculates the SAUC and its confidence interval.
#' Publication bias is not considered.
#'
#' @inheritParams srocsa.t
#'
#' @param data A data frame containing the columns of \code{TP}, \code{FN}, \code{FP}, and \code{TN}.
#' One row corresponds to one study.
#' If \code{data = NULL} (default), \code{TP}, \code{FN}, \code{FP}, and \code{TN} are evaluated
#' in the calling environment (e.g., numeric vectors).
#'
#' @param parset A list of the settings for the estimation, which updates the default settings.
#' Bad initial values and bounds may cause non-convergence.
#' \itemize{
#'     \item \code{inital.values}: The initial values of \code{c(mu1, mu2, tau1, tau2, rho)}.
#'     The default is \code{c(0.1, 0.1, 0.1, 0.1, -0.1)}.
#'     If it is \code{NULL} or not a vector of length 5, the estimates by \code{\link[mixmeta]{mixmeta}} are used.
#'     \item \code{mu.bound}: \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2} are searched within
#'     \code{[-mu.bound, mu.bound]}. The default is 5.
#'     \item \code{tau.bound}: \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2} are searched within
#'     \code{[eps, tau.bound]}. The default is 5.
#'     \item \code{eps}: The lower bound of \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}. The default is \code{1e-5}.
#'   }
#'
#' @details
#' \describe{
#' \item{Bivariate model}{See details in \code{\link{srocsa.t}}.}
#' \item{SAUC}{See \code{\link{sauc}}.}
#' \item{Continuity correction}{\url{https://en.wikipedia.org/wiki/Continuity_correction}}
#' }
#'
#' @return The output of \code{\link[stats]{nlminb}}
#' (\code{par} is \code{(mu1, mu2, tau1, tau2, rho)}; \code{tau1} and \code{tau2} are standard deviations)
#' with the following additional elements.
#' If \code{\link[stats]{nlminb}} failed, an object of class \code{"try-error"}.
#' \describe{
#' \item{\code{var.ml}}{The inverse of the numerical Hessian matrix of the negative log-likelihood.}
#' \item{\code{sauc.ci}}{SAUC and its confidence interval.}
#' \item{\code{mu1.ci}, \code{mu2.ci}}{The estimates and confidence intervals of \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2},
#' and of the sensitivity and specificity.}
#' \item{\code{par.all}}{All the estimates:
#' \code{mu1}, \code{mu2}, \code{tau1^2}, \code{tau2^2}, \code{tau12}, \code{sauc}, \code{sens}, \code{spec}.}
#' }
#'
#' @references
#' Reitsma JB, Glas AS, Rutjes AW, Scholten RJ, Bossuyt PM, Zwinderman AH (2005).
#' Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews.
#' \emph{Journal of Clinical Epidemiology}, 58(10), 982--990.
#' \doi{10.1016/j.jclinepi.2005.02.022}
#'
#' @examples
#' ## load data
#' data(IVD)
#'
#' fit <- sroc(TP = TP, FN = FN, FP = FP, TN = TN, data = IVD)
#' fit$par.all
#' fit$sauc.ci
#'
#' @seealso
#' \code{\link{srocsa.t}},
#' \code{\link{sauc}}
#'
#' @export
sroc <- function(
  TP, FN, FP, TN,
  data = NULL,
  cc.value = 0.5,
  cc.type = c("single", "all"),
  ci.level = 0.95,
  parset = list()
) {
  cc.type <- match.arg(cc.type)
  inputs <- resolve_inputs(
    exprs = list(
      TP  = substitute(TP),
      FN  = substitute(FN),
      FP  = substitute(FP),
      TN  = substitute(TN)
    ),
    data = data,
    env = parent.frame(),
    numeric = FALSE,
    same.length = TRUE
  )


  TP <- inputs$TP
  FN <- inputs$FN
  FP <- inputs$FP
  TN <- inputs$TN


  ## update parameter setting
  parset <- modifyList(default.parset.sroc(), parset)

  data.conv <- convert.dt.sroc(TP, FN, FP, TN, data = NULL, cc.value = cc.value, cc.type = cc.type)

  meta.sroc.mle(data.conv$y1, data.conv$y2, data.conv$v1, data.conv$v2, ci.level = ci.level, parset = parset)
}
