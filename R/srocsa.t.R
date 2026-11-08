#' @title Sensitivity analysis for publication bias on the SROC curve
#'
#' @description This function conducts the likelihood-based sensitivity analysis for publication bias
#' on the summary receiver operating characteristic (SROC) curve and the area under the SROC curve (SAUC)
#' in meta-analysis of diagnostic test accuracy (Zhou et al., 2023).
#' The bivariate random-effects model of the logit-transformed sensitivity and specificity is combined with
#' a selection function on the \eqn{t}-statistic of the linear combination of them
#' with the contrast vector \eqn{(c_1, c_2)}.
#' The model is fitted for each marginal selection probability \eqn{p} in \code{sa.p}.
#'
#' @param TP Specify the column of true positives (TP).
#' @param FN Specify the column of false negatives (FN).
#' @param TN Specify the column of true negatives (TN).
#' @param FP Specify the column of false positives (FP).
#'
#' @param data A data frame containing the columns of \code{TP}, \code{FN}, \code{TN}, and \code{FP}.
#' One row corresponds to one study. It is required.
#'
#' @param cc.value Value added to the cells for the ``continuity correction'', within (0, 1].
#' The default is 0.5.
#'
#' @param cc.type Type of the ``continuity correction'' (\code{c("single", "all")}).
#' The default is \code{"single"}.
#' \code{"single"}: \code{cc.value} is added to all the cells of the studies with at least one zero cell.
#' \code{"all"}: \code{cc.value} is added to all the cells of all the studies
#' if at least one study has a zero cell.
#'
#' @param sa.p A numeric vector of the marginal selection probabilities \eqn{p = P(\mathrm{select})}{p = P(select)},
#' within (0, 1]. The model is fitted for each value.
#' \eqn{p = 1} indicates no publication bias, and \eqn{\beta}{beta} and \eqn{\alpha}{alpha} are not estimated.
#' The default is 0.9.
#'
#' @param c11 Pre-specified value of \eqn{c_1^2}{c1^2} within [0, 1] (\eqn{c_2^2 = 1 - c_1^2}{c2^2 = 1 - c1^2}).
#' If \code{c11 = NULL} (default), \eqn{c_1^2}{c1^2} is estimated together with the other parameters.
#'
#' @param ci.level The confidence level of the confidence intervals, within (0, 1).
#' The default is 0.95.
#'
#' @param parallel.ncores The number of cores used for the parallel computation over \code{sa.p}
#' (PSOCK cluster of the package \pkg{parallel}); at most \code{length(sa.p)} cores are used.
#' The default is \code{1L} (no parallel computation).
#'
#' @param parset A list of the settings for the estimation, which updates the default settings.
#' Bad bounds may cause non-convergence.
#' \itemize{
#'     \item \code{mu.bound}: \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2} are searched within
#'     \code{[-mu.bound, mu.bound]}. The default is 5.
#'     \item \code{tau.bound}: \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2} are searched within
#'     \code{[eps, tau.bound]}. The default is 5.
#'     \item \code{beta.interval}: The lower and upper bounds of \eqn{\beta}{beta}. The default is \code{c(0, 5)}.
#'     \item \code{alpha.bound}: The root of \eqn{\alpha}{alpha} is searched starting from the interval
#'     \code{[-alpha.bound, alpha.bound]}, which is extended if necessary (see \code{\link[stats]{uniroot}}).
#'     The default is 5.
#'     \item \code{eps}: The lower bound of \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}. The default is \code{1e-5}.
#'   }
#' The element \code{inital.values} is not used in this function.
#'
#' @details
#' \describe{
#' \item{Bivariate model}{Let \eqn{y_{1i}}{y1i} and \eqn{y_{2i}}{y2i} be the logit-transformed sensitivity
#' and specificity of study \eqn{i}, and \eqn{v_{1i}}{v1i} and \eqn{v_{2i}}{v2i} be their variances.
#' \eqn{(y_{1i}, y_{2i})}{(y1i, y2i)} follows the bivariate normal distribution with mean \eqn{(\mu_1, \mu_2)}{(mu1, mu2)}
#' and variance-covariance matrix
#' \eqn{\mathrm{diag}(v_{1i}, v_{2i}) + \Omega}{diag(v1i, v2i) + Omega},
#' where \eqn{\Omega}{Omega} has the between-study standard deviations \eqn{\tau_1, \tau_2}{tau1, tau2}
#' and correlation \eqn{\rho}{rho}.}
#' \item{Selection function}{The probability of publishing study \eqn{i} is modeled as
#' \eqn{\Phi(\alpha + \beta t_i)}{Phi(alpha + beta * t_i)}, where
#' \eqn{t_i = (c_1 y_{1i} + c_2 y_{2i}) / \sqrt{c_1^2 v_{1i} + c_2^2 v_{2i}}}{t_i = (c1 * y1i + c2 * y2i) / sqrt(c1^2 * v1i + c2^2 * v2i)},
#' \eqn{c_1, c_2 \ge 0}{c1, c2 >= 0}, and \eqn{c_1^2 + c_2^2 = 1}{c1^2 + c2^2 = 1}.
#' For given \eqn{p} and the other parameters, \eqn{\alpha}{alpha} is determined such that
#' the mean of the inverse of the marginal selection probabilities of the observed studies equals \eqn{1/p}.}
#' \item{SAUC}{See \code{\link{sauc}}.}
#' \item{Initial values}{The initial values of \eqn{(\mu_1, \mu_2, \tau_1, \tau_2, \rho)}{(mu1, mu2, tau1, tau2, rho)}
#' are the estimates of the bivariate model without selection by \code{\link[mixmeta]{mixmeta}}.
#' The initial values of \eqn{\beta}{beta} and \eqn{c_1^2}{c1^2} are 0.1.}
#' \item{Continuity correction}{\url{https://en.wikipedia.org/wiki/Continuity_correction}}
#' }
#'
#' @return An object of class \code{"srocsa"}:
#' a list with one element for each value of \code{sa.p} and the element \code{call}.
#' Each element is the output of \code{\link[stats]{nlminb}}
#' (\code{par} is \code{(mu1, mu2, tau1, tau2, rho, beta)}, followed by \code{c11} when \code{c11 = NULL})
#' with the following additional elements, or \code{NULL} if the estimation failed.
#' \describe{
#' \item{\code{var.ml}}{The inverse of the numerical Hessian matrix of the negative log-likelihood.
#' When \eqn{p = 1}, only for \code{(mu1, mu2, tau1, tau2, rho)}.}
#' \item{\code{sauc.ci}}{SAUC and its confidence interval.}
#' \item{\code{beta}, \code{alpha}}{The estimates of \eqn{\beta}{beta} and \eqn{\alpha}{alpha}
#' (\code{NA} when \eqn{p = 1}).}
#' \item{\code{mu1.ci}, \code{mu2.ci}}{The estimates and confidence intervals of \eqn{\mu_1}{mu1} and \eqn{\mu_2}{mu2},
#' and of the sensitivity and specificity.}
#' \item{\code{par.all}}{All the estimates:
#' \code{mu1}, \code{mu2}, \code{tau1^2}, \code{tau2^2}, \code{tau12}, \code{c1^2}, \code{c2^2},
#' \code{beta}, \code{alpha}, \code{sauc}, \code{sens}, \code{spec}.}
#' }
#' The attribute \code{"data"} is a data frame of the study-level logit-transformed sensitivity (\code{y1})
#' and specificity (\code{y2}) and their variances (\code{v1}, \code{v2}),
#' used by \code{\link{plot.srocsa}} (\code{add.study = TRUE}).
#' The results can be shown by \code{\link{print.srocsa}}, \code{\link{summary.srocsa}}, and \code{\link{plot.srocsa}}.
#'
#' @references
#' Zhou Y, Huang A, Hattori S (2023).
#' A likelihood-based sensitivity analysis for publication bias on the summary receiver operating characteristic in meta-analysis of diagnostic test accuracy.
#' \emph{Statistics in Medicine}, 42(6), 781--798.
#' \doi{10.1002/sim.9643}
#'
#' @examples
#' ## load data
#' data(IVD)
#'
#' ## estimate c1^2
#' fit1 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = c(1, 0.8, 0.6))
#' fit1
#'
#' ## fix c1^2 = 0.5
#' fit2 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD,
#' sa.p = c(1, 0.8, 0.6), c11 = 0.5)
#' print(fit2)
#' summary(fit2)
#' plot(fit2)
#'
#' \donttest{
#' ## parallel computation
#' fit2p <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD,
#' sa.p = c(1, 0.8, 0.6), c11 = 0.5, parallel.ncores = 2L)
#' print(fit2p)
#' }
#'
#' ## continuity correction for all the studies
#' fit3 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = 0.8, cc.type = "all")
#' fit3
#'
#' @seealso
#' \code{\link{print.srocsa}},
#' \code{\link{summary.srocsa}},
#' \code{\link{plot.srocsa}},
#' \code{\link{reitsma.sroc}},
#' \code{\link{sauc}}
#'
#' @export
srocsa.t <- function(
  # data input
  TP, FN, TN, FP,
  data = NULL,
  cc.value = 0.5,
  cc.type = c("single", "all"),
  sa.p = 0.9,
  c11 = NULL,
  ci.level = 0.95,
  parallel.ncores = 1L,
  parset = list()
){

  ## record arguments ----
  this.call <- match.call()

  ## check data input ----
  have_data  <- !missing(data)
  if (!have_data) {
    stop("Data is missing or wrong.")
  }


  ## check parameter input ----
  cc.type <- match.arg(cc.type)

  if (cc.value <=0 | cc.value>1) {
    stop("The value of cc.value must be positive and within (0, 1]",  call. = FALSE)
  }

  if (!is.numeric(sa.p) || length(sa.p) == 0 || any(is.na(sa.p)) || any(sa.p <= 0 | sa.p > 1)) {
    stop("The values of sa.p (marginal selection probabilities) must be within (0, 1]", call. = FALSE)
  }

  if (!is.null(c11) && (c11 <0 | c11>1)) {
    stop("The value of c11 must be within [0, 1]",  call. = FALSE)
  }

  if (ci.level <=0 | ci.level>=1) {
    stop("The value of ci.level must be within (0, 1)",  call. = FALSE)
  }


  ## if data is used: require TP, FN, TN, FP ----
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


      TP = inputs$TP
      FN = inputs$FN
      FP = inputs$FP
      TN = inputs$TN


      ## update parameter setting
  parset <- modifyList(default.parset.sroc(), parset)

  data <- convert.dt.sroc(TP, FN, FP, TN, data=NULL, cc.value=cc.value, cc.type= cc.type)


 ## estimation ----

  n  <- nrow(data)
  y1 <- data$y1
  y2 <- data$y2
  v1 <- data$v1
  v2 <- data$v2


f <- function(pp){

  if(is.null(c11)) {

    sroc.pb.t.unfix(y1=y1, y2=y2, v1=v1, v2=v2, p=pp, ci.level = ci.level, parset = parset)

    } else {

      sroc.pb.t.fix(y1=y1, y2=y2, v1=v1, v2=v2, p=pp, c11 = c11, ci.level = ci.level, parset = parset)
    }

  }

  ## fit for each p (in parallel when parallel.ncores > 1)
  res <- .lapply_parallel(sa.p, f, ncores = parallel.ncores)

  ## LOGIT-DATA  -----
    # opt$l.data <- data
  res$call <- this.call

  class(res) <- "srocsa"

  ## information for print/summary/plot methods
  attr(res, "type")    <- "srocsa.t"
  attr(res, "sa.p")    <- sa.p
  attr(res, "n.study") <- n
  attr(res, "data")    <- data.frame(y1 = y1, y2 = y2, v1 = v1, v2 = v2)

  return(res)

}





