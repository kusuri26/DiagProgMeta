#' Sensitivity Analysis for Publication Bias on the SROC Curve
#'
#' This package implements the likelihood-based sensitivity analysis for publication bias
#' on the summary receiver operating characteristic (SROC) curve and the area under the SROC curve (SAUC)
#' in meta-analysis of diagnostic test accuracy, see Zhou et al. (2023) \doi{10.1002/sim.9643}.
#' It also provides the sensitivity analysis on the time-dependent SROC curve in meta-analysis of prognostic studies,
#' where the time-dependent sensitivity and specificity are reconstructed from the published survival probabilities,
#' see Combescure et al. (2016) \doi{10.1177/0962280212464542}.
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{srocsa.t}} Sensitivity analysis for publication bias on the SROC curve and SAUC.
#'   \item \code{\link{tsrocsa.lrt}} Sensitivity analysis for publication bias on the time-dependent SROC curve and SAUC.
#'   \item \code{\link{sroc}} Bivariate random-effects model for the SROC curve without publication bias.
#'   \item \code{\link{tsroc}} Time-dependent SROC curve without publication bias.
#'   \item \code{\link{sauc}} SAUC and its confidence interval.
#' }
#'
#' Results output functions:
#' \itemize{
#'   \item \code{\link{print.srocsa}} Print the main results from 'srocsa' objects.
#'   \item \code{\link{summary.srocsa}} Summarize the results from 'srocsa' objects.
#'   \item \code{\link{plot.srocsa}} Plot the SROC curves from 'srocsa' objects.
#' }
#'
#' Example data:
#' \itemize{
#'   \item \code{\link{IVD}}, \code{\link{CD64}} Meta-analyses of diagnostic studies.
#' }
#'
#' @importFrom stats integrate nlminb plogis pnorm qlogis qnorm uniroot
#' @importFrom stats complete.cases na.omit setNames
#' @importFrom utils head tail modifyList
#' @importFrom mixmeta mixmeta
#' @importFrom lifecycle deprecated
"_PACKAGE"
