#' @title Diagnostic and Prognostic Meta-Analysis with Sensitivity Analysis
#'
#' @description
#' This package implements the meta-analysis of diagnostic test accuracy (DTA) studies and prognostic studies
#' by the summary receiver operating characteristic (SROC) curve and the area under the SROC curve (SAUC),
#' together with the likelihood-based sensitivity analysis for publication bias.
#' \itemize{
#'   \item \strong{Diagnostic studies}: the bivariate random-effects model of the logit-transformed
#'   sensitivity and specificity (Reitsma et al., 2005) \doi{10.1016/j.jclinepi.2005.02.022},
#'   and the sensitivity analysis for publication bias on the SROC curve and SAUC
#'   (Zhou et al., 2023) \doi{10.1002/sim.9643}.
#'   \item \strong{Prognostic studies}: the time-dependent SROC curve and SAUC, where the time-dependent
#'   sensitivity and specificity are reconstructed from the published survival probabilities
#'   (Hattori and Zhou, 2016) \doi{10.1002/sim.7029}
#'   and the sensitivity analysis for publication bias on the time-dependent SROC curve and SAUC
#'   (Zhou et al., 2025) \doi{10.1017/rsm.2025.14}.
#' }
#'
#' @section Fitting functions:
#' \tabular{lll}{
#'   \tab \strong{Without publication bias} \tab \strong{Sensitivity analysis for publication bias} \cr
#'   Diagnostic studies (SROC) \tab \code{\link{reitsma.sroc}} \tab \code{\link{srocsa.t}} \cr
#'   Prognostic studies (time-dependent SROC) \tab \code{\link{hattori.tsroc}} \tab \code{\link{tsrocsa.lrt}}
#' }
#' \code{\link{sauc}} calculates the SAUC and its confidence interval.
#'
#' @section Methods for the results:
#' \itemize{
#'   \item Objects of class \code{"sroc"} (from \code{\link{reitsma.sroc}} and \code{\link{hattori.tsroc}}):
#'   \code{\link{print.sroc}}, \code{\link{summary.sroc}} (printed by \code{\link{print.summary.sroc}}),
#'   and \code{\link{plot.sroc}}.
#'   \item Objects of class \code{"srocsa"} (from \code{\link{srocsa.t}} and \code{\link{tsrocsa.lrt}}):
#'   \code{\link{print.srocsa}}, \code{\link{summary.srocsa}}, and \code{\link{plot.srocsa}}.
#' }
#' For the time-dependent results, the plot methods draw one plot for each time point
#' (press <Return> for the next in an interactive session), or all the time points in one panel
#' (\code{overlay = TRUE}).
#'
#' @section Example data:
#' \itemize{
#'   \item \code{\link{IVD}}, \code{\link{CD64}}: meta-analyses of diagnostic studies.
#'   \item \code{\link{Ki67}}: meta-analysis of Ki-67 as a prognostic marker in early breast cancer.
#' }
#'
#' @references
#'
#' Hattori S, Zhou XH (2016).
#' Time-dependent summary receiver operating characteristics for meta-analysis of prognostic studies.
#' \emph{Statistics in Medicine}, 35(26), 4746--4763.
#' \doi{10.1002/sim.7029}
#'
#' Reitsma JB, Glas AS, Rutjes AW, Scholten RJ, Bossuyt PM, Zwinderman AH (2005).
#' Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews.
#' \emph{Journal of Clinical Epidemiology}, 58(10), 982--990.
#' \doi{10.1016/j.jclinepi.2005.02.022}
#'
#' Zhou Y, Huang A, Hattori S (2023).
#' A likelihood-based sensitivity analysis for publication bias on the summary receiver operating characteristic
#' in meta-analysis of diagnostic test accuracy.
#' \emph{Statistics in Medicine}, 42(6), 781--798.
#' \doi{10.1002/sim.9643}
#'
#' Zhou Y, Huang A, Hattori S (2025).
#' Sensitivity analysis for reporting bias on the time-dependent summary receiver operating characteristics curve
#' in meta-analysis of prognosis studies with time-to-event outcomes.
#' \emph{Research Synthesis Methods}, 16(3), 528--549.
#' \doi{10.1017/rsm.2025.14}
#'
#'
#' @docType package
#' @name DiagProgMeta-package
#' @importFrom stats integrate nlminb plogis pnorm qlogis qnorm uniroot
#' @importFrom stats complete.cases na.omit setNames
#' @importFrom utils head tail modifyList
#' @importFrom mixmeta mixmeta
#' @importFrom lifecycle deprecated
"_PACKAGE"
