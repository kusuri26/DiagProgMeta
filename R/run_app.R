#' @title Launch the DiagProgMeta Shiny app
#'
#' @description
#' Launches a Shiny app for the meta-analysis of diagnostic test accuracy
#' studies (SROC curve): fitting \code{\link{reitsma.sroc}} (without
#' publication bias) and \code{\link{srocsa.t}} (likelihood-based sensitivity
#' analysis for publication bias), on the built-in example data
#' (\code{\link{IVD}}, \code{\link{CD64}}) or on a user-uploaded CSV file of
#' study-level \code{TP}, \code{FN}, \code{FP}, \code{TN} counts.
#'
#' @param ... Further arguments passed to \code{\link[shiny]{runApp}}
#' (e.g. \code{launch.browser}, \code{port}).
#'
#' @return Called for its side effect (launches the app); does not return.
#'
#' @examples
#' if (interactive()) {
#'   run_app()
#' }
#'
#' @seealso \code{\link{reitsma.sroc}}, \code{\link{srocsa.t}}
#'
#' @export
run_app <- function(...) {

  pkgs <- c("shiny", "DT", "bslib")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "The following packages are required to run the app: ",
      paste(missing, collapse = ", "),
      ".\nInstall them with install.packages(c(", paste0('"', missing, '"', collapse = ", "), ")).",
      call. = FALSE
    )
  }

  app_dir <- system.file("shiny-app", package = "DiagProgMeta")
  if (app_dir == "") {
    stop("Could not find the shiny-app directory. Try re-installing DiagProgMeta.", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
