#' Print method for 'meta.pb' objects
#'
#' Internal. It prints \code{par.all} of a fitted model.
#' No function in this package currently returns an object of class \code{"meta.pb"}.
#'
#' @param x An object of class \code{"meta.pb"}.
#' @param digits The digits of the results.
#' @param ... Other arguments in function \code{\link[base]{print}}.
#'
#' @exportS3Method NULL
#' @keywords internal
#' @noRd
print.meta.pb <- function(x, digits = 3, ...){


  print(list(par.all = x$par.all), digits = digits, ...)

}
