#' Compute NCV
#'
#' @param object An object.
#' @param ... Passed to methods.
#' @export
NCV <- function(object, ...) {
  UseMethod("NCV")
}