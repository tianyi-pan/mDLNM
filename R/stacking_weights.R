#' Compute stacking weights
#' @description This function is from loo package: https://github.com/stan-dev/loo/blob/master/R/loo_model_weights.R
#'
#' @export
stacking_weights <- function(lpd_point,
                             optim_method = "BFGS",
                             optim_control = list()) {
  ##
  stopifnot(is.matrix(lpd_point))
  N <- nrow(lpd_point)
  K <- ncol(lpd_point)
  if (K < 2) {
    stop("At least two models are required for stacking weights.")
  }

  negative_log_score_loo <- function(w) {
    # objective function: log score
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    # avoid over- and underflows using log weights and rowLogSumExps
    sum <- sum(matrixStats::rowLogSumExps(sweep(lpd_point[1:N,], 2, log(w_full), '+')))
    return(-as.numeric(sum))
  }

  gradient <- function(w) {
    # gradient of the objective function
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    grad <- rep(0, K - 1)
    # avoid over- and underflows using log weights, rowLogSumExps,
    # and by subtracting the row maximum of lpd_point
    mlpd <- matrixStats::rowMaxs(lpd_point)
    for (k in 1:(K - 1)) {
      grad[k] <- sum((exp(lpd_point[, k] - mlpd) - exp(lpd_point[, K] - mlpd)) / exp(matrixStats::rowLogSumExps(sweep(lpd_point, 2, log(w_full), '+')) - mlpd))
    }
    return(-grad)
  }

  ui <- rbind(rep(-1, K - 1), diag(K - 1))  # K-1 simplex constraint matrix
  ci <- c(-1, rep(0, K - 1))
  w <- constrOptim(
    theta = rep(1 / K, K - 1),
    f = negative_log_score_loo,
    grad = gradient,
    ui = ui,
    ci = ci,
    method = optim_method,
    control = optim_control
  )$par

  wts <- structure(
    c(w, 1 - sum(w)),
    names = paste0("model", 1:K),
    class = c("stacking_weights")
  )

  return(wts)
}
