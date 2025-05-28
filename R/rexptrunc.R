#### function 4: rexptrunc()

#' @title Truncated exponential distribution
#'
#' @description Random generation for the truncated exponential distribution
#'
#' @param n number of observations
#' @param range numeric vector of length two with bounds for the truncated distribution
#' @param rate rate for distribution
#'
#' @return Vector of length `n` of random samples from truncated exponential distribution
#'
#' @importFrom stats pexp runif qexp
#' @export

rexptrunc <- function(n, range, rate) {
  F.a <- pexp(min(range), rate = rate)
  F.b <- pexp(max(range), rate = rate)
  u <- runif(n, min = F.a, max = F.b)
  qexp(u, rate = rate)
}
