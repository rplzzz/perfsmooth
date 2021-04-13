#' Evaluate polynomials at specified x values using synthetic substitution
#'
#' Polynomials are specified by a vector of coefficients \eqn{a_n}, with \eqn{a_0} in
#' the first vector entry, and counting upward from there.
#'
#' @param x Vector of values to evaluate
#' @param a Vector of polynomial coefficients
#' @return Vector of P(x) values, for all x
#' @export
polyeval <- function(x, a)
{
  px <- rep(0, length(x))
  d <- length(a)              # actually 1+degree
  for (i in seq(d,2)) {
    px <- (px + a[i])*x
  }
  px + a[1]
}

#' @describeIn polyeval Evaluate a sum of polynomials at specified values
#'
#' @param am Matrix of polynomial coefficients; each column is one polynomial
#' @export
polysumeval <- function(x, am)
{
  a <- apply(am, 1, sum)
  polyeval(x, a)
}
