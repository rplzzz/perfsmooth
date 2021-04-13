#' Evaluate a Legendre series at a specified set of points
#'
#' Given the series coefficients, evaluate the Legendre series at the specified
#' points.
#'
#' The Legendre polynomials are defined on the interval [-1,1].  A function defined
#' on a different (finite) interval can be rescaled to work with the Legendre polynomials
#' using the \code{interval} parameter.  If the interval is omitted, it is assumed to
#' be [-1,1]
#'
#' @param x The values at which to evaluate the series
#' @param lc The coefficients of the Legendre polynomials in the series, starting
#' with P_0 and counting upward.  Any omitted coefficients are filled in with zeros.
#' @param interval Boundaries of the x interval (see details)
#' @export
lserieseval <- function(x, lc, interval)
{
  ## Rescale, if necessary
  if(!missing(interval)) {
    stopifnot(length(interval) == 2)
    L <- interval[2] - interval[1]
    S <- 2*L
    x <- (x-interval[1])*S - 1
  }

  ## Pad coefficients, if necessary
  if (length(lc) < ncol(legendrecoef)) {
    lc <- c(lc, rep(0, ncol(legendrecoef)-length(lc)))
  }
  else if(length(lc) > ncol(legendrecoef)) {
    warning('More than ', ncol(legendrecoef),
            ' series coefficients supplied. Excess will be dropped.')
    lc <- lc[seq(1,ncol(legendrecoef))]
  }

  am <- legendrecoef * rep(lc, rep(nrow(legendrecoef), length(lc)))

  polysumeval(x, am)
}
