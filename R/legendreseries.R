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

  am <- lpolycoef(lc)

  polysumeval(x, am)
}

#' Scale Legendre polynomials by the coefficients of a Legendre series
#'
#' @param lc Coefficients of the Legendre series
#' @keywords internal
lpolycoef <- function(lc)
{
  ## Pad coefficients, if necessary
  if (length(lc) < ncol(legendrecoef)) {
    lc <- c(lc, rep(0, ncol(legendrecoef)-length(lc)))
  }
  else if(length(lc) > ncol(legendrecoef)) {
    warning('More than ', ncol(legendrecoef),
            ' series coefficients supplied. Excess will be dropped.')
    lc <- lc[seq(1,ncol(legendrecoef))]
  }

  ## Scale the coefficients of each polynomial by the coefficient corresponding
  ## to that polynomial.
  legendrecoef * rep(lc, rep(nrow(legendrecoef), length(lc)))
}


#' Find the extrema of the derivative of a Legendre series
#'
#' Find the min and max values, over the interval [-1,1] of a Legendre series,
#' given the coefficients of the series.
#'
#' We find these values by finding the real zeros of the second derivative of
#' the series and evaluating the first derivative at those points, plus the
#' end points of the interval.
#'
#' @param lc Coefficients of the Legendre series
#' @return Vector with the minimum and maximum value of the derivative, in that
#' order
#' @keywords internal
legendrederiv_extrema <- function(lc)
{
  am <- lpolycoef(lc)
  a <- apply(am, 1, sum)

  ## First derivative
  ad1 <- a[2:length(a)]
  ad1 <- ad1 * seq_along(ad1)

  ## second derivative
  ad2 <- ad1[2:length(ad1)]
  ad2 <- ad2 * seq_along(ad2)

  ad2roots <- polyroot(ad2)

  ## We only want the real roots, but there is no harm in evaluating some that have
  ## small but slightly nonzero values.
  ad2roots <- ad2roots[abs(Im(ad2roots)) < 0.1 & Re(ad2roots) > -1 & Re(ad2roots) < 1]
  x <- c(-1,1, Re(ad2roots))

  d1 <- polyeval(x, ad1)

  c(min(d1), max(d1))
}
