#' Find an optimal monotonic performance curve.
#'
#' Find a curve that best matches a set of points, subject to the constraint
#' that the curve must be monotonic.
#'
#' The points are assumed to be observations of a binary event, which give an
#' estimate of the probability of the event.  The uncertainty of those probability
#' values is therefore beta distributed with parameters derived from the number
#' of observations.
#'
#' The intended use here is for a case where the x-values are a measure of
#' algorithm selectivity (with 0 being highly selective and 1 being highly permissive),
#' and the y-values are PPV (or similar measure).  A more selective algorithm should
#' generate fewer false positives, but with fewer cases, there is more uncertainty
#' in the measured value.
#'
#' The smoothed curve is represented as a series of Legendre polynomials, up to
#' degree 10 (or less, if specified).  The return value is the vector of series coefficients.
#'
#' @param x x-values of control points.
#' @param y y-values of control points.  These should be fractions between 0 and 1.
#' @param N Number of cases that contributed to each y value.
#' @param maxdegree Maximum degree to use in the Legendre series.  Must be <= 10.
#' @param xinterval Domain of valid values for \code{x} (whether or not the entire
#' interval is actually represented in the x values).  Default is [0,1].
#' @return Vector of Legendre series coefficients.
#' @export
perfsmooth <- function(x, y, N, maxdegree = 5, xinterval = c(0,1))
{
  stopifnot(maxdegree <= 10)
  stopifnot(maxdegree >= 1)

  stopifnot(length(x) == length(y) && length(y) == length(N))

  ## parameters of the beta distributions
  alpha <- 1 + y*N
  beta <- N+2 - alpha

  ## base size of the penalty for having a nondecreasing segment.  This will be
  ## multiplied by the max of the derivative on the (shifted) interval [-1,1], if
  ## that maximum is > 0.  This construct gives us a value that is sufficiently
  ## large that it will override any agreement with the data, ensuring that
  ## nonmonotonic functions are excluded.
  basepenalty <- -2*sum(stats::dbeta(0.9, 1, N, log=TRUE))

  ## optimization target function is a function of the Legendre series coefficients.
  opttarg <- function(lc) {
    yhat <- lserieseval(x, lc, xinterval)
    D <- -2*stats::dbeta(yhat, alpha, beta, log = TRUE)         # Deviance = -2*log(L)

    ## Now apply a penalty for nondecreasing values in yhat.
    dmax <- legendrederiv_extrema(lc)[2]

    if(dmax > 0) sum(D) + basepenalty * dmax else sum(D)
  }

  ## Optimize against the optimization target
  initguess <- rep(0, maxdegree+1)
  initguess[1:2] <- c(0.5, -0.5)         # initial guess is straight line with slope = -1

  ## We could actually supply a gradient function for this optimization, but we'll
  ## wait to see if we really need it before we go down that road.
  rslt <- optim(initguess, opttarg, control=list(maxit=2500))

  if(rslt$convergence != 0) {
    warning('Convergence failed with code ', rslt$convergence,
            ':  ', rslt$message)
  }

  rslt$par
}
