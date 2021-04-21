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
#' algorithm sensitivity (with 0 being highly selective and 1 being highly sensitive),
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
#' @param itmax Maximum number of iterations in the solver
#' @return Vector of Legendre series coefficients.
#' @export
perfsmooth <- function(x, y, N, maxdegree = 5, xinterval = c(0,1), itmax = 2500)
{
  stopifnot(maxdegree <= 10)
  stopifnot(maxdegree >= 1)

  stopifnot(length(x) == length(y) && length(y) == length(N))

  ## make sure x values are in order from low to high
  iperm <- order(x)
  x <- x[iperm]
  y <- y[iperm]
  N <- N[iperm]

  ## Limit N to 5000 to avoid saturating the likelihood.  With very large N, you
  ## risk getting a lot of likelihood values that are pinned at the largest negative
  ## floating point value.  That makes it hard for the algorithm to find the right
  ## direction.
  N <- pmin(N, 5000)

  ## If the sensitivity scale doesn't go all the way up to the boundaries, add
  ## dummy entries to regularize the behavior near the boundaries.
  if(x[length(x)] < 1) {
    x <- c(x, 1)
    y <- c(y, y[length(y)])
    N <- c(N, 1)
  }
  if(x[1] > 0) {
    x <- c(0, x)
    y <- c(y[1], y)
    N <- c(1, N)
  }

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
    D <- pmax(D, -.Machine$double.xmax)                         # Handle zero densities

    ## Now apply a penalty for nondecreasing values in yhat.
    dmax <- legendrederiv_extrema(lc)[2]

    if(dmax > 0) sum(D) + basepenalty * dmax else sum(D)
  }

  ### Optimize against the optimization target
  ## Initial guess will be a straight line going from the largest y-value to the smallest.
  y0 <- max(y)
  y1 <- min(y)
  a <- 0.5*(y1+y0)
  b <- 0.5*(y1-y0)
  initguess <- rep(0, maxdegree+1)
  initguess[1:2] <- c(a, b)

  ## We could actually supply a gradient function for this optimization, but we'll
  ## wait to see if we really need it before we go down that road.
  rslt <- stats::optim(initguess, opttarg, control=list(maxit=itmax))

  if(rslt$convergence != 0) {
    warning('Convergence failed with code ', rslt$convergence,
            ':  ', rslt$message)
  }

  rslt$par
}
