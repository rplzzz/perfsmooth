test_that("perfsmooth works", {
  ## simulate data for use in the algorithm
  x <- seq(0.01,1, 0.1)
  N <- round(250*x)
  p0 <- -0.6 * x^2 + 0.8
  set.seed(867-5309)
  p <- rbinom(length(N), N, p0) / N

  lcrslt <- expect_silent(perfsmooth(x, p, N))

  xtst <- seq(0.01, 1, 0.01)
  ytst <- lserieseval(xtst, lcrslt, c(0,1))
  ## Test to see that f(x) is monotonic
  expect_true(all(diff(ytst) < 0))

  ## Test to see that f(x) is at least reasonably close to the intended values
  alpha <- p0*N + 1
  beta <- N+2-alpha
  ytst2 <- lserieseval(xtst, lcrslt, c(0,1))
  qb <- qbeta(ytst2, alpha, beta)
  lim <- 1/(2*length(x))
  expect_true(all(qb > lim) & all(qb < 1-lim))
})

test_that('perfsmooth handles large N', {
  ## simulate data for use in the algorithm
  x <- seq(0.01,1, 0.1)
  N <- round(250000*x)
  p0 <- -0.6 * x^2 + 0.8
  set.seed(867-5309)
  p <- rbinom(length(N), N, p0) / N

  lcrslt <- expect_silent(perfsmooth(x, p, N))

  xtst <- seq(0.01, 1, 0.01)
  ytst <- lserieseval(xtst, lcrslt, c(0,1))
  ## Test to see that f(x) is monotonic
  expect_true(all(diff(ytst) < 0))

  ## Test to see that f(x) is at least reasonably close to the intended values
  alpha <- p0*N + 1
  beta <- N+2-alpha
  ytst2 <- lserieseval(xtst, lcrslt, c(0,1))
  qb <- qbeta(ytst2, alpha, beta)
  lim <- 1/(2*length(x))
  expect_true(all(qb > lim) & all(qb < 1-lim))
})
