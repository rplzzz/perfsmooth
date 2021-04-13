test_that("Legendre series eval works on default interval", {
  ## Linear combination of P1, P2, and P5
  lc <- c(0, -1, -0.25, 0, 0, -0.1, rep(0,5))
  x <- seq(-1,1,0.25)
  ## calculate expected values manually
  p1 <- -x
  p2 <- -0.25 * (3*x^2 - 1)/2
  p5 <- -0.1 * (63*x^5 - 70*x^3 + 15*x)/8
  expct <- p1 + p2 + p5
  expect_equal(lserieseval(x, lc), expct)
})

test_that('Legendre series eval works on rescaled interval', {
  lc <- c(0, -1, -0.25, 0, 0, -0.1, rep(0,5))
  x <- seq(-1,1,0.25)
  xi <- seq(1,2,0.125)

  expect_equal(lserieseval(xi, lc, c(1,2)), lserieseval(x, lc))

})

test_that('Truncated coefficients are filled out', {
  lct <- c(0, -1, -0.25, 0, 0, -0.1)
  lc <- c(lct, rep(0,5))
  x <- seq(-1,1,0.25)

  expect_equal(lserieseval(x,lct), lserieseval(x,lc))
})

test_that('Excess coefficients generate a warning', {
  lc <- c(0, -1, -0.25, 0, 0, -0.1, rep(0,5))
  lce <- c(lc, 0.1)
  x <- seq(-1,1,0.25)

  expect_warning({y <- lserieseval(x, lce)}, regexp = 'Excess will be dropped')
  expect_equal(y, lserieseval(x,lc))
})
