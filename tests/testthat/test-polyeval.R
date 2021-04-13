test_that('polynomials are evaluated correctly', {
  a <- c(1,-3,2)           # 2x^2 - 3x + 1
  x <- seq(-1,1,0.5)
  expect_equal(polyeval(x,a), c(6,3,1,0,0))
})

test_that('polynomial sums are evaluated correctly', {
  a <- c(1,-3,2)
  b <- c(-1,1,0)
  am <- cbind(b,a)
  x <- seq(-1,1,0.5)
  expect_equal(polysumeval(x,am), c(4,1.5,0,-0.5,0))
})
