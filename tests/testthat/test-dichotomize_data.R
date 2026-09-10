test_that("normal distribution thresholds at qnorm(p0) as documented", {
  x <- c(-2, -1, -0.0001, 0, 0.0001, 1, 2)
  result <- dichotomize_data(x, distribution = "normal", p0 = 0.5)
  expect_equal(result, as.integer(x > 0))
})

test_that("laplace distribution reproduces the exact closed-form threshold", {
  p0 <- 0.3
  expected_threshold <- log(2 * p0)
  x <- c(expected_threshold - 0.5, expected_threshold + 0.5)
  result <- dichotomize_data(x, distribution = "laplace", p0 = p0)
  expect_equal(result, c(0L, 1L))
})

test_that("uniform distribution thresholds at p0 itself", {
  x <- c(0.1, 0.5, 0.9)
  result <- dichotomize_data(x, distribution = "uniform", p0 = 0.5)
  expect_equal(result, c(0L, 0L, 1L))
})

test_that("exponential distribution reproduces the qexp threshold", {
  p0 <- 0.4
  th <- stats::qexp(p0, rate = 1)
  x <- c(th - 0.1, th + 0.1)
  result <- dichotomize_data(x, distribution = "exponential", p0 = p0)
  expect_equal(result, c(0L, 1L))
})

test_that("errors on an unrecognized distribution", {
  expect_error(
    dichotomize_data(rnorm(10), distribution = "not_a_distribution"),
    "Unknown distribution"
  )
})

test_that("the comparison direction actually matters (would fail if reversed)", {
  x <- c(-1, 1)
  correct <- dichotomize_data(x, distribution = "normal", p0 = 0.5)
  deliberately_reversed <- as.integer(x < stats::qnorm(0.5))
  expect_false(isTRUE(all.equal(correct, deliberately_reversed)))
})

