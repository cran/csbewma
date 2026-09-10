test_that("known-limits path flags values outside the stated range, matching dissertation Section 6.4.1 encoding (outside = 1)", {
  x <- c(5.0, 6.5, 7.0, 8.5, 9.0)
  result <- dichotomize_range(x, lower_limit = 6.5, upper_limit = 8.5)
  expect_equal(result, c(1L, 0L, 0L, 0L, 1L))
})

test_that("empirical path matches the 1st/99th percentile rule from dissertation Section 6.4.1", {
  set.seed(42)
  reference <- rnorm(10000)
  lower <- stats::quantile(reference, 0.01, names = FALSE)
  upper <- stats::quantile(reference, 0.99, names = FALSE)
  test_values <- c(lower - 1, (lower + upper) / 2, upper + 1)
  result <- dichotomize_range(test_values, reference = reference)
  expect_equal(result, c(1L, 0L, 1L))
})

test_that("empirical path honors custom percentile bounds", {
  set.seed(7)
  reference <- runif(10000)
  result <- dichotomize_range(c(-1, 0.5, 2), reference = reference,
                               lower_percentile = 0.05, upper_percentile = 0.95)
  expect_equal(result, c(1L, 0L, 1L))
})

test_that("errors when both limits and reference are supplied", {
  expect_error(
    dichotomize_range(rnorm(5), lower_limit = 1, upper_limit = 2, reference = rnorm(100)),
    "not both"
  )
})

test_that("errors when neither limits nor reference are supplied", {
  expect_error(dichotomize_range(rnorm(5)), "Supply either")
})

test_that("errors when only one of lower_limit/upper_limit is supplied", {
  expect_error(dichotomize_range(rnorm(5), lower_limit = 1), "Both")
})

test_that("errors when lower_limit is not less than upper_limit", {
  expect_error(dichotomize_range(rnorm(5), lower_limit = 5, upper_limit = 1), "less than")
})

test_that("errors on invalid percentile bounds", {
  expect_error(
    dichotomize_range(rnorm(5), reference = rnorm(100), lower_percentile = 0.9, upper_percentile = 0.1),
    "lower_percentile"
  )
})

test_that("dichotomize_data is untouched by this addition", {
  x <- rnorm(20)
  expect_equal(
    dichotomize_data(x, distribution = "normal", p0 = 0.5),
    as.integer(x > stats::qnorm(0.5))
  )
})

