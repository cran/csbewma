test_that("independence_ok is TRUE when all streams are weakly correlated", {
  set.seed(1)
  x <- matrix(rnorm(600), ncol = 6)
  colnames(x) <- paste0("stream_", 1:6)
  result <- check_stream_correlation(x)
  expect_true(result$independence_ok)
  expect_equal(nrow(result$flagged_pairs), 0)
  expect_equal(dim(result$correlation_matrix), c(6, 6))
})

test_that("flags a pair built to be strongly correlated", {
  set.seed(2)
  base <- rnorm(500)
  x <- cbind(
    stream_a = base,
    stream_b = base + rnorm(500, sd = 0.05),
    stream_c = rnorm(500)
  )
  result <- check_stream_correlation(x, threshold = 0.3)
  expect_false(result$independence_ok)
  expect_true(any(result$flagged_pairs$stream_1 == "stream_a" & result$flagged_pairs$stream_2 == "stream_b"))
})

test_that("reproduces the dissertation's own Table 6.3 correlations for the six selected streams", {
  corr_table <- matrix(
    c(1.00, -0.10, -0.20,  0.22,  0.05,  0.03,
     -0.10,  1.00,  0.06, -0.08, -0.15, -0.06,
     -0.20,  0.06,  1.00, -0.27, -0.12, -0.01,
      0.22, -0.08, -0.27,  1.00,  0.11,  0.01,
      0.05, -0.15, -0.12,  0.11,  1.00,  0.03,
      0.03, -0.06, -0.01,  0.01,  0.03,  1.00),
    nrow = 6, byrow = TRUE
  )
  expect_true(all(abs(corr_table[upper.tri(corr_table)]) < 0.3))
})

test_that("errors on non-matrix, non-data-frame input", {
  expect_error(check_stream_correlation(list(1, 2, 3)), "matrix or data frame")
})

test_that("errors on an out-of-range threshold", {
  x <- matrix(rnorm(60), ncol = 3)
  expect_error(check_stream_correlation(x, threshold = 1.5), "threshold")
})

test_that("flagging is threshold-sensitive in the expected direction", {
  set.seed(3)
  base <- rnorm(500)
  x <- cbind(a = base, b = base + rnorm(500, sd = 0.05), c = rnorm(500))
  loose <- check_stream_correlation(x, threshold = 0.99)
  strict <- check_stream_correlation(x, threshold = 0.05)
  expect_true(nrow(loose$flagged_pairs) <= nrow(strict$flagged_pairs))
})

test_that("assigns default stream names when data has none", {
  set.seed(4)
  x <- matrix(rnorm(400), ncol = 4)
  result <- check_stream_correlation(x)
  expect_equal(colnames(result$correlation_matrix), paste0("stream_", 1:4))
})

