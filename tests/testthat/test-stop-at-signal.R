test_that("run_csb_ewma stop_at_signal defaults to TRUE and reproduces the exact prior single-signal behavior", {
  set.seed(101)
  bin_matrix <- matrix(rbinom(6 * 60, 1, 0.5), nrow = 6, ncol = 60)
  for (i in 1:2) bin_matrix[i, ] <- rbinom(60, 1, 0.85)
  lambda <- 0.2
  L <- 1.3
  p0 <- 0.5
  var_cache <- precompute_variance(lambda, max_t = 60)

  result_default <- run_csb_ewma(bin_matrix, lambda, L, var_cache, p0 = p0)
  result_explicit <- run_csb_ewma(bin_matrix, lambda, L, var_cache, p0 = p0, stop_at_signal = TRUE)

  expect_identical(result_default, result_explicit)
  expect_true(result_default$stop_at_signal)
  expect_true(result_default$signal_detected)
  expect_equal(result_default$T_sig, result_default$signal_time)
  expect_equal(result_default$signal_times, result_default$signal_time)
  expect_length(result_default$signal_times, 1)
})

test_that("run_csb_ewma stop_at_signal = FALSE matches a hand-computed recursion and does not reset the EWMA statistic", {
  # A small, fully deterministic (non-random) 3-stream, 6-time-point series so
  # every r_t can be independently re-derived from the recursion documented
  # for run_csb_ewma(): cum_sum_t = cum_sum_(t-1) + C_t,
  # W_t = (cum_sum_t - mu0 * t) / sqrt(t * sigma2_0),
  # r_t = lambda * W_t + (1 - lambda) * r_(t-1), with cum_sum_0 = 0, r_0 = 0.
  bin_matrix <- matrix(c(
    1, 1, 1,
    1, 1, 1,
    1, 1, 0,
    1, 1, 1,
    0, 1, 1,
    1, 1, 1
  ), nrow = 3, ncol = 6)

  lambda <- 0.5
  p0 <- 0.5
  k <- nrow(bin_matrix)
  n <- ncol(bin_matrix)

  C <- colSums(bin_matrix)
  mu0 <- k * p0
  sigma2_0 <- k * p0 * (1 - p0)
  cum_sum <- cumsum(C)
  W <- (cum_sum - mu0 * seq_len(n)) / sqrt(seq_len(n) * sigma2_0)
  expected_r <- Reduce(function(r_prev, w) lambda * w + (1 - lambda) * r_prev,
                        W, accumulate = TRUE, init = 0)[-1]

  var_cache <- precompute_variance(lambda, max_t = n)

  # A small L so at least one signal fires partway through; continuation is
  # then verified by comparing r_history across the full series.
  L <- 0.6

  result_continue <- run_csb_ewma(bin_matrix, lambda, L, var_cache,
                                   max_time = n, p0 = p0, stop_at_signal = FALSE)

  expect_equal(result_continue$r_history, expected_r, tolerance = 1e-8)
  expect_true(result_continue$signal_detected)
  expect_true(length(result_continue$signal_times) >= 1)
  expect_equal(result_continue$T_sig, n)
  expect_false(result_continue$stop_at_signal)

  # Compare against stop_at_signal = TRUE on the identical inputs: the
  # trajectory up through the first signal must be identical, which is only
  # possible if continuation mode never resets cum_sum or r_prev.
  result_stop <- run_csb_ewma(bin_matrix, lambda, L, var_cache,
                               max_time = n, p0 = p0, stop_at_signal = TRUE)
  first_signal <- result_stop$T_sig
  expect_equal(result_continue$r_history[1:first_signal],
               result_stop$r_history[1:first_signal])
  expect_equal(result_continue$UCL_history[1:first_signal],
               result_stop$UCL_history[1:first_signal])
  expect_equal(result_continue$LCL_history[1:first_signal],
               result_stop$LCL_history[1:first_signal])
})

test_that("run_csb_ewma signal_times records every signal, not only the first (would fail under the old single-signal behavior)", {
  set.seed(202)
  bin_matrix <- matrix(rbinom(5 * 80, 1, 0.5), nrow = 5, ncol = 80)
  for (i in 1:3) bin_matrix[i, ] <- rbinom(80, 1, 0.8)
  lambda <- 0.2
  L <- 1.2
  p0 <- 0.5
  var_cache <- precompute_variance(lambda, max_t = 80)

  result <- run_csb_ewma(bin_matrix, lambda, L, var_cache, p0 = p0, stop_at_signal = FALSE)

  expect_true(length(result$signal_times) > 1)
  expect_true(all(result$signal_times >= 1 & result$signal_times <= 80))
  expect_true(all(diff(result$signal_times) > 0))
  expect_equal(result$signal_times[1], result$signal_time)
})

test_that("run_csb_ewma stop_at_signal = FALSE with no signal reproduces the no-signal warning and empty signal_times", {
  set.seed(303)
  bin_matrix <- matrix(rbinom(4 * 30, 1, 0.5), nrow = 4, ncol = 30)
  lambda <- 0.2
  L <- 8
  p0 <- 0.5
  var_cache <- precompute_variance(lambda, max_t = 30)

  expect_warning(
    result <- run_csb_ewma(bin_matrix, lambda, L, var_cache, p0 = p0, stop_at_signal = FALSE),
    "No signal detected"
  )
  expect_false(result$signal_detected)
  expect_equal(result$T_sig, 30)
  expect_equal(result$signal_times, integer(0))
})

test_that("csb_ewma stop_at_signal = FALSE performs post-hoc identification separately at each signal, and flagged stays backward compatible", {
  set.seed(404)
  bin_data <- matrix(rbinom(8 * 100, 1, 0.5), nrow = 8, ncol = 100)
  for (i in 1:3) bin_data[i, ] <- rbinom(100, 1, 0.85)

  result <- csb_ewma(bin_data, lambda = 0.2, L = 1.2, stop_at_signal = FALSE)

  expect_true(length(result$signal_times) >= 1)
  expect_identical(names(result$flagged_by_signal), as.character(result$signal_times))
  expect_identical(result$flagged, result$flagged_by_signal[[1]])
})

test_that("csb_ewma default matches explicit stop_at_signal = TRUE (backward compatible)", {
  set.seed(505)
  bin_data <- matrix(rbinom(8 * 60, 1, 0.5), nrow = 8, ncol = 60)
  for (i in 1:3) bin_data[i, ] <- rbinom(60, 1, 0.85)

  result_default <- csb_ewma(bin_data, lambda = 0.2, L = 1.2)
  result_explicit <- csb_ewma(bin_data, lambda = 0.2, L = 1.2, stop_at_signal = TRUE)

  expect_identical(result_default, result_explicit)
})
