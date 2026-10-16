test_that("run_csb_ewma: scalar p0 is unchanged (backward compatibility)", {
  set.seed(1)
  bin_matrix <- matrix(rbinom(5 * 50, 1, 0.5), nrow = 5, ncol = 50)
  var_cache <- precompute_variance(0.2, max_t = 50)
  res_scalar <- run_csb_ewma(bin_matrix, lambda = 0.2, L = 1.5, var_cache, p0 = 0.5)

  k <- 5
  mu0 <- k * 0.5
  sigma2_0 <- k * 0.5 * (1 - 0.5)
  cum_sum <- 0
  r_prev <- 0
  r_manual <- numeric(50)
  for (t in 1:50) {
    C_t <- sum(bin_matrix[, t])
    cum_sum <- cum_sum + C_t
    W_t <- (cum_sum - mu0 * t) / sqrt(t * sigma2_0)
    r_prev <- 0.2 * W_t + (1 - 0.2) * r_prev
    r_manual[t] <- r_prev
  }
  expect_equal(res_scalar$r_history, r_manual[seq_along(res_scalar$r_history)])
})

test_that("run_csb_ewma: vector p0 of length k matches the sum-based formula", {
  set.seed(2)
  k <- 4
  p0_vec <- c(0.5, 0.4, 0.6, 0.5)
  bin_matrix <- matrix(rbinom(k * 40, 1, 0.5), nrow = k, ncol = 40)
  var_cache <- precompute_variance(0.2, max_t = 40)
  res <- run_csb_ewma(bin_matrix, lambda = 0.2, L = 1.5, var_cache, p0 = p0_vec)

  mu0 <- sum(p0_vec)
  sigma2_0 <- sum(p0_vec * (1 - p0_vec))
  cum_sum <- 0
  r_prev <- 0
  r_manual <- numeric(40)
  for (t in 1:40) {
    C_t <- sum(bin_matrix[, t])
    cum_sum <- cum_sum + C_t
    W_t <- (cum_sum - mu0 * t) / sqrt(t * sigma2_0)
    r_prev <- 0.2 * W_t + (1 - 0.2) * r_prev
    r_manual[t] <- r_prev
  }
  expect_equal(res$r_history, r_manual[seq_along(res$r_history)])
})

test_that("run_csb_ewma: a vector p0 of the wrong length errors clearly", {
  bin_matrix <- matrix(rbinom(3 * 20, 1, 0.5), nrow = 3, ncol = 20)
  var_cache <- precompute_variance(0.2, max_t = 20)
  expect_error(
    run_csb_ewma(bin_matrix, lambda = 0.2, L = 1.5, var_cache, p0 = c(0.5, 0.5)),
    "p0 must be a single value or a vector of length k"
  )
})

test_that("calculate_pvalues: scalar p0 is unchanged (backward compatibility)", {
  set.seed(3)
  bin_mat <- matrix(rbinom(200, 1, 0.5), nrow = 4, ncol = 50)
  pvals <- calculate_pvalues(bin_mat, p0 = 0.5)
  successes <- rowSums(bin_mat)
  expected <- sapply(seq_len(4), function(i) {
    1 - stats::pbinom(successes[i] - 1, size = 50, prob = 0.5)
  })
  expect_equal(pvals, expected)
})

test_that("calculate_pvalues: vector p0 tests each stream against its own proportion", {
  set.seed(4)
  bin_mat <- matrix(rbinom(200, 1, 0.5), nrow = 4, ncol = 50)
  p0_vec <- c(0.5, 0.4, 0.6, 0.3)
  pvals <- calculate_pvalues(bin_mat, p0 = p0_vec)
  successes <- rowSums(bin_mat)
  expected <- sapply(seq_len(4), function(i) {
    1 - stats::pbinom(successes[i] - 1, size = 50, prob = p0_vec[i])
  })
  expect_equal(pvals, expected)
})

test_that("calculate_pvalues: a vector p0 of the wrong length errors clearly", {
  bin_mat <- matrix(rbinom(200, 1, 0.5), nrow = 4, ncol = 50)
  expect_error(
    calculate_pvalues(bin_mat, p0 = c(0.5, 0.5)),
    "p0 must be a single value or a vector of length k"
  )
})

test_that("identify_ooc: vector p0 is forwarded to calculate_pvalues", {
  set.seed(5)
  bin_mat <- matrix(rbinom(400, 1, 0.5), nrow = 4, ncol = 100)
  p0_vec <- c(0.5, 0.4, 0.6, 0.3)
  result <- identify_ooc(bin_mat, p0 = p0_vec)
  direct_pvals <- calculate_pvalues(bin_mat, p0 = p0_vec)
  expect_equal(sort(result$p_value), sort(direct_pvals))
})

test_that("dichotomize_data: scalar p0 is unchanged (backward compatibility)", {
  x <- c(-2, -1, -0.0001, 0, 0.0001, 1, 2)
  result <- dichotomize_data(x, distribution = "normal", p0 = 0.5)
  expect_equal(result, as.integer(x > 0))
})

test_that("dichotomize_data: vector p0 gives one threshold per observation (normal)", {
  x <- c(-1, 0, 1, 5)
  p0_vec <- c(0.5, 0.5, 0.1, 0.99)
  result <- dichotomize_data(x, distribution = "normal", p0 = p0_vec)
  expected <- as.integer(x > stats::qnorm(p0_vec))
  expect_equal(result, expected)
})

test_that("dichotomize_data: vector p0 vectorizes correctly for laplace", {
  # Regression test: the laplace branch previously used a scalar if/else on
  # p0, which is only valid for length-1 p0. It was rewritten with
  # ifelse() to vectorize correctly; this test would fail against the old
  # scalar if/else when given a vector p0.
  x <- c(-1, 0, 1, 2)
  p0_vec <- c(0.2, 0.5, 0.7, 0.9)
  result <- dichotomize_data(x, distribution = "laplace", p0 = p0_vec)
  expected_th <- ifelse(p0_vec < 0.5, log(2 * p0_vec), -log(2 * (1 - p0_vec)))
  expect_equal(result, as.integer(x > expected_th))
})

test_that("dichotomize_data: a vector p0 of the wrong length errors clearly", {
  expect_error(
    dichotomize_data(c(1, 2, 3), distribution = "normal", p0 = c(0.5, 0.5)),
    "p0 must be a single value or a vector the same length as data"
  )
})
