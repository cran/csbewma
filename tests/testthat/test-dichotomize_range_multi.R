test_that("applies different known limits per column correctly", {
  df <- data.frame(
    a = c(1, 5, 10, 15),
    b = c(100, 150, 200, 250)
  )
  specs <- list(
    a = list(lower_limit = 2, upper_limit = 12),
    b = list(lower_limit = 120, upper_limit = 220)
  )
  result <- dichotomize_range_multi(df, specs)
  expect_equal(result[, "a"], c(1L, 0L, 0L, 1L))
  expect_equal(result[, "b"], c(1L, 0L, 0L, 1L))
})

test_that("different columns can use different methods (limits vs reference)", {
  set.seed(42)
  ref_b <- rnorm(1000, mean = 100, sd = 5)
  df <- data.frame(
    a = c(1, 20),
    b = c(100, 200)
  )
  specs <- list(
    a = list(lower_limit = 0, upper_limit = 10),
    b = list(reference = ref_b, lower_percentile = 0.01, upper_percentile = 0.99)
  )
  result <- dichotomize_range_multi(df, specs)
  expect_equal(result[, "a"], c(0L, 1L))
  expected_b <- as.integer(
    df$b < stats::quantile(ref_b, 0.01, names = FALSE) |
      df$b > stats::quantile(ref_b, 0.99, names = FALSE)
  )
  expect_equal(result[, "b"], expected_b)
})

test_that("errors when a column's specification is missing", {
  df <- data.frame(a = 1:3, b = 4:6)
  specs <- list(a = list(lower_limit = 0, upper_limit = 10))
  expect_error(dichotomize_range_multi(df, specs), "No range specification")
})

test_that("errors when specs has an entry with no matching column", {
  df <- data.frame(a = 1:3)
  specs <- list(
    a = list(lower_limit = 0, upper_limit = 10),
    z = list(lower_limit = 0, upper_limit = 1)
  )
  expect_error(dichotomize_range_multi(df, specs), "no matching column")
})

test_that("preserves column names and dimensions", {
  df <- data.frame(a = 1:5, b = 6:10)
  specs <- list(
    a = list(lower_limit = 0, upper_limit = 10),
    b = list(lower_limit = 0, upper_limit = 20)
  )
  result <- dichotomize_range_multi(df, specs)
  expect_equal(dim(result), dim(df))
  expect_equal(colnames(result), colnames(df))
})

test_that("the direction of the range check actually matters (would fail if reversed)", {
  df <- data.frame(a = c(1, 50))
  specs <- list(a = list(lower_limit = 0, upper_limit = 10))
  correct <- dichotomize_range_multi(df, specs)
  deliberately_reversed <- as.integer(!as.logical(correct[, "a"]))
  expect_false(isTRUE(all.equal(as.integer(correct[, "a"]), deliberately_reversed)))
})
