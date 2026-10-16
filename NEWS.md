# csbewma 1.2.0

## New features

* `p0` (the in-control proportion) can now be a numeric vector of length
  `k` (one value per stream) instead of a single shared value, in
  `run_csb_ewma()`, `csb_ewma()`, `calculate_pvalues()`, `identify_ooc()`,
  and `dichotomize_data()`. This supports monitoring multiple streams that
  do not all share the same in-control proportion (for example, several
  water-quality variables whose historical out-of-band rates differ).
  * In `run_csb_ewma()`/`csb_ewma()`, the chart's mean and variance become
    `mu0 = sum(p0)` and `sigma2_0 = sum(p0 * (1 - p0))` across streams. A
    scalar `p0` is recycled to length `k` before summing, which reproduces
    the previous `k * p0` / `k * p0 * (1 - p0)` formulas exactly, so this
    is a strict, backward-compatible generalization: every existing call
    with a scalar `p0` returns bit-identical results.
  * In `calculate_pvalues()`/`identify_ooc()`, each stream's exact binomial
    p-value is computed against its own `p0[i]` instead of one shared
    value, again falling back to the original single-value behavior when
    `p0` has length 1.
  * In `dichotomize_data()`, the in-control threshold `th` becomes a vector
    the same length as `data` when `p0` is a vector, so a single call can
    dichotomize one time slice across streams with different in-control
    proportions. The Laplace branch was rewritten with `ifelse()` so it
    vectorizes correctly over a vector `p0` (it previously used a scalar
    `if`/`else`, which is only valid for a single value).
  * The exact time-varying variance formula (`var_rt_exact_single()` /
    `precompute_variance()`, Theorem 2 of the derivations paper) is
    unchanged: it depends only on `lambda` and `t`, never on `p0` or `k`,
    so it already applies to the vector-`p0` case with no modification.
  * All four functions now validate that a vector `p0` has the expected
    length (`k` for the chart/post-hoc functions, `length(data)` for
    `dichotomize_data()`) and error clearly otherwise.

# csbewma 1.1.0

## New features

* Added `dichotomize_range()`: converts continuous measurements to a binary
  indicator based on whether each value falls inside or outside an
  acceptable operating range, for practical real-data use as described in
  dissertation Section 6.4.1. The range can be supplied directly as known
  `lower_limit`/`upper_limit` values, or estimated from historical
  `reference` data using percentiles (default 1st/99th).
* Added `dichotomize_range_multi()`: applies `dichotomize_range()`
  separately to each variable in a dataset, since different variables
  monitored together typically have different acceptable ranges (for
  example, different variables in the same process rarely share one
  tolerance band). Every variable's range must be specified explicitly;
  there is no shared default across variables.
* Added `check_stream_correlation()`: reports the pairwise Pearson
  correlation among candidate streams and flags any pair whose absolute
  correlation exceeds a threshold (default 0.3, per dissertation Section
  6.2.7), as a diagnostic check of the independence assumption underlying
  the CSB-EWMA chart. This function is diagnostic only; it does not remove
  or adjust correlated streams, since no validated correction method for
  correlated streams currently exists for this chart.
* Added a `stop_at_signal` argument to `run_csb_ewma()` and `csb_ewma()`.
  The default (`stop_at_signal = TRUE`) is unchanged from every released
  version through 1.1.0: monitoring stops at the first signal. Setting
  `stop_at_signal = FALSE` instead continues monitoring through `max_time`,
  recording every signal in `signal_times` without resetting the cumulative
  sum or EWMA statistic after a signal. `csb_ewma()` performs post-hoc
  identification separately at each signal in this mode, stored as a named
  list in `flagged_by_signal` (`flagged` is kept as the first signal's
  result for backward compatibility). `plot.csb_ewma()` and
  `plot_csb_ewma_direct()` annotate every signal in `signal_times` when
  present. This continuation mode is a practical convenience for
  real-time or dashboard-style monitoring; it is not itself a procedure
  described in the dissertation, which presents only the
  stop-at-first-signal convention.

* Added a `verbose` argument to `csb_ewma()` (default `FALSE`). When `TRUE`, prints informational messages about whether the input data required dichotomization; these are silent by default.

## Bug fixes

* Fixed a label placement bug in `plot.csb_ewma()` and
  `plot_csb_ewma_direct()`: the "Signal at t = ..." annotation used
  `hjust = 0` instead of `hjust = 1`, which misaligned the label against
  the signal line.

* Fixed a variance-calculation bug in `var_rt_exact_single()`: the loop
  `for (j in 1:(t-1))` evaluated to `1:0` when `t = 1`, which R treats as
  the two-element sequence `c(1, 0)` rather than an empty range, so two
  spurious terms were added to `Var(r_1)`, inflating it to about 4.8 times
  its correct value and widening `UCL_1`/`LCL_1` accordingly. Replaced with
  `seq_len(t - 1)`, which is empty when `t = 1` as required; `Var(r_1)`
  now correctly equals `lambda^2`. All other time points are unaffected.

* Fixed a documentation bug: a stray line in the roxygen `@examples`
  blocks for `csb_ewma()` and `rlaplace()` was being included in the
  generated Rd files' `\examples{}` sections, causing R CMD check's
  example-parsing step to fail with a syntax error. Removed.

## Documentation

* Added `URL` and `BugReports` fields to `DESCRIPTION`, pointing to the
  package's development repository.
* Added `@seealso` cross-references between `run_csb_ewma()` and
  `csb_ewma()` clarifying that they are not duplicates: `run_csb_ewma()`
  is the lower-level monitoring engine, useful directly when a variance
  cache should be precomputed once and reused across many chart runs (for
  example, in simulation studies); `csb_ewma()` is the higher-level
  convenience wrapper that dichotomizes continuous data, precomputes
  variance, calls `run_csb_ewma()`, and performs post-hoc identification
  automatically.
* Added unit tests (`testthat`) covering every exported function.

# csbewma 1.0.1

* First version accepted to CRAN.

# csbewma 1.0.0

* Initial development version.
