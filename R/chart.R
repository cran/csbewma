#' Run CSB-EWMA Chart on Binary Data
#'
#' This function implements the core CSB-EWMA monitoring algorithm exactly as
#' implemented in the original simulation code.
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Initialize cumulative sum (cum_sum = 0) and EWMA (r_prev = 0)
#'   \item For each time point t = 1, 2, ..., max_time:
#'   \itemize{
#'     \item Get binary vector for current time point
#'     \item Calculate C_t = sum of binary indicators
#'     \item Update cumulative sum: cum_sum = cum_sum + C_t
#'     \item Compute standardized statistic: W_t = (cum_sum - mu0*t) / sqrt(t*sigma2_0)
#'     \item Update EWMA: r_t = lambda * W_t + (1 - lambda) * r_prev
#'     \item Get exact variance from precomputed cache: v_t = var_cache[t]
#'     \item Compute control limits: UCL_t = L * sqrt(v_t), LCL_t = -L * sqrt(v_t)
#'     \item If r_t > UCL_t or r_t < LCL_t, a signal is recorded at time t
#'     \item Update r_prev = r_t and continue, unless stop_at_signal ends monitoring here
#'   }
#' }
#'
#' By default (\code{stop_at_signal = TRUE}), monitoring stops at the first
#' signal, matching the behavior of every released version of this function
#' through 1.1.0. Setting \code{stop_at_signal = FALSE} instead continues
#' monitoring through \code{max_time}, recording every time point at which a
#' signal condition is met in \code{signal_times}, without resetting the
#' cumulative sum or the EWMA statistic after a signal; the formula above is
#' unchanged, only the decision to stop early changes. This continuation mode
#' is a practical convenience for real-time or dashboard-style monitoring
#' where the process keeps running after an alarm; it is not itself a
#' procedure described in the dissertation, which presents only the
#' stop-at-first-signal convention.
#'
#' @param bin_matrix A matrix of binary indicators (streams as rows, time as columns)
#' @param lambda Smoothing parameter for EWMA (0 < lambda <= 1)
#' @param L Control limit multiplier
#' @param var_cache Precomputed variance vector from precompute_variance()
#' @param max_time Maximum time points to monitor (default = NULL uses all)
#' @param p0 In-control proportion (default = 0.5). May instead be a numeric
#'   vector of length k (one value per stream/row of bin_matrix) to allow
#'   each stream to have its own in-control proportion; mu0 and sigma2_0 are
#'   then computed as sum(p0) and sum(p0 * (1 - p0)) across streams, which
#'   reduces exactly to the scalar formulas when p0 has length 1.
#' @param stop_at_signal Logical. If TRUE (default), monitoring stops at the
#'   first signal, reproducing the behavior of every version through 1.1.0.
#'   If FALSE, monitoring continues through max_time and every signal time is
#'   recorded in signal_times.
#' @return A list of class "csb_ewma" containing chart results. Includes
#'   signal_time and signal_detected (the first signal, kept for backward
#'   compatibility) and signal_times (an integer vector of every time point
#'   at which a signal was recorded; empty if none).
#' @seealso [csb_ewma()] for the higher-level convenience wrapper that dichotomizes continuous data, precomputes variance, calls this function, and performs post-hoc identification automatically.
#' @export
#'
#' @examples
#' bin_matrix <- matrix(rbinom(10*200, 1, 0.5), nrow = 10, ncol = 200)
#' for(i in 1:3) bin_matrix[i, ] <- rbinom(100, 1, 0.8)
#' var_cache <- precompute_variance(0.175, max_t = 200)
#' result <- run_csb_ewma(bin_matrix, lambda = 0.175, L = 1.375, var_cache)
#' print(paste("Signal at time:", result$signal_time))
#'
#' # Continue monitoring past the first signal
#' result2 <- run_csb_ewma(bin_matrix, lambda = 0.175, L = 1.375, var_cache,
#'                          stop_at_signal = FALSE)
#' print(result2$signal_times)
#'
#' # Heterogeneous streams: one in-control proportion per stream
#' p0_vec <- c(0.5, 0.5, 0.5, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5)
#' result3 <- run_csb_ewma(bin_matrix, lambda = 0.175, L = 1.375, var_cache,
#'                          p0 = p0_vec)
#' print(paste("Signal at time:", result3$signal_time))
# ----------------------------------------------------------------------------
run_csb_ewma <- function(bin_matrix, lambda, L, var_cache,
                         max_time = NULL, p0 = 0.5, stop_at_signal = TRUE) {

  # ========================================================================
  # STEP 1: Input Validation
  # ========================================================================

  # Check that bin_matrix is a matrix
  if (!is.matrix(bin_matrix)) {
    stop("bin_matrix must be a matrix")
  }

  # Check that bin_matrix contains only 0s and 1s
  if (!all(bin_matrix %in% c(0, 1))) {
    stop("bin_matrix must contain only 0s and 1s")
  }

  # Check that lambda is between 0 and 1
  if (lambda <= 0 || lambda > 1) {
    stop("lambda must be between 0 and 1")
  }

  # Check that L is positive
  if (L <= 0) {
    stop("L must be positive")
  }

  # ========================================================================
  # STEP 2: Get Dimensions and Set Parameters
  # ========================================================================

  # Get number of streams (rows) and total time points (columns)
  k <- nrow(bin_matrix)
  n_time <- ncol(bin_matrix)

  # Set max_time to total available if not specified
  if (is.null(max_time)) {
    max_time <- n_time
  } else {
    max_time <- min(max_time, n_time)
  }

  # Check that var_cache has enough length
  if (length(var_cache) < max_time) {
    stop("var_cache length is less than max_time")
  }

  # Compute derived parameters.
  # p0 may be a single in-control proportion shared by every stream (the
  # original, scalar behavior), or a numeric vector of length k giving a
  # separate in-control proportion for each stream (row of bin_matrix).
  # A scalar p0 is recycled to length k before summing, which reproduces
  # the original k * p0 / k * p0 * (1 - p0) formulas exactly. This makes
  # the vector case a strict, backward-compatible generalization:
  #   mu0     = sum(p0_i)
  #   sigma2_0 = sum(p0_i * (1 - p0_i))
  # which reduce to k*p0 and k*p0*(1-p0) when every p0_i equals p0.
  if (!(length(p0) == 1 || length(p0) == k)) {
    stop("p0 must be a single value or a vector of length k (nrow(bin_matrix))")
  }
  if (any(p0 <= 0 | p0 >= 1)) {
    stop("p0 must be strictly between 0 and 1")
  }
  p0_vec <- rep(p0, length.out = k)

  # mu0 = expected value of sum at each time point = sum(p0_i)
  mu0 <- sum(p0_vec)

  # sigma2_0 = variance of sum at each time point = sum(p0_i * (1 - p0_i))
  sigma2_0 <- sum(p0_vec * (1 - p0_vec))

  # Initialize variables (exactly as in original simulation code)
  cum_sum <- 0            # Cumulative sum of C_t across time
  r_prev <- 0             # Previous EWMA value (r_0 = 0)

  # Pre-allocate vectors to store history
  r_history <- numeric(max_time)    # Store EWMA values
  UCL_history <- numeric(max_time)  # Store upper control limits
  LCL_history <- numeric(max_time)  # Store lower control limits

  # Initialize signal tracking variables
  signal_time <- NA
  signal_detected <- FALSE
  signal_times <- integer(0)

  # ========================================================================
  # STEP 3: Main Monitoring Loop
  # ========================================================================

  for (t in 1:max_time) {

    # Get the binary vector for this time point (column t)
    bin_t <- bin_matrix[, t]

    # Calculate C_t = sum of binary indicators across streams at time t
    C_t <- sum(bin_t)

    # Update cumulative sum
    cum_sum <- cum_sum + C_t

    # Compute standardized statistic W_t
    W_t <- (cum_sum - mu0 * t) / sqrt(t * sigma2_0)

    # Update EWMA statistic r_t
    r_t <- lambda * W_t + (1 - lambda) * r_prev

    # Store r_t in history
    r_history[t] <- r_t

    # Get exact variance at time t from precomputed cache
    v_t <- var_cache[t]

    # Compute control limits
    UCL_t <- L * sqrt(v_t)
    LCL_t <- -L * sqrt(v_t)

    # Store control limits for plotting
    UCL_history[t] <- UCL_t
    LCL_history[t] <- LCL_t

    # Check for signal: r_t exceeds either control limit
    if (r_t > UCL_t || r_t < LCL_t) {
      if (!signal_detected) {
        signal_time <- t
      }
      signal_detected <- TRUE
      signal_times <- c(signal_times, t)
      if (stop_at_signal) {
        break
      }
    }

    # Update previous EWMA value for next iteration
    r_prev <- r_t
  }

  # ========================================================================
  # STEP 4: Prepare Output with Proper Class Attribute
  # ========================================================================

  # Determine signal time (if no signal, use max_time)
  if (stop_at_signal) {
    if (signal_detected) {
      T_sig <- signal_time
    } else {
      T_sig <- max_time
      warning("No signal detected within max_time")
    }
  } else {
    T_sig <- max_time
    if (!signal_detected) {
      warning("No signal detected within max_time")
    }
  }

  # Truncate binary matrix to signal time
  bin_matrix_truncated <- bin_matrix[, 1:T_sig, drop = FALSE]

  # Calculate successes per stream
  successes <- rowSums(bin_matrix_truncated)

  # Create the result list
  result_list <- list(
    signal_time = signal_time,
    signal_detected = signal_detected,
    signal_times = signal_times,
    stop_at_signal = stop_at_signal,
    T_sig = T_sig,
    r_history = r_history[1:T_sig],
    UCL_history = UCL_history[1:T_sig],
    LCL_history = LCL_history[1:T_sig],
    bin_matrix = bin_matrix_truncated,
    successes = successes,
    k = k,
    lambda = lambda,
    L = L,
    p0 = p0
  )

  # Set the class attribute for S3 method dispatch
  class(result_list) <- "csb_ewma"

  return(result_list)
}
#' CSB-EWMA Control Chart
#'
#' Runs the Cumulative Standardized Binomial EWMA control chart on multiple stream data.
#'
#' By default (\code{stop_at_signal = TRUE}), monitoring stops at the first
#' signal and post-hoc identification is performed once, matching the
#' behavior of every released version of this function through 1.1.0.
#' Setting \code{stop_at_signal = FALSE} instead continues monitoring
#' through \code{max_time} and performs post-hoc identification separately
#' at every signal recorded in \code{signal_times}, using only the data
#' observed up to and including that signal. See \code{\link{run_csb_ewma}}
#' for the full explanation of continuation mode; it applies identically
#' here since this function calls that one internally.
#'
#' @param data A matrix of binary indicators (0/1) with streams as rows and time as columns
#' @param lambda Smoothing parameter for EWMA (0 < lambda <= 1)
#' @param L Control limit multiplier
#' @param p0 In-control proportion (default = 0.5). May instead be a numeric
#'   vector of length k (one value per stream), passed through to
#'   \code{\link{run_csb_ewma}} and to the post-hoc identification step;
#'   see that function's documentation for the exact generalization used.
#' @param max_time Maximum time points to monitor (default = NULL uses all)
#' @param stop_at_signal Logical. If TRUE (default), monitoring stops at the
#'   first signal and post-hoc identification is performed once, stored in
#'   \code{flagged} (backward compatible with every released version through
#'   1.1.0). If FALSE, monitoring continues through \code{max_time} and
#'   post-hoc identification is performed separately at each signal in
#'   \code{signal_times}, stored as a named list in \code{flagged_by_signal}
#'   (names are the signal times as character strings); \code{flagged} is
#'   then set to the first signal's result for convenience.
#' @param distribution If data is continuous, specify distribution
#' @param posthoc_method Method for post-hoc identification (default = "BH")
#' @param alpha Significance level for post-hoc (default = 0.05)
#' @param verbose Logical. If TRUE, prints informational messages about whether the input data required dichotomization. Default FALSE (silent).
#' @return A list of class "csb_ewma" containing chart results and flagged
#'   streams. When \code{stop_at_signal = FALSE}, also contains
#'   \code{signal_times} (every time point at which a signal was detected)
#'   and \code{flagged_by_signal} (post-hoc identification results at each
#'   of those times).
#' @seealso [run_csb_ewma()] for the lower-level engine this function calls, useful directly when a variance cache should be precomputed once and reused across many chart runs (for example, in simulation studies).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' bin_data <- matrix(rbinom(10*100, 1, 0.5), nrow = 10, ncol = 100)
#' for(i in 1:3) bin_data[i, ] <- rbinom(100, 1, 0.8)
#' result <- csb_ewma(bin_data, lambda = 0.175, L = 1.375)
#' print(result)
#' plot(result)
#'
#' # Continue monitoring past the first signal and see every signal time
#' result2 <- csb_ewma(bin_data, lambda = 0.175, L = 1.375,
#'                      stop_at_signal = FALSE)
#' print(result2$signal_times)
#' }
csb_ewma <- function(data, lambda, L, p0 = 0.5, max_time = NULL,
                     stop_at_signal = TRUE, distribution = NULL,
                     posthoc_method = "BH", alpha = 0.05, verbose = FALSE) {

  # ========================================================================
  # STEP 1: Input Validation
  # ========================================================================

  # Check that data is a matrix
  if (!is.matrix(data)) {
    stop("data must be a matrix")
  }

  # Check that lambda is between 0 and 1
  if (lambda <= 0 || lambda > 1) {
    stop("lambda must be between 0 and 1")
  }

  # Check that L is positive
  if (L <= 0) {
    stop("L must be positive")
  }

  # ========================================================================
  # STEP 2: Determine if data is binary or continuous
  # ========================================================================

  # Check if data contains only 0s and 1s (binary)
  is_binary <- all(data %in% c(0, 1))

  if (is_binary) {
    # Data is already binary - use as-is
    bin_matrix <- as.matrix(data)
    if (verbose) message("Data is already binary; no dichotomization needed.\n")
  } else {
    # Data is continuous - need to dichotomize
    if (is.null(distribution)) {
      stop("For continuous data, please specify 'distribution' parameter")
    }

    # Convert continuous data to binary using dichotomize_data
    bin_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    for (t in 1:ncol(data)) {
      bin_matrix[, t] <- dichotomize_data(data[, t], distribution, p0 = p0)
    }
    if (verbose) message("Continuous data dichotomized using", distribution, "distribution\n")
  }

  # ========================================================================
  # STEP 3: Precompute variance
  # ========================================================================

  # Determine max_time
  if (is.null(max_time)) {
    max_time <- ncol(bin_matrix)
  } else {
    max_time <- min(max_time, ncol(bin_matrix))
  }

  # Precompute variance vector
  var_cache <- precompute_variance(lambda, max_t = max_time, converge_t = 500)

  # ========================================================================
  # STEP 4: Run CSB-EWMA chart
  # ========================================================================

  chart_result <- run_csb_ewma(bin_matrix, lambda, L, var_cache,
                               max_time = max_time, p0 = p0,
                               stop_at_signal = stop_at_signal)

  # ========================================================================
  # STEP 5: Perform post-hoc identification if signal detected
  # ========================================================================

  if (chart_result$signal_detected) {
    if (stop_at_signal) {
      message("Signal detected at time t =", chart_result$signal_time, "\n")

      # Identify out-of-control streams using specified method
      flagged <- identify_ooc(chart_result$bin_matrix,
                               alpha = alpha,
                               method = posthoc_method,
                               p0 = p0)

      # Add flagged results to the output
      chart_result$flagged <- flagged
      chart_result$flagged_by_signal <- stats::setNames(
        list(flagged), as.character(chart_result$signal_time)
      )
      chart_result$posthoc_method <- posthoc_method
      chart_result$alpha <- alpha

      # Print summary of flagged streams
      flagged_summary <- flagged_streams_summary(flagged)
      if (flagged_summary$total_flagged > 0) {
        message("Flagged streams:", flagged_summary$flagged_streams, "\n")
      } else {
        message("No streams were flagged as out-of-control\n")
      }
    } else {
      message("Signal(s) detected at time t =",
              paste(chart_result$signal_times, collapse = ", "), "\n")

      # Post-hoc identification is performed separately at each signal,
      # using only the data observed up to and including that signal time.
      # This is a practical extension of the single-signal identification
      # above to continuation mode; it is not itself a procedure described
      # in the dissertation.
      flagged_by_signal <- lapply(chart_result$signal_times, function(t) {
        identify_ooc(chart_result$bin_matrix[, 1:t, drop = FALSE],
                     alpha = alpha,
                     method = posthoc_method,
                     p0 = p0)
      })
      names(flagged_by_signal) <- as.character(chart_result$signal_times)

      chart_result$flagged_by_signal <- flagged_by_signal
      chart_result$flagged <- flagged_by_signal[[1]]
      chart_result$posthoc_method <- posthoc_method
      chart_result$alpha <- alpha

      for (t in chart_result$signal_times) {
        flagged_summary <- flagged_streams_summary(flagged_by_signal[[as.character(t)]])
        if (flagged_summary$total_flagged > 0) {
          message("Signal at t =", t, ": flagged streams:",
                  flagged_summary$flagged_streams, "\n")
        } else {
          message("Signal at t =", t, ": no streams were flagged as out-of-control\n")
        }
      }
    }
  } else {
    message("No signal detected within", max_time, "time points\n")
    chart_result$flagged <- NULL
    chart_result$flagged_by_signal <- NULL
  }

  # ========================================================================
  # STEP 6: Return results
  # ========================================================================

  return(chart_result)
}


#' Prints a formatted summary of CSB-EWMA chart results.
#'
#' @param x csb_ewma object from csb_ewma() or run_csb_ewma() function
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the object
#' @export
# ----------------------------------------------------------------------------
print.csb_ewma <- function(x, ...) {

  # Print a nice formatted header
  message("\n========================================\n")
  message("CSB-EWMA Control Chart Results\n")
  message("========================================\n")

  # Print chart parameters
  message("\nChart Parameters:\n")
  message("  lambda =", x$lambda, "\n")
  message("  L =", x$L, "\n")
  message("  p0 =", x$p0, "\n")
  message("  Number of streams (k) =", x$k, "\n")

  # Print signal detection results
  message("\nSignal Detection:\n")
  if (x$signal_detected) {
    message("  Signal detected at time t =", x$signal_time, "\n")
    message("  Signal time T_sig =", x$T_sig, "\n")
  } else {
    message("  No signal detected\n")
    message("  Maximum time monitored =", x$T_sig, "\n")
  }

  # Print post-hoc identification results if available
  if (!is.null(x$flagged)) {
    cat("\nPost-hoc Identification (", x$posthoc_method, ", alpha = ", x$alpha, "):\n", sep="")
    flagged_count <- sum(x$flagged$flagged)
    message("  Flagged streams:", flagged_count, "out of", x$k, "\n")

    if (flagged_count > 0) {
      flagged_streams <- x$flagged$stream[x$flagged$flagged]
      message("  Stream numbers:", paste(flagged_streams, collapse = ", "), "\n")
    }
  }

  message("\n========================================\n")
  invisible(x)
}
