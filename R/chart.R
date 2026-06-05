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
#'     \item If r_t > UCL_t or r_t < LCL_t, signal at time t and break
#'     \item Otherwise, update r_prev = r_t and continue
#'   }
#' }
#'
#' @param bin_matrix A matrix of binary indicators (streams as rows, time as columns)
#' @param lambda Smoothing parameter for EWMA (0 < lambda <= 1)
#' @param L Control limit multiplier
#' @param var_cache Precomputed variance vector from precompute_variance()
#' @param max_time Maximum time points to monitor (default = NULL uses all)
#' @param p0 In-control proportion(s). Either a single number (same for all streams) or a numeric vector of length equal to number of rows in bin_matrix.
#'           Example p0 = c(0.5, 0.6, 0.7)   # different p0 for stream1, stream2, stream3.
#' @return A list of class "csb_ewma" containing chart results
#' @export
#'
#' @examples
#' bin_matrix <- matrix(rbinom(10*200, 1, 0.5), nrow = 10, ncol = 200)
#' for(i in 1:3) bin_matrix[i, ] <- rbinom(100, 1, 0.8)
#' var_cache <- precompute_variance(0.175, max_t = 200)
#' result <- run_csb_ewma(bin_matrix, lambda = 0.175, L = 1.375, var_cache)
#' message(paste("Signal at time:", result$signal_time))
# ----------------------------------------------------------------------------
run_csb_ewma <- function(bin_matrix, lambda, L, var_cache,
                         max_time = NULL, p0 = 0.5) {

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

  # ========================================================================
  # STEP 2b: Handle p0 (single value or per‑stream vector)
  # ========================================================================
  if (length(p0) == 1) {
    p0 <- rep(p0, k)
  } else if (length(p0) != k) {
    stop("p0 must be length 1 or length equal to number of streams (k = ", k, ")")
  }
  mu0 <- sum(p0)
  sigma2_0 <- sum(p0 * (1 - p0))

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
      signal_time <- t
      signal_detected <- TRUE
      break
    }

    # Update previous EWMA value for next iteration
    r_prev <- r_t
  }

  # ========================================================================
  # STEP 4: Prepare Output with Proper Class Attribute
  # ========================================================================

  # Determine signal time (if no signal, use max_time)
  if (signal_detected) {
    T_sig <- signal_time
  } else {
    T_sig <- max_time
    warning("No signal detected within max_time")
  }

  # Truncate binary matrix to signal time
  bin_matrix_truncated <- bin_matrix[, 1:T_sig, drop = FALSE]

  # Calculate successes per stream
  successes <- rowSums(bin_matrix_truncated)

  # Create the result list
  result_list <- list(
    signal_time = signal_time,
    signal_detected = signal_detected,
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
#' @param data A matrix of binary indicators (0/1) with streams as rows and time as columns
#' @param lambda Smoothing parameter for EWMA (0 < lambda <= 1)
#' @param L Control limit multiplier
#' @param p0 In-control proportion(s). Either a single number (same for all streams)
#'          or a numeric vector of length equal to number of rows in data.
#' @param max_time Maximum time points to monitor (default = NULL uses all)
#' @param posthoc_method Method for post-hoc identification (default = "BH")
#' @param alpha Significance level for post-hoc (default = 0.05)
#' @param verbose If TRUE, print informational messages (default = TRUE)
#' @return A list of class "csb_ewma" containing chart results and flagged streams
#' @export
#'
#' @examples{
#' set.seed(123)
#' bin_data <- matrix(rbinom(10*100, 1, 0.5), nrow = 10, ncol = 100)
#' for(i in 1:3) bin_data[i, ] <- rbinom(100, 1, 0.8)
#' result <- csb_ewma(bin_data, lambda = 0.175, L = 1.375)
#' message(result)
#' plot(result)
#' }
# ----------------------------------------------------------------------------
csb_ewma <- function(data, lambda, L, p0 = 0.5, max_time = NULL,
                     posthoc_method = "BH", alpha = 0.05, verbose = TRUE) {

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
  # STEP 2: Check data is binary (no continuous conversion)
  # ========================================================================

  if (!all(data %in% c(0, 1))) {
    stop("data must contain only 0s and 1s (binary indicators)")
  }

  bin_matrix <- as.matrix(data)
  if (verbose) message("Data is binary. Using as-is.")

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
                               max_time = max_time, p0 = p0)

  # ========================================================================
  # STEP 5: Perform post-hoc identification if signal detected
  # ========================================================================

  if (chart_result$signal_detected) {
    if (verbose) message("Signal detected at time t = ", chart_result$signal_time)

    # Identify out-of-control streams using specified method
    flagged <- identify_ooc(chart_result$bin_matrix,
                            alpha = alpha,
                            method = posthoc_method,
                            p0 = chart_result$p0,
                            verbose = verbose)

    # Add flagged results to the output
    chart_result$flagged <- flagged
    chart_result$posthoc_method <- posthoc_method
    chart_result$alpha <- alpha

    # message summary of flagged streams
    flagged_summary <- flagged_streams_summary(flagged, verbose = verbose)
    if (flagged_summary$total_flagged > 0) {
      if (verbose) message("Flagged streams: ", paste(flagged_summary$flagged_streams, collapse = ", "))
    } else {
      if (verbose) message("No streams were flagged as out-of-control")
    }
  } else {
    if (verbose) message("No signal detected within ", max_time, " time points")
    chart_result$flagged <- NULL
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
    message("\nPost-hoc Identification (", x$posthoc_method, ", alpha = ", x$alpha, "):\n", sep="")
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
