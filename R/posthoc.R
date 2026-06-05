#' Calculate p-values for all streams
#'
#' Computes exact binomial p-values for each stream based on successes.
#'
#' @param bin_matrix Matrix of binary indicators (streams as rows, time as columns)
#' @param p0 In-control proportion(s). Either a single number or a vector of length = nrow(bin_matrix)
#' @return Numeric vector of p-values
#' @export
#'
#' @examples
#' bin_mat <- matrix(rbinom(500, 1, 0.5), nrow = 10, ncol = 50)
#' pvals <- calculate_pvalues(bin_mat, p0 = 0.5)
calculate_pvalues <- function(bin_matrix, p0 = 0.5) {
  k <- nrow(bin_matrix)
  T_sig <- ncol(bin_matrix)
  successes <- rowSums(bin_matrix)
  # Handle p0: single value or vector
  if (length(p0) == 1) p0 <- rep(p0, k)
  if (length(p0) != k) stop("p0 must be length 1 or length equal to number of streams")
  pvals <- sapply(1:k, function(i) {
    1 - stats::pbinom(successes[i] - 1, size = T_sig, prob = p0[i])
  })
  return(pvals)
}

#' Apply multiple testing corrections
#'
#' Applies Bonferroni, Holm, or Benjamini-Hochberg corrections.
#'
#' @param pvals Numeric vector of raw p-values
#' @param method Correction method: "BH", "bonferroni", "holm", or "raw"
#' @param alpha Significance level (default = 0.05)
#' @return List with adjusted_pvals and flags
#' @export
#'
#' @examples
#' pvals <- c(0.001, 0.01, 0.03, 0.10, 0.50)
#' result <- apply_multiple_testing(pvals, method = "BH", alpha = 0.05)
apply_multiple_testing <- function(pvals, method = "BH", alpha = 0.05) {
  if (method == "raw") {
    flags <- pvals < alpha
    return(list(adjusted_pvals = pvals, flags = flags))
  }
  if (method == "bonferroni") {
    adjusted_pvals <- stats::p.adjust(pvals, method = "bonferroni")
    flags <- adjusted_pvals < alpha
    return(list(adjusted_pvals = adjusted_pvals, flags = flags))
  }
  if (method == "holm") {
    adjusted_pvals <- stats::p.adjust(pvals, method = "holm")
    flags <- adjusted_pvals < alpha
    return(list(adjusted_pvals = adjusted_pvals, flags = flags))
  }
  if (method == "BH") {
    adjusted_pvals <- stats::p.adjust(pvals, method = "BH")
    flags <- adjusted_pvals < alpha
    return(list(adjusted_pvals = adjusted_pvals, flags = flags))
  }
  stop("Unknown method. Choose from: 'BH', 'bonferroni', 'holm', 'raw'")
}


#' Identify out-of-control streams
#' Main function for post-hoc identification using multiple testing corrections.
#'
#' @param bin_matrix Matrix of binary indicators (streams as rows, time as columns)
#' @param alpha Significance level (default = 0.05)
#' @param method Correction method (default = "BH")
#' @param p0 In-control proportion(s). Either a single number or a vector of length = nrow(bin_matrix)
#' @param verbose If TRUE, print informational messages (default = TRUE)
#' @return Data frame with p-values and flags for each stream
#' @export
#'
#' @examples
#' set.seed(123)
#' bin_mat <- matrix(rbinom(10*100, 1, 0.5), nrow = 10, ncol = 100)
#' for(i in 1:3) bin_mat[i, ] <- rbinom(100, 1, 0.8)
#' result <- identify_ooc(bin_mat, p0 = 0.5)
#' print(result[result$flagged, ])
identify_ooc <- function(bin_matrix, alpha = 0.05, method = "BH", p0 = 0.5, verbose = TRUE) {
  if (!is.matrix(bin_matrix)) stop("bin_matrix must be a matrix")
  if (!method %in% c("BH", "bonferroni", "holm", "raw")) {
    stop("method must be one of: 'BH', 'bonferroni', 'holm', 'raw'")
  }
  k <- nrow(bin_matrix)
  T_sig <- ncol(bin_matrix)
  successes <- rowSums(bin_matrix)
  proportion <- successes / T_sig
  pvals <- calculate_pvalues(bin_matrix, p0 = p0)
  corrected <- apply_multiple_testing(pvals, method = method, alpha = alpha)
  results <- data.frame(
    stream = 1:k,
    successes = successes,
    proportion = proportion,
    p_value = pvals,
    adjusted_p_value = corrected$adjusted_pvals,
    flagged = corrected$flags,
    stringsAsFactors = FALSE
  )
  results <- results[order(-results$flagged, results$p_value), ]
  rownames(results) <- NULL
  return(results)
}


#' Summarize flagged streams
#'
#' Creates a summary of which streams were flagged.
#'
#' @param flagged_results Output from identify_ooc()
#' @param verbose If TRUE, print informational messages (default = TRUE)
#' @return List with flagged streams and summary table
#' @export
#'
#' @examples
#' bin_mat <- matrix(rbinom(10*100, 1, 0.5), nrow = 10, ncol = 100)
#' for(i in 1:3) bin_mat[i, ] <- rbinom(100, 1, 0.8)
#' result <- identify_ooc(bin_mat)
#' summary <- flagged_streams_summary(result)
#' print(summary$flagged_streams)
flagged_streams_summary <- function(flagged_results, verbose = TRUE) {
  flagged_streams <- flagged_results$stream[flagged_results$flagged]
  total_flagged <- length(flagged_streams)
  if (total_flagged > 0) {
    summary_table <- flagged_results[flagged_results$flagged,
                                     c("stream", "successes", "proportion",
                                       "p_value", "adjusted_p_value")]
  } else {
    summary_table <- NULL
    if (verbose) message("No streams were flagged as out-of-control")
  }
  return(list(
    flagged_streams = flagged_streams,
    total_flagged = total_flagged,
    summary_table = summary_table
  ))
}
