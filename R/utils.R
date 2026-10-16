#' Generate Laplace (Double Exponential) Random Variables
#'
#' Generates random numbers from the Laplace distribution using inverse transform sampling.
#'
#' @param n Number of observations to generate
#' @param location Location parameter (median) of the distribution (default = 0)
#' @param scale Scale parameter (spread) of the distribution (default = 1)
#' @return A numeric vector of length n containing Laplace random variables
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rlaplace(100, location = 0, scale = 1)
#' }
rlaplace <- function(n, location = 0, scale = 1) {
  u <- stats::runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}
#' Generate Continuous Data from Specified Distribution
#'
#' Generates continuous observations from normal, Laplace, uniform, or exponential distributions
#' with optional shift parameter for out-of-control simulation.
#'
#' @param distribution Character string: "normal", "laplace", "uniform", or "exponential"
#' @param n Number of observations to generate
#' @param shift Amount to shift the proportion (default = 0)
#' @param p0 In-control proportion (default = 0.5)
#' @return A numeric vector of length n containing the generated data
#' @export
#'
#' @examples
#' data <- generate_continuous_data("normal", n = 100, shift = 0.2)
# ----------------------------------------------------------------------------
generate_continuous_data <- function(distribution, n, shift = 0, p0 = 0.5) {

  # ========================================================================
  # CASE 1: Normal Distribution
  # ========================================================================
  # For normal, we shift the mean. The relationship is:
  #   p1 = P(X > threshold) = 1 - Φ((threshold - μ)/σ)
  #   Therefore μ = threshold - σ * Φ^{-1}(1 - p1)
  # ========================================================================
  if (distribution == "normal") {
    # Compute the in-control threshold (p0-th quantile)
    # For p0 = 0.5, this is the median = 0
    threshold_in <- stats::qnorm(p0)

    # Calculate out-of-control proportion, bounded to avoid numerical issues
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)  # Keep between 0.01 and 0.99

    # Compute the mean shift needed to achieve p1 proportion above threshold
    mean_shift <- threshold_in - stats::qnorm(1 - p1)

    # Generate normal data with shifted mean, unit variance
    data <- stats::rnorm(n, mean = mean_shift, sd = 1)
  }

  else if (distribution == "laplace") {
    # Compute in-control threshold (p0-th quantile)
    if (p0 < 0.5) {
      threshold_in <- log(2 * p0)  # For p0 < 0.5
    } else {
      threshold_in <- -log(2 * (1 - p0))  # For p0 > 0.5
    }

    # Calculate out-of-control proportion
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)

    # Compute location shift to achieve p1 proportion above threshold
    if (p1 < 0.5) {
      location <- threshold_in - log(2 * p1)
    } else {
      location <- threshold_in + log(2 * (1 - p1))
    }

    # Generate Laplace data with shifted location
    data <- rlaplace(n, location = location, scale = 1)
  }


  else if (distribution == "uniform") {
    # In-control threshold is the p0-th quantile
    threshold_in <- p0

    # Calculate out-of-control proportion
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)

    # Compute the shift amount for the uniform distribution
    shift_amount <- threshold_in + p1 - 1

    # Generate uniform data with shifted bounds
    data <- stats::runif(n, min = shift_amount, max = 1 + shift_amount)
  }


  else if (distribution == "exponential") {
    # In-control threshold is the p0-th quantile
    threshold_in <- stats::qexp(p0, rate = 1)

    # Calculate out-of-control proportion
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)

    # Compute rate parameter to achieve p1 proportion above threshold
    rate <- -log(1 - p1) / threshold_in

    # Generate exponential data with adjusted rate
    data <- stats::rexp(n, rate = rate)
  }


  else {
    stop("Unknown distribution. Choose from: 'normal', 'laplace', 'uniform', 'exponential'")
  }

  return(data)
}


#' Dichotomize Continuous Data to Binary Indicators
#'
#' Converts continuous observations to binary indicators based on whether each
#' observation exceeds the in-control median or specified quantile.
#'
#' @param data Numeric vector of continuous observations
#' @param distribution Character string specifying the distribution type
#' @param p0 In-control proportion (default = 0.5). May instead be a numeric
#'   vector the same length as data, giving a separate in-control
#'   proportion (and therefore threshold) for each observation; this is
#'   used, for example, to dichotomize one time slice across k streams that
#'   do not all share the same in-control proportion.
#' @return Integer vector of binary indicators (0 or 1)
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' binary <- dichotomize_data(x, distribution = "normal", p0 = 0.5)
#'
#' # Per-observation (e.g. per-stream) in-control proportions
#' x2 <- rnorm(5)
#' binary2 <- dichotomize_data(x2, distribution = "normal",
#'                              p0 = c(0.5, 0.5, 0.4, 0.6, 0.5))
# ----------------------------------------------------------------------------
dichotomize_data <- function(data, distribution, p0 = 0.5) {

  # ========================================================================
  # STEP 1: Compute the threshold based on distribution type
  # ========================================================================

  # p0 may be a single in-control proportion applied to every observation,
  # or a numeric vector the same length as data giving a separate in-control
  # proportion per observation (for example, one value per stream when data
  # is a single time slice across streams). All branches below are written
  # with vectorized R functions (qnorm/qexp/ifelse) so a vector p0 produces
  # a vector of per-observation thresholds th, and a scalar p0 reproduces
  # the original single-threshold behavior exactly.
  if (!(length(p0) == 1 || length(p0) == length(data))) {
    stop("p0 must be a single value or a vector the same length as data")
  }

  if (distribution == "normal") {
    # For normal distribution, use qnorm() to get the p0-th quantile
    # When p0 = 0.5, this returns 0 (the median)
    th <- stats::qnorm(p0)
  }
  else if (distribution == "laplace") {
    # For Laplace distribution with location = 0, scale = 1
    # The p0-th quantile is given by (vectorized over p0):
    th <- ifelse(p0 < 0.5, log(2 * p0), -log(2 * (1 - p0)))
  }
  else if (distribution == "uniform") {
    # For Uniform[0,1], the p0-th quantile is simply p0
    th <- p0
  }
  else if (distribution == "exponential") {
    # For Exponential(rate = 1), use qexp()
    th <- stats::qexp(p0, rate = 1)
  }
  else {
    # If distribution is not recognized, stop with error
    stop("Unknown distribution. Choose from: 'normal', 'laplace', 'uniform', 'exponential'")
  }

  # ========================================================================
  # STEP 2: Convert to binary indicators
  # ========================================================================
  # as.integer() converts TRUE to 1 and FALSE to 0
  # This creates a binary indicator where 1 means "exceeded the threshold"
  as.integer(data > th)
}
