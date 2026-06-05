#' Compute exact CSB-EWMA variance for a single time point
#'
#' This function implements the exact variance formula derived in Theorem 2.
#' The variance is computed using double summation over the covariance structure
#' of the standardized statistics. For a given smoothing parameter lambda and
#' time point t, this returns Var(r_t) as defined in Equation (19) of the
#' supplementary material.
#'
#' @param lambda Smoothing parameter for EWMA. Must be between 0 and 1.
#'               Typical values range from 0.05 to 0.5.
#' @param t Time point (positive integer). The variance is computed for
#'          this specific time point.
#'
#' @return The exact variance Var(r_t) at time t.
#'         For t = 0, returns 0.
#'         
#' @export
#' 
#' @examples
#' # Compute variance at time 10 for lambda = 0.175
#' var_rt_exact_single(lambda = 0.175, t = 10)
#'
#' # Compute variance at time 50 for lambda = 0.15
#' var_rt_exact_single(lambda = 0.15, t = 50)
#'
var_rt_exact_single <- function(lambda, t) {
  
  # Input validation: if t is 0, variance is 0
  # This handles the initial condition before any observations
  if (t == 0) return(0)
  
  # Initialize sum accumulator S
  # This will accumulate the double summation terms from Equation (19)
  S <- 0 
  
  # First term: double summation over j = 1 to t-1 and i = j+1 to t
  # This accounts for the covariance between different time points (i != j)
  for (j in 1:(t-1)) {
    
    # Compute exponential weights for all i > j
    # The exponent is (2t - i - j) as in Equation (19)
    weights_j <- (1 - lambda)^(2*t - j - (j+1):t)
    
    # Compute covariance terms sqrt(j / i) from Equation (15)
    # This comes from Cov(W_i, W_j) = sqrt(min(i,j))/sqrt(max(i,j))
    cov_terms <- sqrt(j / ((j+1):t))
    
    # Add contribution to sum (multiply by 2 for symmetry i < j and i > j)
    # The factor 2 accounts for both triangles of the double summation
    S <- S + 2 * sum(weights_j * cov_terms)
  }
  
  # Second term: diagonal terms where i = j
  # These are the variance terms for each individual time point
  diag_weights <- (1 - lambda)^(2*t - 2*(1:t))
  S <- S + sum(diag_weights)
  
  # Multiply by lambda^2 as in Equation (19)
  # The smoothing parameter squared scales the entire variance
  return(lambda^2 * S)
}


#' Precompute variance vector for all time points
#'
#' This function precomputes the exact variance for all time points from 1 to
#' max_t. For efficiency, it computes exact variances up to a convergence
#' threshold (converge_t) and then sets remaining values to 1 (asymptotic).
#' Based on the derivation, variance reaches 99% of asymptotic value by t=227,
#' so converge_t = 500 is safe and efficient.
#'
#' @param lambda Smoothing parameter for EWMA
#' @param max_t Maximum time point to precompute variance for
#' @param converge_t Time point after which variance is set to 1 (asymptotic)
#'
#' @return A numeric vector of length max_t containing variance at each time point
#'
#' @export
#'
#' @examples
#' # Precompute variance for lambda = 0.175, up to t = 1000
#' var_cache <- precompute_variance(lambda = 0.175, max_t = 1000, converge_t = 500)
#'
precompute_variance <- function(lambda, max_t, converge_t = 500) {
  
  # Initialize a numeric vector of zeros of length max_t
  varvec <- numeric(max_t)
  
  # Compute exact variance for time points up to converge_t (or max_t if smaller)
  # This loop calls var_rt_exact_single for each time point individually
  for (t in 1:min(converge_t, max_t)) {
    varvec[t] <- var_rt_exact_single(lambda, t)
  }
  
  # For time points beyond converge_t, use asymptotic variance of 1
  # This is valid because variance converges to 1 as t -> infinity
  if (max_t > converge_t) {
    varvec[(converge_t+1):max_t] <- 1
  }
  
  # Return the precomputed variance vector
  return(varvec)
}