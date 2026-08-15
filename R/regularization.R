matrix_condition_number <- function(matrix) {
  values <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(values)) || min(values) <= 0) return(Inf)
  max(values) / min(values)
}

regularize_rtmo_moments <- function(moments, threshold = 100) {
  covariance <- moments$third_order_covariance
  diagonal <- diag(covariance)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("The empirical high-order moment covariance has a non-positive diagonal.",
         call. = FALSE)
  }
  scale <- sqrt(diagonal)
  correlation <- covariance / outer(scale, scale)
  correlation <- (correlation + t(correlation)) / 2
  raw_condition <- matrix_condition_number(correlation)

  condition_at <- function(lambda) {
    candidate <- (1 - lambda) * correlation + lambda * diag(nrow(correlation))
    matrix_condition_number(candidate)
  }
  if (is.finite(raw_condition) && raw_condition <= threshold) {
    lambda <- 0
  } else {
    lower <- 0
    upper <- 1
    for (iteration in seq_len(80L)) {
      midpoint <- (lower + upper) / 2
      if (condition_at(midpoint) <= threshold) upper <- midpoint else lower <- midpoint
    }
    lambda <- upper
  }

  regularized_correlation <-
    (1 - lambda) * correlation + lambda * diag(nrow(correlation))
  regularized_covariance <- regularized_correlation * outer(scale, scale)
  list(
    applied = lambda > sqrt(.Machine$double.eps),
    lambda = lambda,
    threshold = threshold,
    raw_correlation = correlation,
    regularized_correlation = regularized_correlation,
    raw_covariance = covariance,
    regularized_covariance = regularized_covariance,
    raw_condition_number = raw_condition,
    regularized_condition_number = matrix_condition_number(regularized_correlation)
  )
}

