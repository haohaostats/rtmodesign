estimate_rtmo_moments <- function(residuals) {
  residuals <- residuals[is.finite(residuals)]
  if (length(residuals) < 10L) {
    stop("At least 10 finite new-treatment residuals are required for moment estimation.",
         call. = FALSE)
  }
  centered <- residuals - mean(residuals)
  values <- vapply(2:6, function(order) mean(centered^order), numeric(1L))
  names(values) <- paste0("mu", 2:6)
  if (!is.finite(values[["mu2"]]) || values[["mu2"]] <= 0) {
    stop("Residual variance must be positive.", call. = FALSE)
  }
  sigma3 <- matrix(
    c(
      values[["mu2"]], values[["mu3"]], values[["mu4"]],
      values[["mu3"]], values[["mu4"]] - values[["mu2"]]^2,
      values[["mu5"]] - values[["mu2"]] * values[["mu3"]],
      values[["mu4"]], values[["mu5"]] - values[["mu2"]] * values[["mu3"]],
      values[["mu6"]] - values[["mu3"]]^2
    ),
    nrow = 3L, byrow = TRUE
  )
  eigenvalues <- eigen(sigma3, symmetric = TRUE, only.values = TRUE)$values
  list(
    central = values,
    standard_deviation = sqrt(values[["mu2"]]),
    skewness = values[["mu3"]] / values[["mu2"]]^(3 / 2),
    kurtosis = values[["mu4"]] / values[["mu2"]]^2,
    third_order_covariance = sigma3,
    eigenvalues = eigenvalues,
    positive_definite = all(eigenvalues > 0),
    condition_number = if (all(eigenvalues > 0)) max(eigenvalues) / min(eigenvalues) else Inf,
    sample_size = length(residuals)
  )
}

