`%||%` <- function(x, y) if (is.null(x)) y else x

predict_rtmo_mean <- function(model, dose, theta) {
  switch(
    model,
    linear = theta[[1L]] * dose,
    mm = theta[[1L]] * dose / (theta[[2L]] + dose),
    emax = theta[[1L]] * dose / (theta[[2L]] + dose),
    sigmoid_emax = theta[[1L]] * dose^theta[[3L]] /
      (theta[[2L]]^theta[[3L]] + dose^theta[[3L]]),
    exponential = theta[[1L]] * (exp(dose / theta[[2L]]) - 1),
    stop("Unsupported model.", call. = FALSE)
  )
}

fit_rtmo_model <- function(model, dose, outcome) {
  if (length(dose) != length(outcome) || !length(dose)) {
    stop("Internal model-fitting inputs are inconsistent.", call. = FALSE)
  }
  if (model == "linear") {
    denominator <- sum(dose^2)
    if (denominator <= 0) stop("The linear model is not identifiable.", call. = FALSE)
    theta <- sum(dose * outcome) / denominator
    if (!is.finite(theta) || theta <= 0) {
      stop("The fitted dose-response trend is not increasing. Check outcome orientation.",
           call. = FALSE)
    }
    fitted <- predict_rtmo_mean(model, dose, theta)
    return(list(
      name = model, coefficients = c(slope = theta), fitted_values = fitted,
      residuals = outcome - fitted, sse = sum((outcome - fitted)^2),
      convergence = 0L, message = "analytic solution"
    ))
  }

  positive_doses <- unique(dose[dose > 0])
  required_positive <- if (model == "sigmoid_emax") 3L else 2L
  if (length(positive_doses) < required_positive) {
    stop(sprintf("Model '%s' requires at least %d distinct positive doses.",
                 model, required_positive), call. = FALSE)
  }
  amplitude <- max(tapply(outcome, dose, mean), na.rm = TRUE)
  if (!is.finite(amplitude) || amplitude <= 0) {
    stop("The fitted dose-response trend is not increasing. Check outcome orientation.",
         call. = FALSE)
  }

  unpack <- switch(
    model,
    mm =, emax = function(z) c(exp(z[[1L]]), exp(z[[2L]])),
    sigmoid_emax = function(z) c(exp(z[[1L]]), exp(z[[2L]]), exp(z[[3L]])),
    exponential = function(z) c(exp(z[[1L]]), exp(z[[2L]]))
  )
  objective <- function(z) {
    theta <- unpack(z)
    prediction <- suppressWarnings(predict_rtmo_mean(model, dose, theta))
    if (any(!is.finite(prediction))) return(.Machine$double.xmax / 100)
    sum((outcome - prediction)^2)
  }

  if (model %in% c("mm", "emax")) {
    starts <- expand.grid(
      amplitude = log(amplitude * c(1.0, 1.5, 3.0)),
      scale = log(c(0.08, 0.25, 0.75))
    )
  } else if (model == "sigmoid_emax") {
    starts <- expand.grid(
      amplitude = log(amplitude * c(1.0, 1.5, 3.0)),
      scale = log(c(0.08, 0.25, 0.75)),
      hill = log(c(0.7, 1.5, 3.0))
    )
  } else {
    starts <- expand.grid(
      amplitude = log(amplitude * c(0.25, 0.75, 1.5)),
      scale = log(c(0.15, 0.5, 1.5))
    )
  }

  fits <- lapply(seq_len(nrow(starts)), function(i) {
    stats::optim(
      par = as.numeric(starts[i, ]), fn = objective, method = "L-BFGS-B",
      lower = rep(log(1e-5), ncol(starts)),
      upper = rep(log(1e3), ncol(starts)),
      control = list(maxit = 5000, factr = 1e7)
    )
  })
  values <- vapply(fits, function(x) x$value, numeric(1L))
  best <- fits[[which.min(values)]]
  if (!is.finite(best$value)) stop("Dose-response model fitting failed.", call. = FALSE)
  theta <- unpack(best$par)
  fitted <- predict_rtmo_mean(model, dose, theta)
  coefficient_names <- switch(
    model,
    mm =, emax = c("maximum", "half_maximum_dose"),
    sigmoid_emax = c("maximum", "half_maximum_dose", "hill"),
    exponential = c("scale", "rate_dose")
  )
  names(theta) <- coefficient_names
  list(
    name = model, coefficients = theta, fitted_values = fitted,
    residuals = outcome - fitted, sse = best$value,
    convergence = best$convergence, message = best$message %||% ""
  )
}

inverse_rtmo_mean <- function(model, mu, theta) {
  if (!is.finite(mu) || mu <= 0) return(NA_real_)
  switch(
    model,
    linear = mu / theta[[1L]],
    mm =, emax = {
      if (mu >= theta[[1L]]) NA_real_ else mu * theta[[2L]] / (theta[[1L]] - mu)
    },
    sigmoid_emax = {
      if (mu >= theta[[1L]]) NA_real_ else
        theta[[2L]] * (mu / (theta[[1L]] - mu))^(1 / theta[[3L]])
    },
    exponential = theta[[2L]] * log1p(mu / theta[[1L]]),
    NA_real_
  )
}

