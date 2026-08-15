gradient_rtmo_mean <- function(model, dose, theta) {
  if (model == "linear") return(c(dose))
  if (model %in% c("mm", "emax")) {
    return(c(
      dose / (theta[[2L]] + dose),
      -theta[[1L]] * dose / (theta[[2L]] + dose)^2
    ))
  }
  if (model == "sigmoid_emax") {
    if (dose <= 0) return(c(0, 0, 0))
    maximum <- theta[[1L]]
    half <- theta[[2L]]
    hill <- theta[[3L]]
    fraction <- dose^hill / (half^hill + dose^hill)
    return(c(
      fraction,
      -maximum * hill * half^(hill - 1) * dose^hill /
        (half^hill + dose^hill)^2,
      maximum * fraction * (1 - fraction) * log(dose / half)
    ))
  }
  if (model == "exponential") {
    exponential <- exp(dose / theta[[2L]])
    return(c(
      exponential - 1,
      -theta[[1L]] * exponential * dose / theta[[2L]]^2
    ))
  }
  stop("Unsupported model.", call. = FALSE)
}

dose_derivative_rtmo_mean <- function(model, dose, theta) {
  switch(
    model,
    linear = theta[[1L]],
    mm =, emax = theta[[1L]] * theta[[2L]] / (theta[[2L]] + dose)^2,
    sigmoid_emax = {
      maximum <- theta[[1L]]
      half <- theta[[2L]]
      hill <- theta[[3L]]
      maximum * hill * half^hill * dose^(hill - 1) /
        (half^hill + dose^hill)^2
    },
    exponential = theta[[1L]] * exp(dose / theta[[2L]]) / theta[[2L]],
    stop("Unsupported model.", call. = FALSE)
  )
}

elemental_rtmo_information <- function(model, dose, theta, moments,
                                       regularization) {
  gradient <- gradient_rtmo_mean(model, dose, theta)
  mu2 <- moments$central[["mu2"]]
  z <- c(1, 0, 3 * mu2)
  nuisance <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 3L, byrow = TRUE)
  derivative <- cbind(outer(z, gradient), nuisance)
  covariance <- regularization$regularized_covariance
  crossprod(derivative, solve(covariance, derivative))
}

target_rtmo_matrix <- function(model, theta, target, target_region, dose_range) {
  standardized <- if (is.null(target_region)) {
    target$standardized_dose
  } else {
    seq(target_region[[1L]], target_region[[2L]], length.out = 81L)
  }
  if (!is.null(target_region)) {
    standardized <- (standardized - dose_range[[1L]]) / diff(dose_range)
  }
  parameter_dimension <- length(theta)
  result <- matrix(0, parameter_dimension + 2L, parameter_dimension + 2L)
  slopes <- numeric(length(standardized))
  for (i in seq_along(standardized)) {
    dose <- standardized[[i]]
    slope <- dose_derivative_rtmo_mean(model, dose, theta)
    if (!is.finite(slope) || abs(slope) <= 1e-10) {
      stop("The fitted dose-response slope is too small in the target region.",
           call. = FALSE)
    }
    target_gradient <- c(gradient_rtmo_mean(model, dose, theta) / slope, 0, 0)
    result <- result + tcrossprod(target_gradient)
    slopes[[i]] <- slope
  }
  list(
    matrix = result / length(standardized),
    standardized_doses = standardized,
    slopes = slopes
  )
}

rtmo_criterion <- function(weights, information, target_matrix) {
  matrix <- Reduce(`+`, Map(function(weight, item) weight * item,
                            weights, information))
  inverse <- tryCatch(solve(matrix), error = function(e) NULL)
  if (is.null(inverse)) return(Inf)
  value <- sum(diag(inverse %*% target_matrix))
  if (!is.finite(value) || value <= 0) Inf else value
}

optimize_rtmo_weights <- function(information, target_matrix, initial = NULL) {
  count <- length(information)
  if (is.null(initial)) initial <- rep(1 / count, count)
  initial <- pmax(initial, 1e-8)
  initial <- initial / sum(initial)
  if (count == 1L) return(list(weights = 1, criterion = rtmo_criterion(1, information, target_matrix)))

  softmax <- function(parameters) {
    logits <- c(parameters, 0)
    shifted <- logits - max(logits)
    values <- exp(shifted)
    values / sum(values)
  }
  objective <- function(parameters) {
    rtmo_criterion(softmax(parameters), information, target_matrix)
  }
  gradient <- function(parameters) {
    weights <- softmax(parameters)
    matrix <- Reduce(`+`, Map(function(weight, item) weight * item,
                              weights, information))
    inverse <- tryCatch(solve(matrix), error = function(e) NULL)
    if (is.null(inverse)) return(rep(1e8, length(parameters)))
    kernel <- inverse %*% target_matrix %*% inverse
    derivatives <- -vapply(information, function(item) sum(diag(item %*% kernel)),
                           numeric(1L))
    average <- sum(weights * derivatives)
    weights[-count] * (derivatives[-count] - average)
  }
  start <- log(initial[-count] / initial[[count]])
  fit <- stats::optim(start, objective, gradient, method = "BFGS",
                      control = list(maxit = 5000, reltol = 1e-11))
  weights <- softmax(fit$par)
  # Multiplicative polishing prevents a softmax optimizer from trapping a
  # scientifically important support point at an effectively zero weight.
  for (iteration in seq_len(4000L)) {
    matrix <- Reduce(`+`, Map(function(weight, item) weight * item,
                              weights, information))
    inverse <- tryCatch(solve(matrix), error = function(e) NULL)
    if (is.null(inverse)) break
    criterion <- sum(diag(inverse %*% target_matrix))
    kernel <- inverse %*% target_matrix %*% inverse
    sensitivity <- vapply(
      information,
      function(item) sum(diag(item %*% kernel)) / criterion,
      numeric(1L)
    )
    updated <- weights * sqrt(pmax(sensitivity, 1e-12))
    updated <- updated / sum(updated)
    if (max(abs(updated - weights)) < 1e-11) {
      weights <- updated
      break
    }
    weights <- updated
  }
  list(
    weights = weights,
    criterion = rtmo_criterion(weights, information, target_matrix),
    convergence = fit$convergence
  )
}

rtmo_sensitivity <- function(information, design_information, target_matrix, criterion) {
  inverse <- solve(design_information)
  kernel <- inverse %*% target_matrix %*% inverse
  vapply(information, function(item) sum(diag(item %*% kernel)) / criterion,
         numeric(1L))
}

compute_rtmo_design <- function(model, theta, moments, regularization, target,
                                target_region, dose_range, control_variance,
                                total_sample_size = NULL) {
  target_system <- target_rtmo_matrix(model, theta, target, target_region, dose_range)
  grid <- seq(0, 1, length.out = 1001L)
  information <- lapply(grid, function(dose) {
    elemental_rtmo_information(
      model = model, dose = dose, theta = theta, moments = moments,
      regularization = regularization
    )
  })
  support <- unique(round(seq(1, length(grid), length.out = 7L)))
  weights <- rep(1 / length(support), length(support))

  for (iteration in seq_len(40L)) {
    optimized <- optimize_rtmo_weights(information[support], target_system$matrix, weights)
    weights <- optimized$weights
    keep <- weights > 1e-5
    if (sum(keep) >= 2L && any(!keep)) {
      candidate_value <- rtmo_criterion(
        weights[keep] / sum(weights[keep]), information[support[keep]], target_system$matrix
      )
      if (is.finite(candidate_value)) {
        support <- support[keep]
        weights <- weights[keep] / sum(weights[keep])
        optimized <- optimize_rtmo_weights(information[support], target_system$matrix, weights)
        weights <- optimized$weights
      }
    }
    design_information <- Reduce(`+`, Map(function(weight, item) weight * item,
                                           weights, information[support]))
    criterion <- rtmo_criterion(weights, information[support], target_system$matrix)
    sensitivity <- rtmo_sensitivity(information, design_information,
                                    target_system$matrix, criterion)
    next_index <- which.max(sensitivity)
    if (max(sensitivity) <= 1 + 2e-5) break
    if (next_index %in% support) break
    support <- sort(c(support, next_index))
    weights <- rep(1 / length(support), length(support))
  }

  optimized <- optimize_rtmo_weights(information[support], target_system$matrix, weights)
  weights <- optimized$weights
  keep <- weights > 5e-3
  if (sum(keep) >= 2L) {
    candidate_weights <- weights[keep] / sum(weights[keep])
    candidate_value <- rtmo_criterion(
      candidate_weights, information[support[keep]], target_system$matrix
    )
    if (is.finite(candidate_value)) {
      support <- support[keep]
      weights <- candidate_weights
    }
  }

  # A dense candidate grid can represent one continuous support point by
  # several adjacent grid points. Merge such numerical clusters before the
  # final weight optimization so the reported design is implementable.
  if (length(support) > 1L) {
    groups <- cumsum(c(TRUE, diff(support) > 30L))
    if (length(unique(groups)) < length(support)) {
      merged_support <- vapply(split(seq_along(support), groups), function(index) {
        as.integer(round(stats::weighted.mean(support[index], weights[index])))
      }, integer(1L))
      merged_weights <- vapply(split(seq_along(support), groups), function(index) {
        sum(weights[index])
      }, numeric(1L))
      merged_support <- pmax(1L, pmin(length(grid), merged_support))
      merged <- !duplicated(merged_support)
      support <- merged_support[merged]
      weights <- merged_weights[merged] / sum(merged_weights[merged])
    }
  }
  optimized <- optimize_rtmo_weights(information[support], target_system$matrix, weights)
  weights <- optimized$weights
  final_keep <- weights > 5e-3
  if (any(!final_keep) && sum(final_keep) >= 1L) {
    candidate_weights <- weights[final_keep] / sum(weights[final_keep])
    candidate_value <- rtmo_criterion(
      candidate_weights, information[support[final_keep]], target_system$matrix
    )
    if (is.finite(candidate_value)) {
      support <- support[final_keep]
      optimized <- optimize_rtmo_weights(
        information[support], target_system$matrix, candidate_weights
      )
      weights <- optimized$weights
    }
  }
  criterion <- optimized$criterion
  design_information <- Reduce(`+`, Map(function(weight, item) weight * item,
                                         weights, information[support]))
  sensitivity <- rtmo_sensitivity(information, design_information,
                                  target_system$matrix, criterion)

  mean_inverse_slope_squared <- mean(1 / target_system$slopes^2)
  control_component <- control_variance * mean_inverse_slope_squared
  control_weight <- sqrt(control_component) /
    (sqrt(control_component) + sqrt(criterion))
  full_new_weights <- (1 - control_weight) * weights
  approximate_variance <- criterion / (1 - control_weight) +
    control_component / control_weight

  approximate <- data.frame(
    component = c(rep("new_treatment", length(support)), "active_control"),
    standardized_dose = c(grid[support], NA_real_),
    dose = c(dose_range[[1L]] + diff(dose_range) * grid[support], NA_real_),
    conditional_weight = c(weights, NA_real_),
    weight = c(full_new_weights, control_weight),
    stringsAsFactors = FALSE
  )

  exact <- NULL
  if (!is.null(total_sample_size)) {
    exact <- exact_rtmo_allocation(
      total_sample_size, approximate, information[support],
      target_system$matrix, control_component, approximate_variance
    )
    approximate$count <- exact$counts
  } else {
    approximate$count <- NA_integer_
  }

  list(
    approximate = approximate,
    exact = exact,
    optimality = list(
      passed = max(sensitivity) <= 1 + 1e-3,
      maximum_sensitivity = max(sensitivity),
      equivalence_gap = max(sensitivity) - 1,
      grid = data.frame(
        standardized_dose = grid,
        dose = dose_range[[1L]] + diff(dose_range) * grid,
        sensitivity = sensitivity
      )
    ),
    treatment_criterion = criterion,
    control_component = control_component,
    integrated_variance = approximate_variance
  )
}

exact_rtmo_allocation <- function(n, approximate, treatment_information,
                                  target_matrix, control_component,
                                  approximate_variance) {
  quotas <- n * approximate$weight
  choices <- lapply(quotas, function(value) unique(c(floor(value), ceiling(value))))
  candidates <- expand.grid(choices, KEEP.OUT.ATTRS = FALSE,
                            stringsAsFactors = FALSE)
  candidates <- as.matrix(candidates)
  candidates <- candidates[rowSums(candidates) == n, , drop = FALSE]
  best_value <- Inf
  best <- NULL
  for (row in seq_len(nrow(candidates))) {
    counts <- as.integer(candidates[row, ])
    if (counts[[length(counts)]] <= 0L || any(counts < 0L)) next
    new_information <- Reduce(`+`, Map(function(count, item) count * item,
                                       counts[-length(counts)], treatment_information))
    inverse <- tryCatch(solve(new_information), error = function(e) NULL)
    if (is.null(inverse)) next
    variance <- sum(diag(inverse %*% target_matrix)) +
      control_component / counts[[length(counts)]]
    if (variance < best_value - 1e-14) {
      best_value <- variance
      best <- counts
    }
  }
  if (is.null(best)) stop("No feasible exact floor-ceiling allocation was found.", call. = FALSE)
  list(
    counts = best,
    quotas = quotas,
    maximum_quota_deviation = max(abs(best - quotas)),
    exact_scaled_variance = n * best_value,
    criterion_change_percent = 100 * (n * best_value - approximate_variance) /
      approximate_variance
  )
}
