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
  # Pairwise simplex polishing can place genuinely inactive support points
  # exactly on the boundary, which a softmax parameterization cannot do.
  if (count > 1L) {
    for (cycle in seq_len(20L)) {
      previous <- weights
      for (i in seq_len(count - 1L)) {
        for (j in (i + 1L):count) {
          lower <- -weights[[j]]
          upper <- weights[[i]]
          if (upper - lower <= 1e-14) next
          along_edge <- function(transfer) {
            candidate <- weights
            candidate[[i]] <- candidate[[i]] - transfer
            candidate[[j]] <- candidate[[j]] + transfer
            rtmo_criterion(candidate, information, target_matrix)
          }
          edge_derivative <- function(transfer) {
            candidate <- weights
            candidate[[i]] <- candidate[[i]] - transfer
            candidate[[j]] <- candidate[[j]] + transfer
            matrix <- Reduce(`+`, Map(function(weight, item) weight * item,
                                      candidate, information))
            inverse <- tryCatch(solve(matrix), error = function(e) NULL)
            if (is.null(inverse)) return(NA_real_)
            criterion <- sum(diag(inverse %*% target_matrix))
            kernel <- inverse %*% target_matrix %*% inverse
            sensitivity_i <- sum(diag(information[[i]] %*% kernel)) / criterion
            sensitivity_j <- sum(diag(information[[j]] %*% kernel)) / criterion
            criterion * (sensitivity_i - sensitivity_j)
          }
          epsilon <- min(1e-9, (upper - lower) / 100)
          scan <- seq(lower + epsilon, upper - epsilon, length.out = 21L)
          derivative <- vapply(scan, edge_derivative, numeric(1L))
          changes <- which(
            is.finite(derivative[-length(derivative)]) &
              is.finite(derivative[-1L]) &
              derivative[-length(derivative)] * derivative[-1L] <= 0
          )
          root <- if (length(changes)) {
            bracket <- scan[c(changes[[1L]], changes[[1L]] + 1L)]
            stats::uniroot(
              edge_derivative, interval = bracket, tol = 1e-15,
              maxiter = 1000L
            )$root
          } else {
            stats::optimize(along_edge, interval = c(lower, upper), tol = 1e-12)$minimum
          }
          locations <- c(lower, root, upper)
          values <- vapply(locations, along_edge, numeric(1L))
          transfer <- locations[[which.min(values)]]
          weights[[i]] <- weights[[i]] - transfer
          weights[[j]] <- weights[[j]] + transfer
          weights[weights < 1e-13] <- 0
          weights <- weights / sum(weights)
        }
      }
      if (max(abs(weights - previous)) < 1e-12) break
    }
  }
  final_criterion <- rtmo_criterion(weights, information, target_matrix)
  final_matrix <- Reduce(`+`, Map(function(weight, item) weight * item,
                                  weights, information))
  final_inverse <- solve(final_matrix)
  final_kernel <- final_inverse %*% target_matrix %*% final_inverse
  final_sensitivity <- vapply(
    information,
    function(item) sum(diag(item %*% final_kernel)) / final_criterion,
    numeric(1L)
  )
  active <- weights > 1e-10
  kkt_residual <- max(
    if (any(active)) max(abs(final_sensitivity[active] - 1)) else 0,
    if (any(!active)) max(pmax(final_sensitivity[!active] - 1, 0)) else 0
  )
  list(
    weights = weights,
    criterion = final_criterion,
    convergence = fit$convergence,
    kkt_residual = kkt_residual
  )
}

rtmo_sensitivity <- function(information, design_information, target_matrix, criterion) {
  inverse <- solve(design_information)
  kernel <- inverse %*% target_matrix %*% inverse
  vapply(information, function(item) sum(diag(item %*% kernel)) / criterion,
         numeric(1L))
}

maximize_rtmo_sensitivity <- function(model, theta, moments, regularization,
                                      design_information, target_matrix,
                                      criterion, grid, grid_sensitivity) {
  inverse <- solve(design_information)
  kernel <- inverse %*% target_matrix %*% inverse
  at_dose <- function(dose) {
    item <- elemental_rtmo_information(
      model = model, dose = dose, theta = theta, moments = moments,
      regularization = regularization
    )
    sum(diag(item %*% kernel)) / criterion
  }

  # Every smooth interior maximizer is bracketed by a discrete local maximum.
  # Refining both adjacent cells also covers maxima close to a grid boundary.
  local <- which(
    grid_sensitivity[2:(length(grid) - 1L)] >=
      grid_sensitivity[1:(length(grid) - 2L)] &
      grid_sensitivity[2:(length(grid) - 1L)] >=
      grid_sensitivity[3:length(grid)]
  ) + 1L
  candidates <- unique(c(1L, length(grid), local, which.max(grid_sensitivity)))
  best_index <- which.max(grid_sensitivity)
  best_dose <- grid[[best_index]]
  best_value <- grid_sensitivity[[best_index]]

  for (index in candidates) {
    lower <- grid[[max(1L, index - 1L)]]
    upper <- grid[[min(length(grid), index + 1L)]]
    if (upper <= lower) next
    refined <- stats::optimize(
      function(dose) -at_dose(dose), interval = c(lower, upper),
      tol = 1e-12
    )
    value <- -refined$objective
    if (is.finite(value) && value > best_value) {
      best_value <- value
      best_dose <- refined$minimum
    }
  }
  list(dose = best_dose, sensitivity = best_value)
}

compute_rtmo_design <- function(model, theta, moments, regularization, target,
                                target_region, dose_range, control_variance,
                                total_sample_size = NULL) {
  target_system <- target_rtmo_matrix(model, theta, target, target_region, dose_range)
  grid <- seq(0, 1, length.out = 10001L)
  information <- lapply(grid, function(dose) {
    elemental_rtmo_information(
      model = model, dose = dose, theta = theta, moments = moments,
      regularization = regularization
    )
  })
  # The manuscript model has two mean parameters and therefore starts at
  # {0, 1/2, 1}. Models with more mean parameters need enough distinct
  # positive doses for a nonsingular initial information matrix.
  support_doses <- seq(0, 1, length.out = max(3L, length(theta) + 1L))
  weights <- rep(1 / length(support_doses), length(support_doses))

  for (iteration in seq_len(100L)) {
    support_information <- lapply(support_doses, function(dose) {
      elemental_rtmo_information(
        model = model, dose = dose, theta = theta, moments = moments,
        regularization = regularization
      )
    })
    optimized <- optimize_rtmo_weights(
      support_information, target_system$matrix, weights
    )
    weights <- optimized$weights
    if (length(support_doses) > 1L) {
      groups <- cumsum(c(TRUE, diff(support_doses) > 2e-2))
      if (length(unique(groups)) < length(support_doses)) {
        positions <- split(seq_along(support_doses), groups)
        support_doses <- vapply(positions, function(index) {
          stats::weighted.mean(support_doses[index], weights[index])
        }, numeric(1L))
        weights <- vapply(positions, function(index) sum(weights[index]), numeric(1L))
        weights <- weights / sum(weights)
        support_information <- lapply(support_doses, function(dose) {
          elemental_rtmo_information(
            model = model, dose = dose, theta = theta, moments = moments,
            regularization = regularization
          )
        })
        optimized <- optimize_rtmo_weights(
          support_information, target_system$matrix, weights
        )
        weights <- optimized$weights
      }
    }
    keep <- weights >= 1e-10
    if (sum(keep) >= 2L && any(!keep)) {
      candidate_value <- rtmo_criterion(
        weights[keep] / sum(weights[keep]), support_information[keep],
        target_system$matrix
      )
      if (is.finite(candidate_value)) {
        support_doses <- support_doses[keep]
        support_information <- support_information[keep]
        weights <- weights[keep] / sum(weights[keep])
        optimized <- optimize_rtmo_weights(
          support_information, target_system$matrix, weights
        )
        weights <- optimized$weights
      }
    }
    if (length(support_doses) > 1L) {
      for (location_cycle in seq_len(5L)) {
        for (index in which(support_doses > 0 & support_doses < 1)) {
          lower <- if (index == 1L) 0 else
            (support_doses[[index - 1L]] + support_doses[[index]]) / 2
          upper <- if (index == length(support_doses)) 1 else
            (support_doses[[index]] + support_doses[[index + 1L]]) / 2
          location_fit <- stats::optimize(function(dose) {
            candidate_information <- support_information
            candidate_information[[index]] <- elemental_rtmo_information(
              model = model, dose = dose, theta = theta, moments = moments,
              regularization = regularization
            )
            rtmo_criterion(weights, candidate_information, target_system$matrix)
          }, interval = c(lower, upper), tol = 1e-12)
          support_doses[[index]] <- location_fit$minimum
          support_information[[index]] <- elemental_rtmo_information(
            model = model, dose = location_fit$minimum, theta = theta,
            moments = moments, regularization = regularization
          )
        }
        optimized <- optimize_rtmo_weights(
          support_information, target_system$matrix, weights
        )
        weights <- optimized$weights
      }
    }
    design_information <- Reduce(`+`, Map(function(weight, item) weight * item,
                                           weights, support_information))
    criterion <- rtmo_criterion(weights, support_information, target_system$matrix)
    sensitivity <- rtmo_sensitivity(information, design_information,
                                    target_system$matrix, criterion)
    maximum <- maximize_rtmo_sensitivity(
      model, theta, moments, regularization, design_information,
      target_system$matrix, criterion, grid, sensitivity
    )
    if (maximum$sensitivity <= 1 + 1e-8) break
    nearest <- which.min(abs(support_doses - maximum$dose))
    if (abs(support_doses[[nearest]] - maximum$dose) <= 2e-2) {
      movement <- abs(support_doses[[nearest]] - maximum$dose)
      support_doses[[nearest]] <- maximum$dose
      ordering <- order(support_doses)
      support_doses <- support_doses[ordering]
      weights <- weights[ordering]
      if (movement <= 1e-12) break
    } else {
      old_doses <- support_doses
      old_weights <- weights
      support_doses <- sort(c(old_doses, maximum$dose))
      weights <- numeric(length(support_doses))
      is_new <- which.min(abs(support_doses - maximum$dose))
      weights[[is_new]] <- 0.01
      weights[-is_new] <- 0.99 * old_weights[order(old_doses)]
    }
  }

  support_information <- lapply(support_doses, function(dose) {
    elemental_rtmo_information(
      model = model, dose = dose, theta = theta, moments = moments,
      regularization = regularization
    )
  })
  optimized <- optimize_rtmo_weights(
    support_information, target_system$matrix, weights
  )
  weights <- optimized$weights
  keep <- weights >= 1e-10
  if (sum(keep) >= 2L) {
    candidate_weights <- weights[keep] / sum(weights[keep])
    candidate_value <- rtmo_criterion(
      candidate_weights, support_information[keep], target_system$matrix
    )
    if (is.finite(candidate_value)) {
      support_doses <- support_doses[keep]
      support_information <- support_information[keep]
      weights <- candidate_weights
    }
  }
  optimized <- optimize_rtmo_weights(
    support_information, target_system$matrix, weights
  )
  weights <- optimized$weights
  final_keep <- weights >= 1e-10
  if (any(!final_keep) && sum(final_keep) >= 1L) {
    candidate_weights <- weights[final_keep] / sum(weights[final_keep])
    candidate_value <- rtmo_criterion(
      candidate_weights, support_information[final_keep], target_system$matrix
    )
    if (is.finite(candidate_value)) {
      support_doses <- support_doses[final_keep]
      support_information <- support_information[final_keep]
      optimized <- optimize_rtmo_weights(
        support_information, target_system$matrix, candidate_weights
      )
      weights <- optimized$weights
    }
  }
  design_information <- Reduce(`+`, Map(function(weight, item) weight * item,
                                         weights, support_information))
  # Recompute from the final matrix so the sensitivity normalization and the
  # reported criterion use exactly the same floating-point factorization on
  # every supported platform.
  criterion <- rtmo_criterion(
    weights, support_information, target_system$matrix
  )
  sensitivity <- rtmo_sensitivity(information, design_information,
                                  target_system$matrix, criterion)
  maximum <- maximize_rtmo_sensitivity(
    model, theta, moments, regularization, design_information,
    target_system$matrix, criterion, grid, sensitivity
  )
  # The weighted average sensitivity is one analytically. Guard against a
  # sub-unit maximum caused only by floating-point inversion error.
  maximum$sensitivity <- max(1, maximum$sensitivity)

  mean_inverse_slope_squared <- mean(1 / target_system$slopes^2)
  control_component <- control_variance * mean_inverse_slope_squared
  control_weight <- sqrt(control_component) /
    (sqrt(control_component) + sqrt(criterion))
  full_new_weights <- (1 - control_weight) * weights
  approximate_variance <- criterion / (1 - control_weight) +
    control_component / control_weight

  approximate <- data.frame(
    component = c(rep("new_treatment", length(support_doses)), "active_control"),
    standardized_dose = c(support_doses, NA_real_),
    dose = c(dose_range[[1L]] + diff(dose_range) * support_doses, NA_real_),
    conditional_weight = c(weights, NA_real_),
    weight = c(full_new_weights, control_weight),
    stringsAsFactors = FALSE
  )

  exact <- NULL
  if (!is.null(total_sample_size)) {
    exact <- exact_rtmo_allocation(
      total_sample_size, approximate, support_information,
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
      passed = maximum$sensitivity <= 1 + 1e-8,
      maximum_sensitivity = maximum$sensitivity,
      maximum_sensitivity_dose = dose_range[[1L]] +
        diff(dose_range) * maximum$dose,
      equivalence_gap = maximum$sensitivity - 1,
      certified_efficiency_lower_bound = 1 / maximum$sensitivity,
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
