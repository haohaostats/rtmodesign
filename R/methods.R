#' @export
print.rtmo_design <- function(x, ...) {
  cat("Regularized Third-Order Moment Optimal Design\n\n")
  cat("Input\n")
  cat("  Outcome variable:      ", x$variables[["outcome"]], "\n", sep = "")
  cat("  Dose variable:         ", x$variables[["dose"]], "\n", sep = "")
  cat("  Treatment variable:    ", x$variables[["treatment"]], "\n", sep = "")
  cat("  Active control:        ", x$design_conditions$active_control, "\n", sep = "")
  cat("  Dose-response model:   ", x$model$name, "\n", sep = "")
  cat("  Analysis sample size:  ", x$data_summary$analysis_n, "\n", sep = "")
  cat("  Observed doses:        ",
      paste(format(x$data_summary$observed_doses), collapse = ", "), "\n", sep = "")
  cat("  Allowed dose range:    ",
      paste(format(x$design_conditions$dose_range), collapse = " to "), "\n", sep = "")
  if (!is.null(x$design_conditions$total_sample_size)) {
    cat("  Future sample size:    ", x$design_conditions$total_sample_size, "\n", sep = "")
  }
  cat("\nStatus\n")
  cat("  Input validation:      Passed\n")
  cat("  Model fitting:         Passed\n")
  cat("  Estimated target dose: ", format(signif(x$target$dose, 5)), "\n", sep = "")
  cat("  Target status:         ", x$target$status, "\n", sep = "")
  cat("  Target-planning range: ",
      paste(format(signif(x$target$planning_region, 5)), collapse = " to "),
      " (", x$target$planning_region_source, ")\n", sep = "")
  cat("  Residual skewness:     ", format(signif(x$moments$skewness, 4)), "\n", sep = "")
  cat("  Moment regularization: ",
      if (x$regularization$applied) "Applied" else "Not required", "\n", sep = "")
  cat("  Raw condition no.:     ",
      format(signif(x$regularization$raw_condition_number, 4)), "\n", sep = "")
  cat("  Final condition no.:   ",
      format(signif(x$regularization$regularized_condition_number, 4)), "\n", sep = "")
  cat("  Design computation:    Completed\n")
  cat("  Optimality status:     ",
      if (x$optimality$passed) "Passed" else "Warning", "\n", sep = "")
  cat("  Maximum sensitivity:   ",
      format(signif(x$optimality$maximum_sensitivity, 6)), "\n", sep = "")
  cat("\nRecommended allocation\n")
  for (i in seq_len(nrow(x$approximate_design))) {
    row <- x$approximate_design[i, ]
    label <- if (row$component == "active_control") {
      x$design_conditions$active_control
    } else {
      paste0("Dose ", format(signif(row$dose, 5)))
    }
    count_text <- if (is.na(row$count)) "" else paste0("; n = ", row$count)
    cat("  ", label, ": weight = ", format(round(row$weight, 4), nsmall = 4),
        count_text, "\n", sep = "")
  }
  if (length(x$notices)) {
    cat("\nNotices\n")
    for (notice in x$notices) cat("  - ", notice, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.rtmo_design <- function(object, ...) {
  structure(
    list(
      call = object$call,
      status = object$status,
      variables = object$variables,
      model = object$model$name,
      coefficients = object$model$coefficients,
      design_conditions = object$design_conditions,
      data_summary = object$data_summary,
      target = object$target,
      control = object$control,
      moments = object$moments,
      regularization = object$regularization,
      approximate_design = object$approximate_design,
      exact_allocation = object$exact_allocation,
      optimality = object$optimality,
      criterion = object$criterion,
      notices = object$notices,
      warnings = object$warnings
    ),
    class = "summary_rtmo_design"
  )
}

#' @export
as.data.frame.rtmo_design <- function(x, row.names = NULL, optional = FALSE,
                                      what = c("allocation", "data_summary"), ...) {
  what <- match.arg(what)
  if (what == "allocation") {
    if (is.null(x$approximate_design)) {
      return(data.frame(
        component = character(), dose = numeric(), weight = numeric(),
        count = integer(), stringsAsFactors = FALSE
      ))
    }
    return(as.data.frame(x$approximate_design, stringsAsFactors = FALSE))
  }
  data.frame(
    supplied_n = x$data_summary$supplied_n,
    analysis_n = x$data_summary$analysis_n,
    excluded_n = x$data_summary$excluded_n,
    new_treatment_n = x$data_summary$new_treatment_n,
    active_control_n = x$data_summary$active_control_n
  )
}

#' @export
plot.rtmo_design <- function(x, type = c("data", "dose_response", "design", "diagnostic"), ...) {
  type <- match.arg(type)
  if (!type %in% c("data", "dose_response") && is.null(x$approximate_design)) {
    stop("This plot becomes available after the computational engine has fitted the design.",
         call. = FALSE)
  }
  if (type == "data") {
    outcome <- x$analysis_data[[x$variables[["outcome"]]]]
    dose <- x$analysis_data[[x$variables[["dose"]]]]
    treatment <- as.character(x$analysis_data[[x$variables[["treatment"]]]])
    keep <- treatment != x$design_conditions$active_control
    graphics::plot(dose[keep], outcome[keep],
                   xlab = "Dose", ylab = x$variables[["outcome"]],
                   pch = 19, col = "#0072B2", bty = "l", ...)
    return(invisible(x))
  }
  if (type == "dose_response") {
    dose_range <- x$design_conditions$dose_range
    original_grid <- seq(dose_range[[1L]], dose_range[[2L]], length.out = 301L)
    standardized_grid <- (original_grid - dose_range[[1L]]) / diff(dose_range)
    fitted <- predict_rtmo_mean(x$model$name, standardized_grid, x$model$coefficients)
    graphics::plot(original_grid, fitted, type = "l", lwd = 2, col = "#0072B2",
                   xlab = "Dose", ylab = "Reference-adjusted mean response",
                   bty = "l", ...)
    graphics::points(x$target$dose, x$control$adjusted_mean,
                     pch = 19, col = "#D55E00")
    return(invisible(x))
  }
  if (type == "design") {
    labels <- ifelse(
      x$approximate_design$component == "active_control",
      x$design_conditions$active_control,
      paste0("Dose ", format(signif(x$approximate_design$dose, 4)))
    )
    graphics::barplot(x$approximate_design$weight, names.arg = labels,
                      col = "#0072B2", border = NA, bty = "l",
                      ylab = "Allocation weight", ...)
    return(invisible(x))
  }
  if (type == "diagnostic") {
    diagnostic <- x$optimality$grid
    graphics::plot(diagnostic$dose, diagnostic$sensitivity, type = "l",
                   lwd = 2, col = "#D55E00", bty = "l",
                   xlab = "Dose", ylab = "Normalized sensitivity", ...)
    graphics::abline(h = 1, lty = 2, col = "#555555")
    return(invisible(x))
  }
  stop("Requested design plot is not yet implemented.", call. = FALSE)
}
