#' Construct a regularized third-order moment optimal design
#'
#' Validates patient-level pilot data and the clinical target-dose design
#' conditions through the stable ordinary-user interface of rtmodesign.
#'
#' @param formula A formula of the form `outcome ~ dose | treatment`.
#' @param data A patient-level data frame in long format.
#' @param model One of `"linear"`, `"mm"`, `"emax"`, `"sigmoid_emax"`, or
#'   `"exponential"`.
#' @param active_control The treatment-column value identifying the active
#'   control group.
#' @param dose_range Two increasing values giving the allowed new-treatment
#'   dose range in original units.
#' @param n Optional future total sample size for exact allocation.
#' @param target_region Optional two-value target-dose planning interval in
#'   original dose units.
#'
#' @return An object of class `rtmo_design`.
#' @export
rtmo_design <- function(formula, data, model, active_control, dose_range,
                        n = NULL, target_region = NULL) {
  call <- match.call()
  variables <- parse_rtmo_formula(formula)
  model <- match.arg(model, supported_rtmo_models())
  active_control <- validate_scalar_character(active_control, "active_control")
  dose_range <- validate_dose_range(dose_range)
  n <- validate_total_sample_size(n)
  target_region <- validate_target_region(target_region, dose_range)
  prepared <- prepare_rtmo_data(data, variables, active_control, dose_range)

  notices <- character()
  if (prepared$removed > 0L) {
    notices <- c(notices, sprintf("%d incomplete observation%s excluded.",
                                  prepared$removed,
                                  if (prepared$removed == 1L) " was" else "s were"))
  }
  if (prepared$control_dose_values_ignored > 0L) {
    notices <- c(notices, "Dose values in the active-control group were ignored.")
  }

  reference <- abs(prepared$standardized_dose) <= sqrt(.Machine$double.eps)
  if (!any(reference & prepared$is_new, na.rm = TRUE)) {
    stop("The new-treatment data must include the lower endpoint of `dose_range` for reference centering.",
         call. = FALSE)
  }
  reference_mean <- mean(prepared$outcome[reference & prepared$is_new], na.rm = TRUE)
  adjusted_outcome <- prepared$outcome - reference_mean
  model_fit <- fit_rtmo_model(
    model,
    prepared$standardized_dose[prepared$is_new],
    adjusted_outcome[prepared$is_new]
  )
  control_mean <- mean(adjusted_outcome[prepared$is_control])
  control_variance <- stats::var(prepared$outcome[prepared$is_control])
  standardized_target <- inverse_rtmo_mean(model, control_mean, model_fit$coefficients)
  if (!is.finite(standardized_target)) {
    stop("The active-control mean does not define a finite matched dose under the fitted model.",
         call. = FALSE)
  }
  target_dose <- dose_range[[1L]] + diff(dose_range) * standardized_target
  tolerance <- 1e-8 * max(1, diff(dose_range))
  if (target_dose < dose_range[[1L]] - tolerance ||
      target_dose > dose_range[[2L]] + tolerance) {
    stop("The estimated matched dose lies outside `dose_range`.", call. = FALSE)
  }
  target_status <- if (min(target_dose - dose_range[[1L]],
                           dose_range[[2L]] - target_dose) <= 0.05 * diff(dose_range)) {
    "near_boundary"
  } else {
    "interior"
  }
  if (target_status == "near_boundary") {
    notices <- c(notices, "The estimated matched dose is within 5% of a dose-range boundary.")
  }
  target_region_source <- "user"
  if (is.null(target_region)) {
    half_width <- 0.15 * diff(dose_range)
    target_region <- c(
      max(dose_range[[1L]], target_dose - half_width),
      min(dose_range[[2L]], target_dose + half_width)
    )
    target_region_source <- "automatic"
    notices <- c(
      notices,
      "The target-planning region was set automatically to the estimated target plus or minus 15% of the dose range."
    )
  }
  moments <- estimate_rtmo_moments(model_fit$residuals)
  regularization <- regularize_rtmo_moments(moments)
  if (regularization$applied) {
    notices <- c(notices, "Adaptive high-order moment regularization was applied.")
  }
  design_result <- compute_rtmo_design(
    model = model,
    theta = model_fit$coefficients,
    moments = moments,
    regularization = regularization,
    target = list(standardized_dose = standardized_target),
    target_region = target_region,
    dose_range = dose_range,
    control_variance = control_variance,
    total_sample_size = n
  )
  design_warnings <- character()
  if (!design_result$optimality$passed) {
    design_warnings <- "The computed design did not satisfy the optimality tolerance."
  }

  object <- list(
    call = call,
    interface_version = "1.0.0",
    status = if (design_result$optimality$passed) "complete" else "diagnostic_warning",
    variables = variables,
    model = c(model_fit, list(fitted = TRUE, reference_mean = reference_mean)),
    design_conditions = list(
      active_control = active_control,
      dose_range = dose_range,
      target_region = target_region,
      target_region_source = target_region_source,
      total_sample_size = n
    ),
    data_summary = list(
      supplied_n = nrow(data),
      analysis_n = nrow(prepared$data),
      excluded_n = prepared$removed,
      new_treatment_n = sum(prepared$is_new),
      active_control_n = sum(prepared$is_control),
      observed_doses = prepared$observed_doses,
      treatment_levels = sort(unique(prepared$treatment))
    ),
    analysis_data = prepared$data,
    standardized_dose = prepared$standardized_dose,
    control = list(
      adjusted_mean = control_mean,
      variance = control_variance,
      sample_size = sum(prepared$is_control)
    ),
    target = list(
      dose = target_dose,
      standardized_dose = standardized_target,
      status = target_status,
      planning_region = target_region,
      planning_region_source = target_region_source
    ),
    moments = moments,
    regularization = regularization,
    approximate_design = design_result$approximate,
    exact_allocation = design_result$exact,
    optimality = design_result$optimality,
    criterion = list(
      integrated_variance = design_result$integrated_variance,
      treatment_component = design_result$treatment_criterion,
      control_component = design_result$control_component
    ),
    notices = notices,
    warnings = design_warnings
  )
  class(object) <- "rtmo_design"
  object
}
