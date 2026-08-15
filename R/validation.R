supported_rtmo_models <- function() {
  c("linear", "mm", "emax", "sigmoid_emax", "exponential")
}

validate_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be one non-missing character value.", name), call. = FALSE)
  }
  x
}

validate_dose_range <- function(dose_range) {
  if (!is.numeric(dose_range) || length(dose_range) != 2L ||
      any(!is.finite(dose_range)) || dose_range[[1L]] >= dose_range[[2L]]) {
    stop("`dose_range` must contain two finite increasing numeric values.", call. = FALSE)
  }
  as.numeric(dose_range)
}

validate_total_sample_size <- function(n) {
  if (is.null(n)) return(NULL)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 4 || abs(n - round(n)) > sqrt(.Machine$double.eps)) {
    stop("`n` must be NULL or one integer-valued number of at least 4.", call. = FALSE)
  }
  as.integer(round(n))
}

validate_target_region <- function(target_region, dose_range) {
  if (is.null(target_region)) return(NULL)
  if (!is.numeric(target_region) || length(target_region) != 2L ||
      any(!is.finite(target_region)) || target_region[[1L]] >= target_region[[2L]]) {
    stop("`target_region` must be NULL or two finite increasing dose values.", call. = FALSE)
  }
  if (target_region[[1L]] < dose_range[[1L]] ||
      target_region[[2L]] > dose_range[[2L]]) {
    stop("`target_region` must lie inside `dose_range`.", call. = FALSE)
  }
  as.numeric(target_region)
}

prepare_rtmo_data <- function(data, variables, active_control, dose_range) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  missing_columns <- setdiff(unname(variables), names(data))
  if (length(missing_columns)) {
    stop(sprintf("Columns not found in `data`: %s.",
                 paste(missing_columns, collapse = ", ")), call. = FALSE)
  }

  outcome <- data[[variables[["outcome"]]]]
  dose <- data[[variables[["dose"]]]]
  treatment <- data[[variables[["treatment"]]]]
  if (!is.numeric(outcome)) stop("The outcome column must be numeric.", call. = FALSE)
  if (!is.numeric(dose)) stop("The dose column must be numeric.", call. = FALSE)
  if (is.list(treatment)) stop("The treatment column must be atomic.", call. = FALSE)
  treatment <- as.character(treatment)

  if (!active_control %in% treatment) {
    stop(sprintf("Active-control level '%s' was not found in the treatment column.",
                 active_control), call. = FALSE)
  }
  is_control <- !is.na(treatment) & treatment == active_control
  is_new <- !is.na(treatment) & !is_control

  usable <- !is.na(outcome) & !is.na(treatment) & (is_control | !is.na(dose))
  removed <- sum(!usable)
  if (!any(usable)) stop("No complete observations remain after validation.", call. = FALSE)

  analysis_data <- data[usable, , drop = FALSE]
  outcome <- outcome[usable]
  dose <- dose[usable]
  treatment <- treatment[usable]
  is_control <- treatment == active_control
  is_new <- !is_control

  if (sum(is_control) < 2L) {
    stop("The active-control group must contain at least two usable observations.", call. = FALSE)
  }
  if (sum(is_new) < 3L) {
    stop("The new-treatment dose-response data must contain at least three usable observations.",
         call. = FALSE)
  }
  new_doses <- dose[is_new]
  distinct_doses <- sort(unique(new_doses))
  if (length(distinct_doses) < 2L) {
    stop("At least two distinct new-treatment dose levels are required.", call. = FALSE)
  }
  if (any(new_doses < dose_range[[1L]] | new_doses > dose_range[[2L]])) {
    stop("Observed new-treatment doses must lie inside `dose_range`.", call. = FALSE)
  }
  if (anyNA(outcome) || any(!is.finite(outcome))) {
    stop("Usable outcome values must be finite.", call. = FALSE)
  }
  if (any(!is.finite(new_doses))) {
    stop("Usable new-treatment doses must be finite.", call. = FALSE)
  }

  standardized_dose <- rep(NA_real_, length(dose))
  standardized_dose[is_new] <-
    (dose[is_new] - dose_range[[1L]]) / diff(dose_range)

  list(
    data = analysis_data,
    outcome = outcome,
    dose = dose,
    standardized_dose = standardized_dose,
    treatment = treatment,
    is_control = is_control,
    is_new = is_new,
    removed = removed,
    observed_doses = distinct_doses,
    control_dose_values_ignored = sum(is_control & !is.na(dose))
  )
}

