print.summary_rtmo_design <- function(x, ...) {
  cat("rtmodesign input summary\n\n")
  cat("Status: ", x$status, "\n", sep = "")
  cat("Model: ", x$model, "\n", sep = "")
  cat("Target dose: ", format(signif(x$target$dose, 6)), "\n", sep = "")
  cat("Target status: ", x$target$status, "\n", sep = "")
  cat("Residual skewness: ", format(signif(x$moments$skewness, 5)), "\n", sep = "")
  cat("Moment regularization: ",
      if (x$regularization$applied) "Applied" else "Not required", "\n", sep = "")
  cat("Regularized condition number: ",
      format(signif(x$regularization$regularized_condition_number, 5)), "\n", sep = "")
  cat("Optimality status: ", if (x$optimality$passed) "Passed" else "Warning", "\n", sep = "")
  cat("Maximum sensitivity: ",
      format(signif(x$optimality$maximum_sensitivity, 6)), "\n", sep = "")
  if (!is.null(x$exact_allocation)) {
    cat("Exact allocation available for n = ",
        x$design_conditions$total_sample_size, "\n", sep = "")
  }
  cat("Analysis observations: ", x$data_summary$analysis_n, "\n", sep = "")
  cat("New-treatment observations: ", x$data_summary$new_treatment_n, "\n", sep = "")
  cat("Active-control observations: ", x$data_summary$active_control_n, "\n", sep = "")
  invisible(x)
}
