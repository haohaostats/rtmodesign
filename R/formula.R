# Internal formula parser for outcome ~ dose | treatment.
parse_rtmo_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula of the form outcome ~ dose | treatment.",
         call. = FALSE)
  }

  outcome <- formula[[2L]]
  rhs <- formula[[3L]]
  if (!is.symbol(outcome) || !is.call(rhs) || !identical(rhs[[1L]], as.name("|")) ||
      length(rhs) != 3L || !is.symbol(rhs[[2L]]) || !is.symbol(rhs[[3L]])) {
    stop("`formula` must use untransformed column names: outcome ~ dose | treatment.",
         call. = FALSE)
  }

  variables <- c(
    outcome = as.character(outcome),
    dose = as.character(rhs[[2L]]),
    treatment = as.character(rhs[[3L]])
  )
  if (anyDuplicated(unname(variables))) {
    stop("Outcome, dose, and treatment must refer to three different columns.",
         call. = FALSE)
  }
  variables
}

