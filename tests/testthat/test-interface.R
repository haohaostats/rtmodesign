test_that("ordinary-user interface validates a supported design", {
  data(rtmo_example)
  fit <- rtmo_design(
    improvement ~ dose | treatment,
    data = rtmo_example,
    model = "mm",
    active_control = "control",
    dose_range = c(0, 100),
    n = 240
  )

  expect_s3_class(fit, "rtmo_design")
  expect_identical(fit$status, "complete")
  expect_equal(fit$data_summary$analysis_n, 150)
  expect_equal(fit$data_summary$active_control_n, 30)
  expect_equal(fit$design_conditions$total_sample_size, 240L)
  expect_true(is.finite(fit$target$dose))
  expect_true(fit$target$dose > 0 && fit$target$dose < 100)
  expect_length(fit$moments$central, 5)
  expect_true(is.finite(fit$moments$condition_number))
  expect_lte(fit$regularization$regularized_condition_number, 100 + 1e-8)
  expect_true(fit$regularization$lambda >= 0 && fit$regularization$lambda <= 1)
  expect_true(fit$optimality$passed)
  expect_equal(sum(fit$approximate_design$weight), 1, tolerance = 1e-8)
  expect_equal(sum(fit$approximate_design$count), 240)
  expect_lte(fit$exact_allocation$maximum_quota_deviation, 1)
})

test_that("formula must follow the documented contract", {
  data(rtmo_example)
  expect_error(
    rtmo_design(improvement ~ dose + treatment, rtmo_example, "mm",
                "control", c(0, 100)),
    "outcome ~ dose \\| treatment"
  )
})

test_that("unknown active control is rejected", {
  data(rtmo_example)
  expect_error(
    rtmo_design(improvement ~ dose | treatment, rtmo_example, "mm",
                "standard", c(0, 100)),
    "was not found"
  )
})

test_that("target region must lie inside dose range", {
  data(rtmo_example)
  expect_error(
    rtmo_design(improvement ~ dose | treatment, rtmo_example, "mm",
                "control", c(0, 100), target_region = c(20, 120)),
    "must lie inside"
  )
})

test_that("all advertised built-in models complete the ordinary workflow", {
  data(rtmo_example)
  for (model in c("linear", "mm", "emax", "sigmoid_emax", "exponential")) {
    fit <- rtmo_design(
      improvement ~ dose | treatment,
      data = rtmo_example,
      model = model,
      active_control = "control",
      dose_range = c(0, 100)
    )
    expect_identical(fit$status, "complete", info = model)
    expect_true(fit$optimality$passed, info = model)
    expect_equal(sum(fit$approximate_design$weight), 1, tolerance = 1e-8,
                 info = model)
  }
})
