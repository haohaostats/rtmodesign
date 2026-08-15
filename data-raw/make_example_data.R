set.seed(1042)

doses <- rep(c(0, 20, 50, 100), each = 30)
mean_response <- 9 * doses / (28 + doses)
latent <- rnorm(length(doses))
error <- 1.4 * (0.82 * latent + 0.18 * (latent^2 - 1))

new_data <- data.frame(
  improvement = mean_response + error,
  dose = doses,
  treatment = "new",
  stringsAsFactors = FALSE
)

control_data <- data.frame(
  improvement = 5.4 + rnorm(30, sd = 1.4),
  dose = NA_real_,
  treatment = "control",
  stringsAsFactors = FALSE
)

rtmo_example <- rbind(new_data, control_data)
save(rtmo_example, file = file.path("data", "rtmo_example.rda"), compress = "xz")

