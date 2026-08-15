library(rtmodesign)

data(rtmo_example)

fit <- rtmo_design(
  improvement ~ dose | treatment,
  data = rtmo_example,
  model = "mm",
  active_control = "control",
  dose_range = c(0, 100),
  n = 240
)

dir.create(file.path("man", "figures"), recursive = TRUE, showWarnings = FALSE)

rtmo_theme <- function(mar = c(4.5, 4.8, 3.2, 1.2)) {
  par(
    bg = "white",
    fg = "#263238",
    col.axis = "#455A64",
    col.lab = "#263238",
    col.main = "#18262C",
    family = "sans",
    font.main = 2,
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.9,
    las = 1,
    lend = "round",
    ljoin = "round",
    mar = mar,
    mgp = c(2.6, 0.75, 0),
    tcl = -0.25,
    xaxs = "r",
    yaxs = "r"
  )
}

png(
  file.path("man", "figures", "rtmodesign-overview.png"),
  width = 2200,
  height = 1000,
  res = 200,
  type = "cairo-png"
)
par(mfrow = c(1, 2))
rtmo_theme(c(4.5, 4.8, 3.2, 2.2))
plot(fit, type = "dose_response", main = "Active-control-matched target")
abline(v = fit$target$planning_region, col = "#9AA7AD", lty = 3, lwd = 1.5)
abline(v = fit$target$dose, col = "#D55E00", lty = 2, lwd = 1.5)
legend(
  "bottomright",
  legend = c("Fitted dose-response", "Matched target", "Planning limits"),
  col = c("#0072B2", "#D55E00", "#9AA7AD"),
  lty = c(1, 2, 3),
  lwd = c(2, 1.5, 1.5),
  bty = "n",
  cex = 0.84
)
rtmo_theme(c(4.5, 5.2, 3.2, 1.2))
plot(fit, type = "design", main = "Recommended allocation", ylim = c(0, 0.6))
text(
  x = c(0.7, 1.9, 3.1),
  y = fit$approximate_design$weight + 0.025,
  labels = paste0("n = ", fit$approximate_design$count),
  cex = 0.82,
  col = "#455A64"
)
mtext("Approximate weights; exact total n = 240", side = 3, line = 0.35,
      cex = 0.78, col = "#607D8B")
dev.off()

png(
  file.path("man", "figures", "rtmodesign-diagnostics.png"),
  width = 2200,
  height = 1000,
  res = 200,
  type = "cairo-png"
)
par(mfrow = c(1, 2))
rtmo_theme(c(4.5, 4.8, 3.2, 4.0))
plot(fit, type = "data", main = "Observed new-treatment responses")
rtmo_theme(c(4.5, 6.2, 3.2, 1.2))
plot(fit, type = "diagnostic", main = "Optimality certificate")
legend(
  "bottomright",
  legend = c("Sensitivity", "Optimality bound"),
  col = c("#D55E00", "#555555"),
  lty = c(1, 2),
  lwd = c(2, 1),
  bty = "n",
  cex = 0.84
)
dev.off()

print(fit)
