\name{rtmo_design}
\alias{rtmo_design}
\title{Construct a Regularized Third-Order Moment Optimal Design}
\usage{
rtmo_design(formula, data, model, active_control, dose_range,
  n = NULL, target_region = NULL)
}
\arguments{
  \item{formula}{A formula of the form \code{outcome ~ dose | treatment}.}
  \item{data}{A patient-level data frame in long format.}
  \item{model}{A supported dose-response model name.}
  \item{active_control}{The treatment value identifying the active control.}
  \item{dose_range}{Two increasing allowed dose values in original units.}
  \item{n}{Optional future total sample size.}
  \item{target_region}{Optional two-value target-dose planning interval. If
  omitted, the package uses the estimated target plus or minus 15 percent of
  the allowed dose range, truncated to the dose boundaries.}
}
\value{An object of class \code{rtmo_design}.}
\description{
Fits the target-dose problem from patient-level pilot data, estimates and
regularizes the required high-order moments, constructs an approximate RTMO
design, checks its optimality certificate, and optionally returns an exact
integer allocation.
}
\examples{
data(rtmo_example)
fit <- rtmo_design(
  improvement ~ dose | treatment,
  data = rtmo_example,
  model = "mm",
  active_control = "control",
  dose_range = c(0, 100),
  n = 240
)
fit
}
