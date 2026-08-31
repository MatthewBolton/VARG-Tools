args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: test_stratigraphic_loo.R <app-directory>", call. = FALSE)
}

app_directory <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
source(file.path(app_directory, "R", "functions", "sc_helpers.R"))

ties <- data.frame(
  id = 1:5,
  ref_sample = paste0("R", 1:5),
  ref_z = c(0, 10, 25, 30, 40),
  target_sample = paste0("T", 1:5),
  target_z = c(0, 10, 20, 30, 40),
  use_in_warp = TRUE,
  stringsAsFactors = FALSE
)

fit <- sc_calculate_warp(ties, method = "monotonic")
loo <- fit$diagnostics$loo_residuals

stopifnot(
  identical(fit$model_type, "Monotonic piecewise-linear"),
  identical(fit$diagnostics$diagnostic_type, "leave-one-out"),
  identical(
    loo$prediction_type,
    c(
      "Endpoint extrapolation",
      "Interior interpolation",
      "Interior interpolation",
      "Interior interpolation",
      "Endpoint extrapolation"
    )
  ),
  all(loo$status == "Estimated"),
  isTRUE(all.equal(loo$residual[loo$id == 3], 5, tolerance = 1e-12)),
  is.finite(fit$diagnostics$loo_rmse),
  is.finite(fit$diagnostics$loo_mae),
  isTRUE(all.equal(fit$diagnostics$anchor_rmse, 0, tolerance = 0))
)

cat("stratigraphic_loo_ok\n")
