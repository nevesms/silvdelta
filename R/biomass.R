#' Calibrate removed biomass against field measurements
#'
#' The voxel volume lost by a tree is a geometric quantity; converting it to
#' mass needs a calibration with material weighed in the field. This fits a
#' simple linear model and reports its performance by leave-one-out
#' cross-validation, which is the honest figure for small samples.
#'
#' @param volume Lost volume per tree (m3), e.g. `lost_volume` from
#'   [tree_change()].
#' @param mass Removed biomass per tree weighed in the field (kg).
#' @param intercept Fit an intercept? The default (`FALSE`) forces zero mass
#'   at zero lost volume.
#' @return An object of class `biomass_calibration` with the fitted `model`,
#'   the leave-one-out predictions (`loo`) and the cross-validated `rmse`
#'   (kg), `rmse_pct` (% of the mean), `bias` (kg) and `r2`. Use
#'   `predict()` to apply it to new volumes.
#' @export
calibrate_biomass <- function(volume, mass, intercept = FALSE) {
  d <- data.frame(volume = volume, mass = mass)
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 4L) stop("At least four trees with volume and mass are needed.")
  fm <- if (intercept) mass ~ volume else mass ~ 0 + volume
  fit <- lm(fm, data = d)
  loo <- vapply(seq_len(nrow(d)), function(r)
    predict(lm(fm, data = d[-r, ]), newdata = d[r, , drop = FALSE]), numeric(1))
  err <- loo - d$mass
  structure(list(
    model = fit, data = d, loo = loo,
    rmse = sqrt(mean(err^2)),
    rmse_pct = 100 * sqrt(mean(err^2)) / mean(d$mass),
    bias = mean(err),
    r2 = 1 - sum(err^2) / sum((d$mass - mean(d$mass))^2)
  ), class = "biomass_calibration")
}

#' @export
print.biomass_calibration <- function(x, ...) {
  cat("<biomass_calibration> n =", nrow(x$data), "\n")
  cf <- coef(x$model)
  cat("  coefficients:", paste(names(cf), signif(cf, 4), sep = " = ", collapse = "; "), "\n")
  cat(sprintf("  leave-one-out: RMSE %.2f kg (%.1f%%), bias %.2f kg, R2 %.2f\n",
              x$rmse, x$rmse_pct, x$bias, x$r2))
  invisible(x)
}

#' @export
predict.biomass_calibration <- function(object, volume, ...) {
  as.numeric(predict(object$model, newdata = data.frame(volume = volume)))
}
