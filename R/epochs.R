#' Coerce a point cloud to the internal representation
#'
#' Accepts a file path (`.las`/`.laz` through the suggested package lidR, or
#' any delimited text file with X, Y, Z in the first three columns), a lidR
#' `LAS` object, a `data.frame`/`data.table` or a numeric matrix.
#'
#' @param x Point cloud (see above).
#' @param intensity Name or position of the intensity column. By default a
#'   column whose name contains "intens" (any case) is used when present,
#'   as in LAS files and CloudCompare exports.
#' @return A `data.table` with numeric columns `X`, `Y`, `Z` and, when
#'   available, `Intensity` (other columns are dropped). The input is never
#'   modified.
#' @export
as_cloud <- function(x, intensity = NULL) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("File not found: ", x)
    ext <- tolower(tools::file_ext(x))
    if (ext %in% c("las", "laz")) {
      if (!requireNamespace("lidR", quietly = TRUE))
        stop("Reading LAS/LAZ files requires the 'lidR' package.")
      x <- lidR::readLAS(x, select = "xyzi")
    } else {
      x <- data.table::fread(file = x)
    }
  }
  if (inherits(x, "LAS")) x <- methods::slot(x, "data")
  if (is.matrix(x)) x <- data.table::as.data.table(x)
  if (!is.data.frame(x)) stop("Unsupported point cloud input.")
  nm <- toupper(names(x))
  pos <- match(c("X", "Y", "Z"), nm)
  if (anyNA(pos)) {
    if (ncol(x) < 3L) stop("A point cloud needs at least three columns (X, Y, Z).")
    pos <- 1:3
  }
  out <- data.table::data.table(
    X = as.numeric(x[[pos[1]]]),
    Y = as.numeric(x[[pos[2]]]),
    Z = as.numeric(x[[pos[3]]])
  )
  ic <- if (is.null(intensity)) grep("intens", names(x), ignore.case = TRUE)[1]
        else if (is.numeric(intensity)) intensity else match(intensity, names(x))
  if (length(ic) && !is.na(ic)) out[, Intensity := as.numeric(x[[ic]])]
  out[is.finite(X) & is.finite(Y) & is.finite(Z)]
}

#' Pair two point clouds as epochs t0 and t1
#'
#' @param t0,t1 Point clouds before and after the intervention. Anything
#'   accepted by [as_cloud()].
#' @param intensity Intensity column, see [as_cloud()].
#' @return An object of class `epochs`: a list with the two clouds (`t0`,
#'   `t1`) and the processing state (`registration`, `ground`, `normalized`).
#' @examples
#' sim <- simulate_stand(nx = 2, ny = 2, seed = 1)
#' sim
#' @export
epochs <- function(t0, t1, intensity = NULL) {
  if (inherits(t0, "epochs")) return(t0)
  structure(
    list(t0 = as_cloud(t0, intensity), t1 = as_cloud(t1, intensity),
         registration = NULL, ground = NULL, normalized = FALSE),
    class = "epochs"
  )
}

#' @rdname epochs
#' @export
read_epochs <- function(t0, t1, intensity = NULL) epochs(t0, t1, intensity)

#' @export
print.epochs <- function(x, ...) {
  cat("<epochs> bi-temporal point cloud\n")
  for (e in c("t0", "t1")) {
    d <- x[[e]]
    cat(sprintf("  %s: %s points | X [%.1f, %.1f] Y [%.1f, %.1f] Z [%.1f, %.1f]\n",
                e, format(nrow(d), big.mark = ","),
                min(d$X), max(d$X), min(d$Y), max(d$Y), min(d$Z), max(d$Z)))
  }
  if (!is.null(x$registration))
    cat(sprintf("  registered: yes, %s (residual RMSE %.3f m)\n",
                x$registration$method, x$registration$rmse))
  else cat("  registered: no\n")
  cat("  normalized:", if (isTRUE(x$normalized)) "yes" else "no", "\n")
  invisible(x)
}

.check_epochs <- function(ep, normalized = FALSE) {
  if (!inherits(ep, "epochs")) stop("'ep' must be an 'epochs' object; see epochs().")
  if (normalized && !isTRUE(ep$normalized))
    warning("Epochs are not height-normalized; run normalize_epochs() first ",
            "(or set ep$normalized <- TRUE if Z already is height above ground).",
            call. = FALSE)
  invisible(TRUE)
}
