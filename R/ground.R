#' Simple gridded terrain model
#'
#' Builds a terrain model from the lowest points of each grid cell. Cells
#' without ground returns (e.g. occupied by a stem) are detected as steps
#' relative to their neighbours, removed and filled, and the surface is
#' smoothed. It is intentionally simple and dependency-free; it suits plot
#' clouds from terrestrial/mobile scanning, where the ground is densely
#' sampled. For difficult terrain, normalize the clouds with a dedicated
#' ground classifier beforehand and skip [normalize_epochs()].
#'
#' @param cloud Point cloud ([as_cloud()] input or `data.table` with X, Y, Z).
#' @param res Cell size (m).
#' @param q Low quantile of Z taken as ground level in each cell (robust to
#'   a few points below the ground).
#' @param max_step Maximum height step (m) of a cell relative to the median
#'   of its neighbours before it is treated as non-ground.
#' @return A list with the matrix `z`, the origin (`ox`, `oy`) and `res`.
#'   Use [ground_height()] to interpolate.
#' @export
ground_model <- function(cloud, res = 0.5, q = 0.01, max_step = 0.3) {
  if (!data.table::is.data.table(cloud)) cloud <- as_cloud(cloud)
  ox <- min(cloud$X); oy <- min(cloud$Y)
  g <- data.table::data.table(i = floor((cloud$X - ox) / res) + 1L,
                              j = floor((cloud$Y - oy) / res) + 1L,
                              Z = cloud$Z)
  g <- g[, .(z = quantile(Z, q, names = FALSE)), by = .(i, j)]
  ni <- max(2L, max(g$i)); nj <- max(2L, max(g$j))
  M <- matrix(NA_real_, ni, nj)
  M[cbind(g$i, g$j)] <- g$z
  for (it in 1:3) {
    med <- .focal(M, "median")
    bad <- !is.na(M) & !is.na(med) & (M - med > max_step | med - M > 2 * max_step)
    if (!any(bad)) break
    M[bad] <- NA_real_
  }
  guard <- 0L
  while (anyNA(M) && guard < 200L) {
    f <- .focal(M, "mean")
    fill <- is.na(M) & !is.na(f)
    if (!any(fill)) break
    M[fill] <- f[fill]
    guard <- guard + 1L
  }
  if (anyNA(M)) M[is.na(M)] <- median(M, na.rm = TRUE)
  M <- .focal(M, "mean", include_self = TRUE)
  list(z = M, ox = ox, oy = oy, res = res)
}

#' Terrain height at given positions
#'
#' Bilinear interpolation of a [ground_model()].
#'
#' @param g A terrain model from [ground_model()].
#' @param X,Y Coordinates.
#' @return Numeric vector of terrain elevations.
#' @export
ground_height <- function(g, X, Y) {
  ni <- nrow(g$z); nj <- ncol(g$z)
  fx <- pmin(pmax((X - g$ox) / g$res - 0.5, 0), ni - 1)
  fy <- pmin(pmax((Y - g$oy) / g$res - 0.5, 0), nj - 1)
  i0 <- pmin(floor(fx), ni - 2); j0 <- pmin(floor(fy), nj - 2)
  tx <- fx - i0; ty <- fy - j0
  z00 <- g$z[cbind(i0 + 1, j0 + 1)]; z10 <- g$z[cbind(i0 + 2, j0 + 1)]
  z01 <- g$z[cbind(i0 + 1, j0 + 2)]; z11 <- g$z[cbind(i0 + 2, j0 + 2)]
  (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 +
    (1 - tx) * ty * z01 + tx * ty * z11
}

.focal <- function(M, fun = c("mean", "median"), include_self = FALSE) {
  fun <- match.arg(fun)
  ni <- nrow(M); nj <- ncol(M)
  P <- matrix(NA_real_, ni + 2L, nj + 2L)
  P[2:(ni + 1L), 2:(nj + 1L)] <- M
  sh <- expand.grid(di = -1:1, dj = -1:1)
  if (!include_self) sh <- sh[!(sh$di == 0 & sh$dj == 0), ]
  A <- vapply(seq_len(nrow(sh)), function(s)
    as.vector(P[(2:(ni + 1L)) + sh$di[s], (2:(nj + 1L)) + sh$dj[s]]),
    numeric(ni * nj))
  out <- if (fun == "mean") rowMeans(A, na.rm = TRUE)
         else apply(A, 1L, median, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  matrix(out, ni, nj)
}

#' Normalize both epochs with a common terrain model
#'
#' Heights of t0 and t1 are computed against the *same* terrain model so that
#' no artificial vertical change is introduced by two independent ground
#' estimates. Run [register_epochs()] first.
#'
#' @param ep An [epochs()] object.
#' @param res,q,max_step Passed to [ground_model()].
#' @param source Which epoch(s) feed the terrain model: `"both"` (default,
#'   lowest surface seen in either epoch), `"t0"` or `"t1"` (useful when
#'   harvest residues cover the ground in one epoch).
#' @return The `epochs` object with `Z` as height above ground, the original
#'   elevation kept in `Zabs`, and the model stored in `ep$ground`.
#' @export
normalize_epochs <- function(ep, res = 0.5, q = 0.01, max_step = 0.3,
                             source = c("both", "t0", "t1")) {
  .check_epochs(ep)
  source <- match.arg(source)
  if (isTRUE(ep$normalized)) {
    warning("Epochs are already normalized; returning unchanged.", call. = FALSE)
    return(ep)
  }
  src <- switch(source, both = rbind(ep$t0[, .(X, Y, Z)], ep$t1[, .(X, Y, Z)]),
                t0 = ep$t0, t1 = ep$t1)
  g <- ground_model(src, res = res, q = q, max_step = max_step)
  for (e in c("t0", "t1")) {
    d <- data.table::copy(ep[[e]])
    d[, Zabs := Z]
    d[, Z := Zabs - ground_height(g, X, Y)]
    ep[[e]] <- d
  }
  ep$ground <- g
  ep$normalized <- TRUE
  ep
}
