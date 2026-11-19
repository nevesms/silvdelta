#' Plot area from the cloud footprint
#'
#' @param cloud A point cloud with X and Y.
#' @return Area (ha) of the convex hull of the points.
#' @export
plot_area <- function(cloud) {
  h <- grDevices::chull(cloud$X, cloud$Y)
  x <- cloud$X[h]; y <- cloud$Y[h]
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2 / 1e4
}

#' Stand-level summary of a thinning
#'
#' Computes stand density and basal area before and after the intervention
#' and the usual descriptors of a thinning.
#'
#' @param trees Either the output of [match_trees()] or a classified
#'   [tree_change()] table.
#' @param area_ha Plot area (ha); see [plot_area()].
#' @return A one-row `data.table`: number of trees (`n_t0`, `n_t1`,
#'   `n_removed`), density (`N_ha_*`, trees/ha), basal area (`G_ha_*`,
#'   m2/ha), quadratic mean diameter (`dg_*`, cm), thinning intensity in
#'   number of trees and in basal area (`intensity_N`, `intensity_G`, %)
#'   and `thinning_ratio` (dg removed / dg before; < 1 indicates thinning
#'   from below, > 1 from above). Trees without a DBH estimate count for
#'   density but not for basal area (`n_no_dbh`).
#' @export
stand_change <- function(trees, area_ha) {
  tr <- data.table::as.data.table(trees)
  if ("status" %in% names(tr)) {
    tr <- tr[status != "new"]
    d <- tr$dbh_t0; rem <- tr$status == "removed"
  } else if ("change" %in% names(tr)) {
    d <- tr$dbh; rem <- tr$change == "removed"
  } else stop("'trees' must come from match_trees() or classify_intervention().")
  g <- function(x) sum(pi * (x / 2)^2, na.rm = TRUE)
  dg <- function(x) if (all(is.na(x))) NA_real_ else 100 * sqrt(mean(x^2, na.rm = TRUE))
  n0 <- length(d); nr <- sum(rem)
  data.table::data.table(
    n_t0 = n0, n_t1 = n0 - nr, n_removed = nr,
    N_ha_t0 = n0 / area_ha, N_ha_t1 = (n0 - nr) / area_ha,
    G_ha_t0 = g(d) / area_ha, G_ha_t1 = g(d[!rem]) / area_ha,
    dg_t0 = dg(d), dg_t1 = dg(d[!rem]), dg_removed = dg(d[rem]),
    intensity_N = 100 * nr / n0,
    intensity_G = 100 * g(d[rem]) / g(d),
    thinning_ratio = dg(d[rem]) / dg(d),
    n_no_dbh = sum(is.na(d))
  )
}

#' Area-based change on a regular grid
#'
#' Wall-to-wall description of the change that needs no tree segmentation:
#' canopy height, canopy cover and occupied volume per grid cell in each
#' epoch and their differences.
#'
#' @param ep A registered, normalized [epochs()] object.
#' @param res Grid cell size (m).
#' @param hcut Height (m) above which a sub-cell counts as canopy cover.
#' @param voxel Voxel size (m) for the occupied volume.
#' @param zmin Minimum height (m) for the occupied volume.
#' @return A `data.table` with cell centres `X`, `Y`, and for each epoch
#'   `chm_*` (99th height percentile, m), `cover_*` (fraction of 16
#'   sub-cells with canopy) and `volume_*` (occupied voxel volume, m3),
#'   plus `d_chm`, `d_cover`, `d_volume` (t1 - t0) and `rel_volume`
#'   (relative change of occupied volume).
#' @export
grid_change <- function(ep, res = 1, hcut = 2, voxel = 0.1, zmin = 0.3) {
  .check_epochs(ep, normalized = TRUE)
  ox <- min(ep$t0$X, ep$t1$X); oy <- min(ep$t0$Y, ep$t1$Y)
  one <- function(d) {
    d <- d[, .(X, Y, Z, gi = floor((X - ox) / res), gj = floor((Y - oy) / res))]
    a <- d[, .(chm = quantile(Z, 0.99, names = FALSE)), by = .(gi, gj)]
    sub <- unique(d[Z > hcut, .(gi, gj, i = floor((X - ox) / (res / 4)),
                                j = floor((Y - oy) / (res / 4)))])
    cv <- sub[, .(cover = .N / 16), by = .(gi, gj)]
    v <- unique(d[Z >= zmin, .(gi, gj, i = floor(X / voxel), j = floor(Y / voxel),
                               k = floor(Z / voxel))])
    vl <- v[, .(volume = .N * voxel^3), by = .(gi, gj)]
    a <- merge(a, cv, by = c("gi", "gj"), all.x = TRUE)
    a <- merge(a, vl, by = c("gi", "gj"), all.x = TRUE)
    a[is.na(cover), cover := 0]; a[is.na(volume), volume := 0]
    a
  }
  g <- merge(one(ep$t0), one(ep$t1), by = c("gi", "gj"), all = TRUE,
             suffixes = c("_t0", "_t1"))
  g[, `:=`(X = ox + (gi + 0.5) * res, Y = oy + (gj + 0.5) * res,
           d_chm = chm_t1 - chm_t0, d_cover = cover_t1 - cover_t0,
           d_volume = volume_t1 - volume_t0)]
  g[, rel_volume := ifelse(volume_t0 > 0, d_volume / volume_t0, NA_real_)]
  g[, `:=`(gi = NULL, gj = NULL)]
  data.table::setcolorder(g, c("X", "Y"))
  g[]
}

#' Vertical profile of occupied volume in both epochs
#'
#' Plot-level vertical distribution of the occupied voxel volume, the
#' area-based signature of an intervention (pruning empties the lower
#' strata; thinning reduces all strata).
#'
#' @param ep A registered, normalized [epochs()] object.
#' @param bin Height bin (m).
#' @param voxel Voxel size (m).
#' @param zmin Minimum height (m).
#' @return A `data.table` with `z_low`, `volume_t0`, `volume_t1`,
#'   `d_volume` and `rel_volume`.
#' @export
vertical_profile <- function(ep, bin = 0.5, voxel = 0.1, zmin = 0.3) {
  .check_epochs(ep, normalized = TRUE)
  one <- function(d) {
    v <- unique(d[Z >= zmin, .(i = floor(X / voxel), j = floor(Y / voxel),
                               k = floor(Z / voxel))])
    v[, .(volume = .N * voxel^3), by = .(z_low = floor((k + 0.5) * voxel / bin) * bin)]
  }
  p <- merge(one(ep$t0), one(ep$t1), by = "z_low", all = TRUE,
             suffixes = c("_t0", "_t1"))
  p[is.na(volume_t0), volume_t0 := 0]; p[is.na(volume_t1), volume_t1 := 0]
  p[, d_volume := volume_t1 - volume_t0]
  p[, rel_volume := ifelse(volume_t0 > 0, d_volume / volume_t0, NA_real_)]
  data.table::setorder(p, z_low)
  p[]
}
