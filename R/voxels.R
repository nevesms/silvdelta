#' Voxel-based change between epochs
#'
#' Both epochs are voxelized on the same grid. A voxel occupied in t0 is
#' `"lost"` when no occupied t1 voxel lies within `tol`; a voxel of t1 is
#' `"gained"` when no t0 voxel lies within `tol`; the others are `"stable"`.
#' Using a distance tolerance instead of a strict voxel-to-voxel comparison
#' absorbs the residual registration error and the sampling differences
#' between surveys. Connected groups of lost (or gained) voxels smaller than
#' `min_component` are relabelled `"noise"`.
#'
#' When `stems` are given, voxels are attributed to trees. Lost voxels are
#' attributed by *connectivity*: a connected group of lost voxels (a removed
#' branch) that touches one single stem (within `attach_dist`) is assigned
#' entirely to that tree, even where its tip crosses into a neighbour's
#' space. Groups touching several stems, or none, fall back to the nearest
#' stem voxel by voxel.
#'
#' Note that `"gained"` voxels between two surveys a few days apart are
#' mostly an occlusion effect (parts not seen in t0). They are reported so
#' that the asymmetry of occlusion between epochs can be evaluated.
#'
#' @param ep A registered and normalized [epochs()] object.
#' @param res Voxel size (m).
#' @param tol Distance tolerance (m). Default: the larger of one voxel and
#'   twice the registration RMSE.
#' @param zmin,zmax Height range (m) analysed (excludes ground and litter).
#' @param min_component Minimum size (voxels) of a connected changed group.
#' @param stems Optional stem table from [detect_stems()].
#' @param attach_dist Horizontal distance (m) from the stem axis within
#'   which a lost group is considered attached to that stem.
#' @param max_radius Maximum distance (m) to a stem for attribution.
#' @return An object of class `voxel_change`: a list with `voxels` (a
#'   `data.table` with voxel centres `X`, `Y`, `Z`, `status`, `component`
#'   and, with stems, `tree_id` and `r_stem`), `res`, `tol` and `zmin`.
#' @export
voxel_change <- function(ep, res = 0.05, tol = NULL, zmin = 0.3, zmax = Inf,
                         min_component = 5L, stems = NULL, attach_dist = 0.5,
                         max_radius = Inf) {
  .check_epochs(ep, normalized = TRUE)
  a <- ep$t0[Z >= zmin & Z <= zmax]; b <- ep$t1[Z >= zmin & Z <= zmax]
  ox <- min(a$X, b$X); oy <- min(a$Y, b$Y)
  vox <- function(d) unique(data.table::data.table(
    i = floor((d$X - ox) / res), j = floor((d$Y - oy) / res),
    k = floor((d$Z - zmin) / res)))
  v0 <- vox(a); v1 <- vox(b)
  ctr <- function(v) cbind(ox + (v$i + 0.5) * res, oy + (v$j + 0.5) * res,
                           zmin + (v$k + 0.5) * res)
  c0 <- ctr(v0); c1 <- ctr(v1)
  if (is.null(tol)) {
    rm <- ep$registration$rmse
    tol <- max(1.01 * res, if (is.null(rm)) 0 else 2 * rm)
  }
  if (nrow(c0) == 0L) stop("No t0 points in the analysed height range.")
  d01 <- if (nrow(c1)) RANN::nn2(c1, c0, k = 1)$nn.dists[, 1] else rep(Inf, nrow(c0))
  d10 <- if (nrow(c1)) RANN::nn2(c0, c1, k = 1)$nn.dists[, 1] else numeric(0)
  V <- rbind(
    data.table::data.table(X = c0[, 1], Y = c0[, 2], Z = c0[, 3],
                           status = ifelse(d01 > tol, "lost", "stable")),
    data.table::data.table(X = c1[, 1], Y = c1[, 2], Z = c1[, 3],
                           status = "gained")[d10 > tol]
  )
  V[, component := NA_integer_]
  off <- 0L
  for (st in c("lost", "gained")) {
    w <- which(V$status == st)
    if (!length(w)) next
    cl <- dbscan::dbscan(as.matrix(V[w, .(X, Y, Z)]), eps = 1.75 * res,
                         minPts = 1)$cluster
    sz <- tabulate(cl)
    small <- sz[cl] < min_component
    data.table::set(V, w, "component", cl + off)
    data.table::set(V, w[small], "status", "noise")
    off <- off + max(cl)
  }
  if (!is.null(stems) && nrow(stems) > 0L) {
    sxy <- as.matrix(stems[, .(X, Y)])
    nn <- RANN::nn2(sxy, as.matrix(V[, .(X, Y)]), k = 1)
    V[, `:=`(si = nn$nn.idx[, 1], r_stem = nn$nn.dists[, 1])]
    L <- V[status == "lost"]
    if (nrow(L) > 0L) {
      m <- L[, .(d = min(r_stem)), by = .(component, si)]
      data.table::setorder(m, component, d)
      best <- m[, .(sj = si[1], d1 = d[1], d2 = if (.N > 1L) d[2] else Inf),
                by = component]
      best <- best[d1 <= attach_dist & d2 > attach_dist]
      if (nrow(best) > 0L) {
        w <- which(V$status == "lost" & V$component %in% best$component)
        new_si <- best$sj[match(V$component[w], best$component)]
        data.table::set(V, w, "si", new_si)
        data.table::set(V, w, "r_stem",
                        sqrt((V$X[w] - sxy[new_si, 1])^2 + (V$Y[w] - sxy[new_si, 2])^2))
      }
    }
    V[, tree_id := stems$tree_id[si]]
    V[r_stem > max_radius, tree_id := NA]
    V[, si := NULL]
  }
  structure(list(voxels = V[], res = res, tol = tol, zmin = zmin),
            class = "voxel_change")
}

#' @export
print.voxel_change <- function(x, ...) {
  tb <- table(x$voxels$status)
  vv <- x$res^3
  cat(sprintf("<voxel_change> voxel %.3f m, tolerance %.3f m\n", x$res, x$tol))
  for (s in names(tb))
    cat(sprintf("  %-7s %9s voxels  %10.3f m3\n", s, format(tb[[s]], big.mark = ","),
                tb[[s]] * vv))
  invisible(x)
}

#' Vertical distribution of lost and gained volume
#'
#' @param vc A [voxel_change()] result.
#' @param bin Height bin (m).
#' @param by_tree Summarise per tree (requires stems in [voxel_change()]).
#' @return A `data.table` with `z_low` (lower limit of the stratum, m) and
#'   the occupied voxel volume (m3) `lost`, `gained` and `stable`.
#' @export
change_profile <- function(vc, bin = 0.5, by_tree = FALSE) {
  V <- vc$voxels[status != "noise"]
  vv <- vc$res^3
  V[, z_low := floor(Z / bin) * bin]
  keys <- if (by_tree && "tree_id" %in% names(V)) c("tree_id", "z_low") else "z_low"
  out <- V[, .(lost = sum(status == "lost") * vv,
               gained = sum(status == "gained") * vv,
               stable = sum(status == "stable") * vv), by = keys]
  data.table::setorderv(out, keys)
  out[]
}
