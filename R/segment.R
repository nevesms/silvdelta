#' Assign points to trees
#'
#' Each point is assigned to the nearest stem in the horizontal plane
#' (Voronoi partition), optionally limited to `max_radius`. The partition is
#' purely geometric and identical for both epochs, so the *difference*
#' between epochs inside each tree cell is consistent even where crowns
#' interlock. Removed material is re-attributed more carefully in
#' [voxel_change()].
#'
#' @param cloud A height-normalized cloud.
#' @param stems Stem table from [detect_stems()].
#' @param max_radius Maximum horizontal distance (m) to the stem; farther
#'   points get `tree_id = NA`.
#' @return A copy of `cloud` with the columns `tree_id` and `r_stem`
#'   (horizontal distance to the stem axis, m).
#' @export
segment_trees <- function(cloud, stems, max_radius = Inf) {
  cloud <- data.table::copy(cloud)
  if (nrow(stems) == 0L) stop("No stems to segment with.")
  nn <- RANN::nn2(as.matrix(stems[, .(X, Y)]), as.matrix(cloud[, .(X, Y)]), k = 1)
  cloud[, `:=`(tree_id = stems$tree_id[nn$nn.idx[, 1]], r_stem = nn$nn.dists[, 1])]
  cloud[r_stem > max_radius, tree_id := NA]
  cloud[]
}

#' @rdname segment_trees
#' @param ep An [epochs()] object; both epochs are segmented with the same
#'   stems.
#' @export
segment_epochs <- function(ep, stems, max_radius = Inf) {
  .check_epochs(ep, normalized = TRUE)
  ep$t0 <- segment_trees(ep$t0, stems, max_radius)
  ep$t1 <- segment_trees(ep$t1, stems, max_radius)
  ep
}

#' Structural metrics per tree for one epoch
#'
#' @param cloud A segmented, height-normalized cloud (see [segment_trees()]).
#' @param stems Stem table from [detect_stems()].
#' @param res Voxel size (m) for the crown volume.
#' @param bin Height bin (m) of the vertical profile used for the crown base.
#' @param stem_buffer Points closer than stem radius + `stem_buffer` to the
#'   stem axis are treated as stem and excluded from crown metrics.
#' @param min_frac A bin belongs to the crown when its occupied volume is at
#'   least `min_frac` of the fullest bin of that tree.
#' @return A `data.table` per tree with `height` (99.9th percentile, m),
#'   `n_points`, `crown_volume` (occupied voxel volume, m3), `cbh` (crown
#'   base height: lowest of two consecutive crown bins, m; it is a geometric
#'   crown base and includes dead branches) and `crown_radius` (95th
#'   percentile of the distance to the stem, m).
#' @export
tree_metrics <- function(cloud, stems, res = 0.1, bin = 0.25, stem_buffer = 0.15,
                         min_frac = 0.15) {
  if (!"tree_id" %in% names(cloud)) cloud <- segment_trees(cloud, stems)
  cl <- cloud[!is.na(tree_id)]
  rad <- ifelse(is.na(stems$dbh), 0.1, stems$dbh / 2) + stem_buffer
  cl[, rs := rad[match(tree_id, stems$tree_id)]]
  H <- cl[, .(height = quantile(Z, 0.999, names = FALSE), n_points = .N), by = tree_id]
  cr <- cl[r_stem > rs & Z > 0.3]
  out <- merge(stems[, .(tree_id)], H, by = "tree_id", all.x = TRUE, sort = FALSE)
  if (nrow(cr) > 0L) {
    v <- unique(cr[, .(tree_id, i = floor(X / res), j = floor(Y / res),
                       k = floor(Z / res))])
    vol <- v[, .(crown_volume = .N * res^3), by = tree_id]
    pr <- v[, .N, by = .(tree_id, zb = floor((k + 0.5) * res / bin))]
    data.table::setorder(pr, tree_id, zb)
    pr[, thr := min_frac * max(N), by = tree_id]
    pr[, ok := N >= thr]
    cb <- pr[, {
      nxt <- (zb + 1L) %in% zb[ok]
      w <- which(ok & nxt)
      .(cbh = if (length(w)) zb[w[1]] * bin else NA_real_)
    }, by = tree_id]
    rad95 <- cr[, .(crown_radius = quantile(r_stem, 0.95, names = FALSE)), by = tree_id]
    out <- merge(out, vol, by = "tree_id", all.x = TRUE, sort = FALSE)
    out <- merge(out, cb, by = "tree_id", all.x = TRUE, sort = FALSE)
    out <- merge(out, rad95, by = "tree_id", all.x = TRUE, sort = FALSE)
  } else {
    out[, `:=`(crown_volume = NA_real_, cbh = NA_real_, crown_radius = NA_real_)]
  }
  out[is.na(crown_volume) & !is.na(height), crown_volume := 0]
  out[match(stems$tree_id, tree_id)]
}
