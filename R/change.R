#' Change per tree between epochs
#'
#' Combines the per-epoch tree metrics, the stem presence test in t1 and the
#' voxel change attributed to each tree.
#'
#' The pruning height is read from the vertical pattern of the change: for
#' each tree, the fraction of off-stem t0 voxels that were lost is computed
#' per height bin; pruning removes branches from the bottom up, so
#' `prune_height` is the top of the continuous run of bins, starting near
#' the crown base, in which most voxels were lost. Trees without such a run
#' get `NA`, which is what separates pruning from scattered losses caused by
#' occlusion.
#'
#' @param ep A registered, normalized [epochs()] object.
#' @param stems Stem table of t0 from [detect_stems()].
#' @param vc A [voxel_change()] result computed with the same `stems`
#'   (computed with defaults if `NULL`).
#' @param bin Height bin (m) for the pruning height.
#' @param stem_buffer Distance (m) beyond the stem radius excluded from
#'   crown computations.
#' @param min_ratio Minimum lost fraction of a bin to belong to the pruned
#'   zone.
#' @param presence Named list of arguments passed to [stem_presence()].
#' @return A `data.table` per tree with position, `dbh`, `height_t0/t1`,
#'   `cbh_t0/t1`, `crown_volume_t0/t1`, `present_t1`, `occupancy_t1`,
#'   `lost_volume`, `gained_volume` (m3 of occupied voxels), `lost_frac`
#'   (lost / t0 occupied volume) and `prune_height` (m).
#' @export
tree_change <- function(ep, stems, vc = NULL, bin = 0.25, stem_buffer = 0.15,
                        min_ratio = 0.5, presence = list()) {
  .check_epochs(ep, normalized = TRUE)
  if (is.null(vc)) vc <- voxel_change(ep, stems = stems)
  if (!"tree_id" %in% names(vc$voxels))
    stop("'vc' has no tree attribution; run voxel_change(ep, stems = stems).")
  if (!"tree_id" %in% names(ep$t0)) ep <- segment_epochs(ep, stems)
  m0 <- tree_metrics(ep$t0, stems, stem_buffer = stem_buffer)
  m1 <- tree_metrics(ep$t1, stems, stem_buffer = stem_buffer)
  pr <- do.call(stem_presence, c(list(cloud = ep$t1, stems = stems), presence))
  vv <- vc$res^3
  V <- vc$voxels[!is.na(tree_id) & status != "noise"]
  agg <- V[, .(lost_volume = sum(status == "lost") * vv,
               gained_volume = sum(status == "gained") * vv,
               t0_volume = sum(status != "gained") * vv), by = tree_id]
  rad <- ifelse(is.na(stems$dbh), 0.1, stems$dbh / 2) + stem_buffer
  V[, rs := rad[match(tree_id, stems$tree_id)]]
  P <- V[r_stem > rs & status != "gained",
         .(lost = sum(status == "lost"), n = .N), by = .(tree_id, zb = floor(Z / bin))]
  P <- P[n >= 3L]
  data.table::setorder(P, tree_id, zb)
  ph <- P[, {
    okb <- lost / n >= min_ratio
    first <- which(okb)[1]
    if (is.na(first) || first > 4L) .(prune_height = NA_real_)
    else {
      last <- first
      while (last < .N && okb[last + 1L]) last <- last + 1L
      .(prune_height = (zb[last] + 1) * bin)
    }
  }, by = tree_id]
  out <- stems[, .(tree_id, X, Y, dbh, found)]
  out[, `:=`(height_t0 = m0$height, height_t1 = m1$height,
             cbh_t0 = m0$cbh, cbh_t1 = m1$cbh,
             crown_volume_t0 = m0$crown_volume, crown_volume_t1 = m1$crown_volume,
             present_t1 = pr$present, occupancy_t1 = pr$occupancy)]
  out <- merge(out, agg, by = "tree_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, ph, by = "tree_id", all.x = TRUE, sort = FALSE)
  for (cc in c("lost_volume", "gained_volume"))
    data.table::set(out, which(is.na(out[[cc]])), cc, 0)
  out[, lost_frac := ifelse(is.na(t0_volume) | t0_volume == 0, NA_real_,
                            lost_volume / t0_volume)]
  out[, t0_volume := NULL]
  out[match(stems$tree_id, tree_id)]
}

#' Classify the intervention applied to each tree
#'
#' Transparent rules, in this order: a tree whose stem is absent in t1 was
#' `"removed"` (thinning/harvest); a standing tree with a bottom-up pruned
#' zone and a relevant loss was `"pruned"`; the others are `"unchanged"`.
#'
#' @param tc A [tree_change()] table.
#' @param min_lost_frac Minimum fraction of the t0 occupied volume lost to
#'   call a tree pruned.
#' @return `tc` with the column `change`; `prune_height` is set to `NA` for
#'   trees that were not pruned.
#' @export
classify_intervention <- function(tc, min_lost_frac = 0.03) {
  tc <- data.table::copy(tc)
  tc[, change := ifelse(!present_t1, "removed",
                 ifelse(!is.na(prune_height) & !is.na(lost_frac) &
                          lost_frac >= min_lost_frac, "pruned", "unchanged"))]
  tc[change != "pruned", prune_height := NA_real_]
  tc[]
}
