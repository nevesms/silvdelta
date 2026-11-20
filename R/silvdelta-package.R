#' silvdelta: structural change detection in forests from bi-temporal point clouds
#'
#' The package is organised as a sequence of steps that can be run one by one
#' or through the wrapper [silvdelta()]:
#'
#' 1. [epochs()] / [read_epochs()]: pair the t0 and t1 clouds.
#' 2. [register_targets()]: automatic registration on high-intensity
#'    (retro-reflective) targets, and/or [register_epochs()]: fine
#'    co-registration on stable structures. [target_plot()] turns the
#'    targets into the plot polygon, area and local frame.
#' 3. [normalize_epochs()]: heights above a terrain model common to both epochs.
#' 4. [detect_stems()]: stem positions and DBH, with or without seeds.
#' 5. Individual-tree level: [segment_epochs()], [voxel_change()],
#'    [tree_change()], [classify_intervention()].
#' 6. Area level: [match_trees()], [stand_change()], [grid_change()],
#'    [vertical_profile()].
#' 7. [classify_operation()]: which operation was carried out (pruning,
#'    thinning and its type, clearcut).
#' 8. [calibrate_biomass()]: link removed volume to field-weighed biomass.
#'
#' @keywords internal
#' @import data.table
#' @importFrom utils globalVariables
#' @importFrom stats quantile median lm coef predict rnorm runif cor
"_PACKAGE"

utils::globalVariables(c(
  ".", ".N", ":=", "X", "Y", "Z", "Zabs", "i", "j", "k", "z", "N", "d", "d1", "d2",
  "tree_id", "r_stem", "status", "component", "size", "zb", "ci", "cj", "nz",
  "cx", "cy", "cl", "dbh", "found", "present", "occupancy", "n_pts", "id",
  "height", "n_points", "crown_volume", "cbh", "crown_radius", "rs", "thr",
  "lost", "stable", "ratio", "ok", "lost_volume", "gained_volume", "t0_volume",
  "lost_frac", "prune_height", "present_t1", "change", "gi", "gj", "chm",
  "cover", "volume", "chm_t0", "chm_t1", "cover_t0", "cover_t1", "volume_t0",
  "volume_t1", "d_chm", "d_cover", "d_volume", "rel_volume", "z_low",
  "id_t1", "dist", "dbh_t0", "dbh_t1", "part", "h", "pruned", "removed",
  "prune_h", "len", "tree", "whole", "si", "sj", "nb", "n", "idx", "tid",
  "new_tree", "height_t0", "occupancy_t1", "L", "zw", "az", "inc", "tt",
  "rr", "jit", "base", "bx", "by", "bz", "rw", "branch",
  "Intensity", "target_id", "intensity", "id_t0", "residual", "cl2", "row",
  "volume_t0", "volume_t1", "operation", "thinning_type", "in_plot"
))
