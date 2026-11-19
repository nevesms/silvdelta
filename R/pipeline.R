#' Run the complete bi-temporal change pipeline
#'
#' Wrapper that chains the steps of the package. Every step can also be run
#' on its own for full control; see the package help page.
#'
#' Two levels of analysis are available:
#'
#' * `level = "tree"` (individual tree): stems are located in t0 (from
#'   `seeds` or automatically), both epochs are partitioned by tree, the
#'   voxel change is attributed to trees and each tree is classified as
#'   pruned, removed or unchanged. Suited to pruning studies and to any
#'   case in which results must be linked to trees measured in the field.
#' * `level = "area"` (area-based): stems are detected independently in
#'   each epoch and matched; the output is the stand summary (density,
#'   basal area, thinning intensity and type), the gridded change and the
#'   vertical profile. Suited to thinning at stand level; no seeds needed.
#'
#' The stand summary, grid and profile are returned at both levels.
#'
#' Registration can use the intensity targets (`registration = "targets"`):
#' targets are detected, matched by their arrangement and t1 is moved onto
#' t0 without any prior alignment; the plot polygon formed by the targets
#' then clips both epochs, gives the plot area and, optionally, the local
#' plot frame. Without targets, `registration = "icp"` aligns roughly
#' pre-aligned clouds on stable structures.
#'
#' @param t0,t1 Point clouds before and after the intervention (anything
#'   accepted by [as_cloud()]), or `t0` an [epochs()] object.
#' @param seeds Optional stem seeds ([as_seeds()]); `NULL` for automatic
#'   stem detection.
#' @param level `"tree"` or `"area"`.
#' @param registration `"targets"` ([register_targets()]), `"icp"`
#'   ([register_epochs()]) or `"none"`.
#' @param normalize Normalize heights with [normalize_epochs()]? Set to
#'   `FALSE` when Z already is height above ground.
#' @param voxel Voxel size (m) for [voxel_change()].
#' @param max_radius Maximum distance (m) from a stem for tree attribution.
#' @param area_ha Plot area (ha). Default: the target polygon (with
#'   targets) or the convex hull of the t0 cloud.
#' @param grid Cell size (m) for [grid_change()].
#' @param clip Clip the epochs to the target polygon?
#' @param clip_buffer Buffer (m) added to the target polygon before clipping.
#' @param plot_frame Express the results in the local plot frame defined by
#'   the targets (see [target_plot()])?
#' @param stem_args,register_args,target_args,voxel_args,operation_args
#'   Named lists of further arguments for [detect_stems()],
#'   [register_epochs()], [register_targets()], [voxel_change()] and
#'   [classify_operation()].
#' @param verbose Print progress?
#' @return An object of class `silvdelta`: a list with `level`, `epochs`,
#'   `stems` (t0), `trees` (per-tree change; at area level the output of
#'   [match_trees()]), `voxels` (tree level only), `stand`, `grid`,
#'   `profile`, `operation` ([classify_operation()]), `plot` (target plot
#'   geometry, if any), `area_ha` and `registration`.
#' @examples
#' \donttest{
#' sim <- simulate_stand(nx = 3, ny = 2, seed = 1)
#' res <- silvdelta(sim, seeds = sim$seeds, level = "tree")
#' res
#' res$trees
#' }
#' @export
silvdelta <- function(t0, t1 = NULL, seeds = NULL, level = c("tree", "area"),
                       registration = c("icp", "targets", "none"), normalize = TRUE,
                       voxel = 0.05, max_radius = Inf, area_ha = NULL, grid = 1,
                       clip = TRUE, clip_buffer = 0, plot_frame = FALSE,
                       stem_args = list(), register_args = list(),
                       target_args = list(), voxel_args = list(),
                       operation_args = list(), verbose = TRUE) {
  level <- match.arg(level)
  registration <- match.arg(registration)
  say <- function(...) if (verbose) message(...)
  ep <- if (inherits(t0, "epochs")) t0 else epochs(t0, t1)
  ep <- structure(ep[intersect(names(ep), c("t0", "t1", "registration", "ground",
                                            "normalized", "targets"))],
                  class = "epochs")
  geom <- NULL
  if (registration != "none" && is.null(ep$registration) && !isTRUE(ep$normalized)) {
    if (registration == "targets") {
      say("[1/6] Registering on intensity targets")
      ep <- do.call(register_targets, c(list(ep = ep), target_args))
      rg <- ep$registration
      say(sprintf("      %d targets matched (RMSE %.3f m); final RMSE %.3f m%s",
                  rg$n_targets, rg$rmse_targets, rg$rmse,
                  if (isTRUE(rg$ambiguous)) " [AMBIGUOUS]" else ""))
      geom <- target_plot(ep$targets$t0)
      ep <- remove_target_points(ep, ep$targets$t0)
      if (!is.null(geom$polygon)) {
        if (clip) ep <- clip_epochs(ep, geom$polygon, clip_buffer)
        if (is.null(area_ha)) {
          pg <- .buffer_poly(geom$polygon$X, geom$polygon$Y, if (clip) clip_buffer else 0)
          area_ha <- .poly_area(pg$X, pg$Y) / 1e4
        }
        say(sprintf("      plot: %d vertices, %.4f ha%s", nrow(geom$polygon),
                    geom$area_ha, if (is.null(geom$frame)) " (no middle target)" else ""))
        if (plot_frame && !is.null(geom$frame)) {
          ep <- to_plot_frame(ep, geom$frame)
          xy <- .frame_xy(geom$polygon$X, geom$polygon$Y, geom$frame)
          geom$polygon <- data.table::data.table(X = xy$X, Y = xy$Y)
        }
      }
    } else {
      say("[1/6] Registering t1 onto t0 (stable layer)")
      ep <- do.call(register_epochs, c(list(ep = ep), register_args))
      say(sprintf("      residual RMSE %.3f m; vertical correction %.3f m",
                  ep$registration$rmse, ep$registration$dz))
    }
  }
  if (normalize && !isTRUE(ep$normalized)) {
    say("[2/6] Normalizing heights with a common terrain model")
    ep <- normalize_epochs(ep)
  } else ep$normalized <- TRUE
  if (is.null(area_ha)) area_ha <- plot_area(ep$t0)

  say("[3/6] Detecting stems in t0", if (is.null(seeds)) " (automatic)" else " (seeded)")
  stems <- do.call(detect_stems, c(list(cloud = ep$t0, seeds = seeds), stem_args))
  if (nrow(stems) == 0L) stop("No stems detected in t0; adjust 'stem_args' or give seeds.")
  say(sprintf("      %d trees (%d with DBH)", nrow(stems), sum(stems$found)))

  vc <- NULL
  if (level == "tree") {
    say("[4/6] Voxel change and attribution to trees")
    ep <- segment_epochs(ep, stems, max_radius)
    vc <- do.call(voxel_change, c(list(ep = ep, res = voxel, stems = stems,
                                       max_radius = max_radius), voxel_args))
    trees <- classify_intervention(tree_change(ep, stems, vc))
  } else {
    say("[4/6] Detecting stems in t1 and matching epochs")
    sa <- stem_args; sa$seed_radius <- NULL
    stems1 <- do.call(detect_stems, c(list(cloud = ep$t1), sa))
    trees <- match_trees(stems, stems1)
    pr <- stem_presence(ep$t1, stems)
    still <- trees$status == "removed" & !is.na(trees$tree_id) &
      pr$present[match(trees$tree_id, pr$tree_id)]
    trees[which(still), status := "kept"]
  }
  say("[5/6] Stand summary, grid and vertical profile")
  out <- list(level = level, epochs = ep, stems = stems, trees = trees, voxels = vc,
              stand = stand_change(trees, area_ha),
              grid = grid_change(ep, res = grid),
              profile = vertical_profile(ep), plot = geom,
              area_ha = area_ha, registration = ep$registration)
  say("[6/6] Identifying the operation")
  out$operation <- do.call(classify_operation, c(list(x = out), operation_args))
  structure(out, class = "silvdelta")
}

#' @export
print.silvdelta <- function(x, ...) {
  cat(sprintf("<silvdelta> %s-level analysis | plot area %.4f ha\n", x$level, x$area_ha))
  rg <- x$registration
  if (!is.null(rg)) {
    cat(sprintf("  registration: %s, RMSE %.3f m", rg$method, rg$rmse))
    if (!is.null(rg$n_targets))
      cat(sprintf(" (%d targets, %.3f m%s)", rg$n_targets, rg$rmse_targets,
                  if (isTRUE(rg$ambiguous)) ", AMBIGUOUS" else ""))
    cat("\n")
  }
  op <- x$operation
  if (!is.null(op))
    cat(sprintf("  OPERATION: %s%s\n", op$operation,
                if (is.na(op$thinning_type)) "" else paste0(" (", op$thinning_type, ")")))
  s <- x$stand
  cat(sprintf("  trees: %d in t0, %d in t1 (%d removed)\n", s$n_t0, s$n_t1, s$n_removed))
  cat(sprintf("  density: %.0f -> %.0f trees/ha | basal area: %.1f -> %.1f m2/ha\n",
              s$N_ha_t0, s$N_ha_t1, s$G_ha_t0, s$G_ha_t1))
  if (s$n_removed > 0)
    cat(sprintf("  thinning intensity: %.1f%% of trees, %.1f%% of basal area (ratio %.2f)\n",
                s$intensity_N, s$intensity_G, s$thinning_ratio))
  if (x$level == "tree") {
    tb <- table(x$trees$change)
    cat("  classification:", paste(names(tb), tb, sep = " = ", collapse = ", "), "\n")
    p <- x$trees[change == "pruned"]
    if (nrow(p))
      cat(sprintf("  pruned trees: mean pruning height %.2f m; lost volume %.3f m3 (sum)\n",
                  mean(p$prune_height), sum(p$lost_volume)))
  }
  invisible(x)
}

#' Plot a silvdelta result
#'
#' @param x A [silvdelta()] result.
#' @param type `"map"`: stem map coloured by change class, over the gridded
#'   relative change of occupied volume; `"profile"`: vertical profile of
#'   occupied volume in t0 and t1.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
plot.silvdelta <- function(x, type = c("map", "profile"), ...) {
  type <- match.arg(type)
  if (type == "profile") {
    p <- x$profile
    graphics::plot(p$volume_t0, p$z_low, type = "l", lwd = 2, col = "grey30",
                   xlim = c(0, max(p$volume_t0, p$volume_t1)),
                   xlab = "Occupied volume (m3)", ylab = "Height (m)")
    graphics::lines(p$volume_t1, p$z_low, lwd = 2, col = "firebrick")
    graphics::legend("topright", c("t0", "t1"), lwd = 2,
                     col = c("grey30", "firebrick"), bty = "n")
    return(invisible(x))
  }
  g <- x$grid
  rel <- pmin(pmax(g$rel_volume, -1), 0); rel[is.na(rel)] <- 0
  pal <- grDevices::colorRampPalette(c("white", "orange", "firebrick"))(11)
  graphics::plot(g$X, g$Y, pch = 15, cex = 2.2, col = pal[round(-rel * 10) + 1],
                 asp = 1, xlab = "X (m)", ylab = "Y (m)",
                 main = "Relative loss of occupied volume")
  tr <- x$trees
  cls <- if ("change" %in% names(tr)) tr$change else tr$status
  cols <- c(pruned = "blue", removed = "black", unchanged = "darkgreen",
            kept = "darkgreen", new = "grey50")
  graphics::points(tr$X, tr$Y, pch = ifelse(cls == "removed", 4, 19),
                   col = cols[cls], cex = 1.2, lwd = 2)
  graphics::legend("topright", legend = unique(cls), col = cols[unique(cls)],
                   pch = ifelse(unique(cls) == "removed", 4, 19), bty = "n",
                   bg = "white")
  invisible(x)
}

#' Export the tables of a result to CSV files
#'
#' @param x A [silvdelta()] result.
#' @param dir Output directory (created if needed). Existing files with the
#'   same names are not overwritten unless `overwrite = TRUE`.
#' @param prefix File name prefix.
#' @param overwrite Overwrite existing files?
#' @return The paths written, invisibly.
#' @export
export_results <- function(x, dir, prefix = "silvdelta", overwrite = FALSE) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tabs <- list(operation = x$operation, trees = x$trees, stems = x$stems,
               stand = x$stand, grid = x$grid, profile = x$profile)
  if (!is.null(x$epochs$targets)) tabs$targets <- x$epochs$targets$t0
  if (!is.null(x$voxels))
    tabs$changed_voxels <- x$voxels$voxels[status %in% c("lost", "gained")]
  paths <- character(0)
  for (nm in names(tabs)) {
    f <- file.path(dir, sprintf("%s_%s.csv", prefix, nm))
    if (file.exists(f) && !overwrite) {
      warning("Not overwriting existing file: ", f, call. = FALSE); next
    }
    data.table::fwrite(tabs[[nm]], f)
    paths <- c(paths, f)
  }
  invisible(paths)
}
