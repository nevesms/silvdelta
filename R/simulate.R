#' Simulate a small plantation scanned before and after an intervention
#'
#' Generates synthetic conifer-like trees (tapered stem, whorls of branches
#' with foliage) on sloping ground, and a second epoch after pruning and/or
#' thinning. The second epoch is an independent sample of the scene with
#' its own noise and a rigid misalignment, as in real repeated surveys.
#' Points carry an `Intensity` value; optional retro-reflective targets
#' (very high intensity) are placed on masts at the plot vertices and at
#' the middle of one edge. Intended for examples, tests and for learning
#' the workflow before using real data.
#'
#' @param nx,ny Number of trees per row and number of rows.
#' @param spacing Spacing (m) between trees within a row (X) and between
#'   rows (Y).
#' @param intervention `"pruning"`, `"thinning"` or `"both"`.
#' @param thin_mode `"random"`, `"below"` (smallest trees) or
#'   `"systematic"` (every third row removed).
#' @param prune_height Mean pruning height (m); varies +-0.3 m among trees.
#' @param thin_fraction Fraction of trees removed (random/below thinning).
#' @param height,dbh Mean tree height (m) and DBH (m).
#' @param slope Ground slope along X (m/m).
#' @param noise Standard deviation (m) of the ranging noise.
#' @param offset Translation (m) applied to t1 (X, Y, Z).
#' @param rotation Rotation (degrees) about the vertical axis applied to t1.
#' @param targets Place intensity targets at the plot vertices?
#' @param midpoint Also place the asymmetric target at the middle of one
#'   edge?
#' @param ground_density Ground points per m2.
#' @param seed Random seed.
#' @return An [epochs()] object with extra elements: `truth` (per-tree
#'   table with position, `row`, `dbh`, `height`, `pruned`, `prune_h`,
#'   `removed` and `removed_length`, the length of pruned branches in m),
#'   `seeds` (stem positions with 10 cm of error) and `targets_truth`
#'   (target centres in t0 coordinates).
#' @examples
#' sim <- simulate_stand(nx = 2, ny = 2, seed = 1)
#' sim$truth
#' @export
simulate_stand <- function(nx = 4, ny = 3, spacing = c(3, 2),
                           intervention = c("pruning", "thinning", "both"),
                           thin_mode = c("random", "below", "systematic"),
                           prune_height = 2.5, thin_fraction = 1 / 3,
                           height = 9, dbh = 0.18, slope = 0.05, noise = 0.01,
                           offset = c(0.04, -0.03, 0.06), rotation = 0.3,
                           targets = TRUE, midpoint = TRUE,
                           ground_density = 300, seed = NULL) {
  intervention <- match.arg(intervention)
  thin_mode <- match.arg(thin_mode)
  if (!is.null(seed)) set.seed(seed)
  n <- nx * ny
  gz <- function(x, y) slope * x + 0.02 * y + 0.05 * sin(x / 3)
  tr <- data.table::data.table(
    tree_id = seq_len(n), row = rep(seq_len(ny), each = nx),
    X = rep(seq_len(nx) * spacing[1], ny) + rnorm(n, 0, 0.15),
    Y = rep(seq_len(ny) * spacing[2], each = nx) + rnorm(n, 0, 0.15),
    height = rnorm(n, height, height * 0.06),
    dbh = pmax(0.08, rnorm(n, dbh, dbh * 0.12)))
  tr[, base := gz(X, Y)]
  tr[, removed := FALSE]
  if (intervention != "pruning") {
    nrem <- max(1L, round(thin_fraction * n))
    idx <- switch(thin_mode,
                  random = sample.int(n, nrem),
                  below = order(tr$dbh)[seq_len(nrem)],
                  systematic = which(tr$row %% 3L == 0L))
    if (!length(idx)) stop("Systematic thinning needs at least three rows (ny >= 3).")
    tr[idx, removed := TRUE]
  }
  tr[, pruned := intervention != "thinning" & !removed]
  tr[, prune_h := ifelse(pruned, prune_height + runif(n, -0.3, 0.3), NA_real_)]

  M <- data.table::rbindlist(lapply(seq_len(n), function(t) .sim_tree(tr[t])))
  M[, `:=`(prune_h = tr$prune_h[tree_id], removed = tr$removed[tree_id])]
  rl <- unique(M[part == "branch" & !is.na(prune_h) & zw < prune_h,
                 .(tree_id, branch, len)])[, .(removed_length = sum(len)), by = tree_id]
  tr[, removed_length := 0]
  tr[rl$tree_id, removed_length := rl$removed_length]

  xr <- c(spacing[1] - 1.5, nx * spacing[1] + 1.5)
  yr <- c(spacing[2] - 1.5, ny * spacing[2] + 1.5)
  ng <- round(ground_density * (diff(xr) + 1) * (diff(yr) + 1))
  G <- data.table::data.table(X = runif(ng, xr[1] - 0.5, xr[2] + 0.5),
                              Y = runif(ng, yr[1] - 0.5, yr[2] + 0.5))
  G[, `:=`(Z = gz(X, Y), Intensity = rnorm(.N, 55, 10))]

  TT <- NULL; Tpts <- NULL
  if (targets) {
    TT <- data.table::data.table(X = c(xr[1], xr[2], xr[2], xr[1]),
                                 Y = c(yr[1], yr[1], yr[2], yr[2]))
    if (midpoint) TT <- rbind(TT, data.table::data.table(X = mean(xr), Y = yr[1]))
    TT[, Z := gz(X, Y) + 1.35]
    Tpts <- data.table::rbindlist(lapply(seq_len(nrow(TT)), function(q) {
      th <- runif(250, 0, 2 * pi); zh <- runif(250, -0.05, 0.05)
      head <- data.table::data.table(X = TT$X[q] + 0.05 * cos(th),
                                     Y = TT$Y[q] + 0.05 * sin(th),
                                     Z = TT$Z[q] + zh, Intensity = rnorm(250, 240, 5))
      tm <- runif(150, 0, 2 * pi); zm <- runif(150, TT$Z[q] - 1.35, TT$Z[q] - 0.05)
      mast <- data.table::data.table(X = TT$X[q] + 0.015 * cos(tm),
                                     Y = TT$Y[q] + 0.015 * sin(tm),
                                     Z = zm, Intensity = rnorm(150, 60, 8))
      rbind(head, mast)
    }))
  }

  keep1 <- !(M$part == "branch" & !is.na(M$prune_h) & M$zw < M$prune_h) &
    !(M$removed & M$h > 0.15)
  sample_scene <- function(tree_pts) {
    stray <- data.table::data.table(X = runif(30, xr[1], xr[2]), Y = runif(30, yr[1], yr[2]),
                                    Z = runif(30, 0, height), Intensity = 235)
    d <- rbind(tree_pts[, .(X, Y, Z, Intensity)], G, Tpts, stray)
    d <- d[sample.int(nrow(d), round(0.9 * nrow(d)))]
    d[, `:=`(X = X + rnorm(.N, 0, noise), Y = Y + rnorm(.N, 0, noise),
             Z = Z + rnorm(.N, 0, noise))]
    d
  }
  t0 <- sample_scene(M); t1 <- sample_scene(M[keep1])
  a <- rotation * pi / 180
  cx <- mean(t1$X); cy <- mean(t1$Y)
  x1 <- cx + (t1$X - cx) * cos(a) - (t1$Y - cy) * sin(a) + offset[1]
  y1 <- cy + (t1$X - cx) * sin(a) + (t1$Y - cy) * cos(a) + offset[2]
  t1[, `:=`(X = x1, Y = y1, Z = Z + offset[3])]

  ep <- epochs(t0, t1)
  ep$truth <- tr[, .(tree_id, row, X, Y, dbh, height, pruned, prune_h, removed,
                     removed_length)]
  ep$seeds <- tr[, .(tree_id, X = X + rnorm(n, 0, 0.1), Y = Y + rnorm(n, 0, 0.1))]
  ep$targets_truth <- TT
  ep
}

.sim_tree <- function(t) {
  hh <- t$height; r0 <- t$dbh / 2
  rstem <- function(z) r0 * pmax((hh - z) / (hh - 1.3), 0)^0.8
  ns <- 3000L
  z <- runif(ns, 0, hh); th <- runif(ns, 0, 2 * pi)
  stem <- data.table::data.table(
    X = t$X + rstem(z) * cos(th), Y = t$Y + rstem(z) * sin(th), h = z,
    part = "stem", branch = 0L, zw = z, len = 0, Intensity = rnorm(ns, 70, 15))
  wh <- seq(0.8, hh - 0.4, by = 0.4)
  B <- data.table::data.table(zw = rep(wh, each = 4L))
  B[, `:=`(branch = seq_len(.N),
           az = rep(runif(length(wh), 0, 2 * pi), each = 4L) +
             rep((0:3) * pi / 2, length(wh)) + rnorm(.N, 0, 0.2),
           L = 0.25 + 1.3 * (1 - zw / hh)^0.7,
           inc = runif(.N, 15, 25) * pi / 180)]
  wood <- B[rep(seq_len(nrow(B)), each = 60L)]
  wood[, `:=`(tt = runif(.N), jit = 0.012)]
  leaf <- B[rep(seq_len(nrow(B)), each = 180L)]
  leaf[, `:=`(tt = runif(.N, 0.35, 1), jit = 0.05)]
  P <- rbind(wood, leaf)
  P[, rr := rstem(zw) + tt * L * cos(inc)]
  br <- P[, .(X = t$X + rr * cos(az) + rnorm(.N, 0, jit),
              Y = t$Y + rr * sin(az) + rnorm(.N, 0, jit),
              h = zw + tt * L * sin(inc) + rnorm(.N, 0, jit),
              part = "branch", branch, zw, len = L,
              Intensity = rnorm(.N, 45, 12))]
  out <- rbind(stem, br)
  out[, `:=`(Z = t$base + h, tree_id = t$tree_id)]
  out
}
