#' Identify the management operation carried out
#'
#' Reads the change detected between epochs and names the operation, with
#' transparent rules and every metric behind the decision:
#'
#' * **thinning**: at least `min_removed_frac` of the trees were removed;
#'   **clearcut** above `clearcut_frac`.
#' * **pruning**: at tree level, at least `min_pruned_frac` of the standing
#'   trees were classified as pruned; at area level, the lower crown lost
#'   clearly more volume than the upper crown (`pruning_index`).
#' * **thinning type**: `"systematic"` when most removed trees form
#'   complete planting rows (row removal); otherwise from the ratio between
#'   the quadratic mean diameter of removed trees and of all trees:
#'   `"from below"` (< `below_ratio`), `"from above"` (> `above_ratio`) or
#'   `"selective"`.
#'
#' The `pruning_index` is the relative volume loss in the lower crown
#' (1 m to 40% of top height) minus that in the upper crown (60% to 90%).
#' Thinning reduces all strata alike (index near 0); pruning empties the
#' lower crown only (high index).
#'
#' @param x A [silvdelta()] result.
#' @param min_removed_frac,clearcut_frac Fractions of removed trees.
#' @param min_pruned_frac Fraction of standing trees pruned (tree level).
#' @param min_pruning_index Pruning index threshold (area level).
#' @param below_ratio,above_ratio Thinning-ratio limits.
#' @param systematic_frac Minimum fraction of removed trees in complete
#'   rows to call the thinning systematic.
#' @return A one-row `data.table` with `operation`, `thinning_type`, and
#'   the metrics: `removed_frac`, `intensity_G`, `thinning_ratio`,
#'   `rows_removed`, `systematic_frac`, `pruned_frac`, `prune_height`
#'   (median, m) and `pruning_index`.
#' @export
classify_operation <- function(x, min_removed_frac = 0.05, clearcut_frac = 0.9,
                               min_pruned_frac = 0.3, min_pruning_index = 0.3,
                               below_ratio = 0.9, above_ratio = 1.1,
                               systematic_frac = 0.8) {
  tr <- x$trees; s <- x$stand
  if (x$level == "tree") {
    removed <- tr$change == "removed"; X <- tr$X; Y <- tr$Y
    standing <- !removed
    pruned_frac <- if (any(standing)) mean(tr$change[standing] == "pruned") else NA_real_
    ph <- if (any(tr$change == "pruned")) median(tr$prune_height, na.rm = TRUE) else NA_real_
  } else {
    k <- tr$status != "new"
    removed <- tr$status[k] == "removed"; X <- tr$X[k]; Y <- tr$Y[k]
    pruned_frac <- NA_real_; ph <- NA_real_
  }
  pidx <- .pruning_index(x$profile)
  rf <- mean(removed)
  rows <- .row_removal(X, Y, removed)
  thinning <- sum(removed) >= 1L && rf >= min_removed_frac
  pruning <- if (!is.na(pruned_frac)) pruned_frac >= min_pruned_frac
             else !is.na(pidx) && pidx >= min_pruning_index
  op <- if (rf >= clearcut_frac) "clearcut"
        else if (thinning && pruning) "thinning + pruning"
        else if (thinning) "thinning" else if (pruning) "pruning" else "none"
  tt <- NA_character_
  if (thinning && op != "clearcut") {
    r <- s$thinning_ratio
    tt <- if (rows$frac >= systematic_frac) "systematic"
          else if (is.na(r)) "selective"
          else if (r < below_ratio) "from below"
          else if (r > above_ratio) "from above" else "selective"
  }
  data.table::data.table(operation = op, thinning_type = tt, removed_frac = rf,
                         intensity_G = s$intensity_G, thinning_ratio = s$thinning_ratio,
                         rows_removed = rows$rows, systematic_frac = rows$frac,
                         pruned_frac = pruned_frac, prune_height = ph,
                         pruning_index = pidx)
}

.pruning_index <- function(p) {
  p <- p[volume_t0 > 0]
  if (!nrow(p)) return(NA_real_)
  H <- max(p$z_low[p$volume_t0 > 0.02 * max(p$volume_t0)]) + 0.5
  lo <- p[z_low >= 1 & z_low < 0.4 * H]; hi <- p[z_low >= 0.6 * H & z_low < 0.9 * H]
  if (!nrow(lo) || !nrow(hi)) return(NA_real_)
  (1 - sum(lo$volume_t1) / sum(lo$volume_t0)) - (1 - sum(hi$volume_t1) / sum(hi$volume_t0))
}

.row_removal <- function(X, Y, removed, min_trees = 3L) {
  n <- length(X)
  if (n < 6L || !any(removed)) return(list(rows = 0L, frac = 0))
  P <- cbind(X, Y); k <- min(5L, n)
  nn <- RANN::nn2(P, P, k = k)
  ii <- nn$nn.idx[, -1, drop = FALSE]
  ang <- atan2(Y[ii] - rep(Y, k - 1), X[ii] - rep(X, k - 1)) %% pi
  h <- tabulate(floor(ang / (pi / 36)) + 1L, 36L)
  th <- (which.max(h) - 0.5) * pi / 36
  sp <- median(nn$nn.dists[, 2])
  best <- list(rows = 0L, frac = 0)
  for (a in c(th, th + pi / 2)) {
    u <- -X * sin(a) + Y * cos(a)
    o <- order(u)
    g <- integer(n); g[o] <- cumsum(c(TRUE, diff(u[o]) > 0.35 * sp))
    full <- vapply(split(removed, g), function(r) length(r) >= min_trees && all(r),
                   logical(1))
    insys <- removed & g %in% as.integer(names(full)[full])
    fr <- sum(insys) / sum(removed)
    if (fr > best$frac) best <- list(rows = sum(full), frac = fr)
  }
  best
}
