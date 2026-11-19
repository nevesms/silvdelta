#' Detect high-intensity targets
#'
#' Retro-reflective targets return a much stronger signal than bark, foliage
#' or soil, whatever their shape. Points above an intensity threshold are
#' grouped by proximity and each compact group is a target candidate;
#' isolated bright points (specular reflections) are discarded by
#' `min_points`. Extra candidates are tolerated: the geometric matching in
#' [match_targets()] keeps only those that form the same arrangement in both
#' epochs.
#'
#' @param cloud A point cloud with an `Intensity` column (see [as_cloud()]).
#' @param threshold Intensity threshold. Default: `median + k_mad * MAD` of
#'   the cloud intensity, which adapts to sensors with different intensity
#'   scales and to differences between surveys.
#' @param k_mad Multiplier of the MAD for the automatic threshold.
#' @param eps Grouping distance (m) between bright points of the same target.
#' @param min_points Minimum number of bright points in a target.
#' @param max_size Maximum extent (m) of a target.
#' @param max_targets Keep at most this many candidates (largest first).
#' @return A `data.table` with `target_id`, the target centre `X`, `Y`, `Z`
#'   (median of its bright points), `n`, `size` (m) and mean `intensity`.
#'   The threshold used is stored in the attribute `"threshold"`.
#' @export
detect_targets <- function(cloud, threshold = NULL, k_mad = 8, eps = 0.1,
                           min_points = 10L, max_size = 0.6, max_targets = 15L) {
  if (!data.table::is.data.table(cloud)) cloud <- as_cloud(cloud)
  if (!"Intensity" %in% names(cloud))
    stop("The cloud has no intensity values; target detection needs them. ",
         "See the 'intensity' argument of as_cloud().")
  I <- cloud$Intensity
  if (is.null(threshold)) threshold <- median(I) + k_mad * stats::mad(I)
  empty <- data.table::data.table(target_id = integer(), X = numeric(), Y = numeric(),
                                  Z = numeric(), n = integer(), size = numeric(),
                                  intensity = numeric())
  b <- cloud[Intensity >= threshold]
  if (nrow(b) >= min_points) {
    b[, cl := dbscan::dbscan(as.matrix(b[, .(X, Y, Z)]), eps = eps, minPts = 1)$cluster]
    tg <- b[, .(X = median(X), Y = median(Y), Z = median(Z), n = .N,
                size = max(diff(range(X)), diff(range(Y)), diff(range(Z))),
                intensity = mean(Intensity)), by = cl]
    tg <- tg[n >= min_points & size <= max_size]
    data.table::setorder(tg, -n)
    if (nrow(tg) > max_targets) tg <- tg[seq_len(max_targets)]
    empty <- tg[, .(target_id = seq_len(.N), X, Y, Z, n, size, intensity)]
  }
  attr(empty, "threshold") <- threshold
  empty
}

#' Match targets between epochs by their arrangement
#'
#' Finds the rigid transformation that maps the targets of t1 onto those of
#' t0 using only their mutual distances, which do not depend on how each
#' survey was oriented. Every pair of candidates and a third consistent
#' candidate define a trial transformation; the one that brings the largest
#' number of targets into agreement wins, and it is refined with all
#' matched targets.
#'
#' A regular arrangement can be matched in more than one way: the four
#' corners of a rectangle fit equally well after a 180 degree turn (a square,
#' after any 90 degree turn). This is why one extra target placed on one
#' side only (e.g. at the middle of one edge) is recommended: it makes the
#' arrangement asymmetric and the solution unique. When two different
#' transformations are equally good, the result is flagged `ambiguous` and
#' a warning is issued.
#'
#' @param t0,t1 Target tables (columns `X`, `Y`, `Z`), e.g. from
#'   [detect_targets()].
#' @param tol Distance tolerance (m) for two targets to agree.
#' @param max_tilt Maximum tilt (degrees) between the vertical axes of the
#'   two surveys. Surveys are levelled by gravity, so solutions that turn a
#'   cloud upside down (which would fit a planar arrangement equally well)
#'   are rejected.
#' @param prefer_closest When the arrangement is symmetric, choose among the
#'   equally good solutions the one that moves t1 the least. Safe only when
#'   the epochs are already roughly aligned (e.g. georeferenced).
#' @return A list with the 4x4 `matrix` (t1 -> t0), `pairs` (`id_t0`, `id_t1`,
#'   `residual` in m), `rmse` (m), `n_matched` and `ambiguous`.
#' @export
match_targets <- function(t0, t1, tol = 0.1, max_tilt = 15, prefer_closest = TRUE) {
  cmin <- cos(max_tilt * pi / 180)
  P0 <- as.matrix(t0[, c("X", "Y", "Z")]); P1 <- as.matrix(t1[, c("X", "Y", "Z")])
  n0 <- nrow(P0); n1 <- nrow(P1)
  if (n0 < 3L || n1 < 3L) stop("At least three targets are needed in each epoch.")
  D0 <- as.matrix(stats::dist(P0)); D1 <- as.matrix(stats::dist(P1))
  sols <- list()
  score <- function(R, tv) {
    Q <- sweep(P1 %*% t(R), 2, tv, "+")
    nn <- RANN::nn2(P0, Q, k = 1)
    d <- nn$nn.dists[, 1]; idx <- nn$nn.idx[, 1]
    ok <- which(d <= tol)
    ok <- ok[order(d[ok])]
    ok <- ok[!duplicated(idx[ok])]
    list(i0 = idx[ok], i1 = ok, d = d[ok])
  }
  pr <- which(upper.tri(D0), arr.ind = TRUE)
  for (r in seq_len(nrow(pr))) {
    a <- pr[r, 1]; b <- pr[r, 2]; dab <- D0[a, b]
    if (dab < 5 * tol) next
    cand <- which(abs(D1 - dab) < tol, arr.ind = TRUE)
    cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
    for (s in seq_len(nrow(cand))) {
      c1 <- cand[s, 1]; d1 <- cand[s, 2]
      for (e in setdiff(seq_len(n0), c(a, b))) {
        ab <- P0[b, ] - P0[a, ]; ae <- P0[e, ] - P0[a, ]
        cr <- c(ab[2] * ae[3] - ab[3] * ae[2], ab[3] * ae[1] - ab[1] * ae[3],
                ab[1] * ae[2] - ab[2] * ae[1])
        if (sqrt(sum(cr^2)) / dab < 3 * tol) next
        f <- which(abs(D1[c1, ] - D0[a, e]) < tol & abs(D1[d1, ] - D0[b, e]) < tol)
        for (ff in setdiff(f, c(c1, d1))) {
          kb <- .kabsch(P1[c(c1, d1, ff), ], P0[c(a, b, e), ])
          if (kb$R[3, 3] < cmin) next
          sc <- score(kb$R, kb$t)
          if (length(sc$i0) < 3L) next
          kb <- .kabsch(P1[sc$i1, , drop = FALSE], P0[sc$i0, , drop = FALSE])
          sc <- score(kb$R, kb$t)
          if (length(sc$i0) < 3L || kb$R[3, 3] < cmin) next
          sols[[length(sols) + 1L]] <- list(R = kb$R, t = kb$t, sc = sc,
                                            n = length(sc$i0),
                                            rmse = sqrt(mean(sc$d^2)))
        }
      }
    }
  }
  if (!length(sols))
    stop("No consistent target arrangement found in both epochs; ",
         "check the detection or increase 'tol'.")
  nn <- vapply(sols, `[[`, numeric(1), "n")
  rm <- vapply(sols, `[[`, numeric(1), "rmse")
  best <- sols[[order(-nn, rm)[1]]]
  ctr <- colMeans(P1)
  pos <- function(s) as.vector(s$R %*% ctr + s$t)
  pb <- pos(best)
  amb <- FALSE
  for (s in sols[nn == best$n]) {
    ang <- acos(max(-1, min(1, (sum(diag(crossprod(s$R, best$R))) - 1) / 2)))
    if (ang > 5 * pi / 180 || sqrt(sum((pos(s) - pb)^2)) > 3 * tol) { amb <- TRUE; break }
  }
  if (amb && prefer_closest) {
    mv <- vapply(sols[nn == best$n], function(s)
      acos(max(-1, min(1, (sum(diag(s$R)) - 1) / 2))) * 10 +
        sqrt(sum((pos(s) - ctr)^2)), numeric(1))
    cands <- sols[nn == best$n]
    best <- cands[[which.min(mv)]]
  }
  Tm <- diag(4); Tm[1:3, 1:3] <- best$R; Tm[1:3, 4] <- best$t
  list(matrix = Tm,
       pairs = data.table::data.table(id_t0 = best$sc$i0, id_t1 = best$sc$i1,
                                      residual = best$sc$d),
       rmse = best$rmse, n_matched = best$n, ambiguous = amb)
}

#' Automatic registration of the epochs on intensity targets
#'
#' Detects high-intensity targets in both epochs ([detect_targets()]),
#' matches them by their arrangement ([match_targets()]) and moves t1 onto
#' t0. The clouds do not need to be roughly aligned beforehand: the
#' matching is independent of the original position and orientation of
#' each survey. Optionally, the result is refined by [register_epochs()]
#' on stable structures, which also corrects the vertical offset from the
#' terrain.
#'
#' @param ep An [epochs()] object, not yet normalized.
#' @param targets_t0,targets_t1 Optional target tables (X, Y, Z) to use
#'   instead of automatic detection (e.g. centres picked by hand).
#' @param detect_args Named list of arguments for [detect_targets()].
#' @param tol,max_tilt,prefer_closest See [match_targets()].
#' @param refine Refine with [register_epochs()]?
#' @param refine_args Named list of arguments for [register_epochs()].
#' @return The `epochs` object with t1 transformed, `ep$targets` (the
#'   matched targets of each epoch, both in t0 coordinates) and
#'   `ep$registration` (`method`, `matrix`, `rmse`, `rmse_targets`,
#'   `n_targets`, `ambiguous`).
#' @export
register_targets <- function(ep, targets_t0 = NULL, targets_t1 = NULL,
                             detect_args = list(), tol = 0.1, max_tilt = 15,
                             prefer_closest = TRUE, refine = TRUE,
                             refine_args = list()) {
  .check_epochs(ep)
  if (isTRUE(ep$normalized))
    stop("register_targets() must be run before normalize_epochs().")
  if (is.null(targets_t0)) targets_t0 <- do.call(detect_targets, c(list(cloud = ep$t0), detect_args))
  if (is.null(targets_t1)) targets_t1 <- do.call(detect_targets, c(list(cloud = ep$t1), detect_args))
  targets_t0 <- data.table::as.data.table(targets_t0)
  targets_t1 <- data.table::as.data.table(targets_t1)
  if (!"target_id" %in% names(targets_t0)) targets_t0[, target_id := seq_len(.N)]
  if (!"target_id" %in% names(targets_t1)) targets_t1[, target_id := seq_len(.N)]
  m <- match_targets(targets_t0, targets_t1, tol = tol, max_tilt = max_tilt,
                     prefer_closest = prefer_closest)
  if (m$ambiguous)
    warning("The target arrangement is symmetric: more than one alignment fits ",
            "equally well. ",
            if (prefer_closest) "The solution that moves t1 the least was kept (valid only if the epochs were already roughly aligned). "
            else "",
            "Add a target on one side only (e.g. at the middle of one edge) to make it unique.",
            call. = FALSE)
  Tm <- m$matrix
  ep$t1 <- .transform_xyz(ep$t1, Tm)
  tt0 <- targets_t0[m$pairs$id_t0]
  tt1 <- .transform_xyz(targets_t1[m$pairs$id_t1], Tm)
  tt0[, target_id := seq_len(.N)]; tt1[, target_id := seq_len(.N)]
  ep$targets <- list(t0 = tt0, t1 = tt1, detected_t0 = targets_t0,
                     detected_t1 = targets_t1)
  ep$registration <- list(method = "targets", matrix = Tm, rmse = m$rmse,
                          rmse_targets = m$rmse, n_targets = m$n_matched,
                          ambiguous = m$ambiguous, dz = 0, iterations = 0L,
                          dof = "rigid")
  if (refine) {
    ra <- utils::modifyList(list(max_dist = 0.15), refine_args)
    ep <- do.call(register_epochs, c(list(ep = ep), ra))
  }
  ep
}

#' Plot geometry from the targets
#'
#' Interprets targets placed at the plot vertices plus one target at the
#' middle of one edge. The middle target is recognised as the one lying
#' between two others; the remaining targets on the convex hull are the
#' vertices. The edge holding the middle target defines a local plot
#' frame: origin at one end of that edge, X axis along it and Y axis
#' pointing into the plot. The frame is the same in every survey, so
#' results from different dates and plots share one orientation.
#'
#' @param targets Target table (X, Y, Z), usually `ep$targets$t0`.
#' @param tol Tolerance (m) to recognise the middle target.
#' @return A list with `polygon` (vertices in order, `data.table` X, Y),
#'   `vertices` and `midpoint` (row numbers in `targets`), `area_ha` and
#'   `frame` (`origin`, unit vectors `ex`, `ey`; `NULL` without a middle
#'   target).
#' @export
target_plot <- function(targets, tol = 0.15) {
  P <- as.matrix(targets[, c("X", "Y")]); n <- nrow(P)
  mid <- NA_integer_; edge <- NULL
  if (n >= 3L) for (m in seq_len(n)) for (a in seq_len(n)) for (b in seq_len(n)) {
    if (a >= b || m == a || m == b || !is.na(mid)) next
    dam <- sqrt(sum((P[a, ] - P[m, ])^2)); dmb <- sqrt(sum((P[m, ] - P[b, ])^2))
    dab <- sqrt(sum((P[a, ] - P[b, ])^2))
    if (abs(dam + dmb - dab) < tol && dam / dab > 0.25 && dam / dab < 0.75) {
      mid <- m; edge <- c(a, b)
    }
  }
  vtx <- setdiff(seq_len(n), mid)
  if (length(vtx) < 3L)
    return(list(polygon = NULL, vertices = vtx, midpoint = mid, area_ha = NA_real_,
                frame = NULL))
  h <- vtx[grDevices::chull(P[vtx, , drop = FALSE])]
  poly <- data.table::data.table(X = P[h, 1], Y = P[h, 2])
  frame <- NULL
  if (!is.na(mid)) {
    o <- P[edge[1], ]; ex <- P[edge[2], ] - o; ex <- ex / sqrt(sum(ex^2))
    ey <- c(-ex[2], ex[1])
    if (sum((colMeans(P[h, , drop = FALSE]) - o) * ey) < 0) {
      o <- P[edge[2], ]; ex <- -ex; ey <- -ey
    }
    frame <- list(origin = o, ex = ex, ey = ey)
  }
  list(polygon = poly, vertices = h, midpoint = mid,
       area_ha = .poly_area(poly$X, poly$Y) / 1e4, frame = frame)
}

#' Clip both epochs to a polygon
#'
#' @param ep An [epochs()] object.
#' @param polygon Polygon vertices (columns X, Y), e.g. from [target_plot()].
#' @param buffer Distance (m) by which the polygon is enlarged (negative
#'   values shrink it).
#' @return The clipped `epochs` object.
#' @export
clip_epochs <- function(ep, polygon, buffer = 0) {
  pg <- .buffer_poly(polygon$X, polygon$Y, buffer)
  for (e in c("t0", "t1"))
    ep[[e]] <- ep[[e]][.pip(X, Y, pg$X, pg$Y)]
  ep
}

#' Express both epochs in the local plot frame
#'
#' @param ep An [epochs()] object (not yet normalized).
#' @param frame A frame from [target_plot()].
#' @return The `epochs` object with X, Y in plot coordinates (origin at a
#'   plot vertex, X along the edge marked by the middle target).
#' @export
to_plot_frame <- function(ep, frame) {
  for (e in c("t0", "t1")) {
    d <- data.table::copy(ep[[e]])
    xy <- .frame_xy(d$X, d$Y, frame)
    d[, `:=`(X = xy$X, Y = xy$Y)]
    ep[[e]] <- d
  }
  if (!is.null(ep$targets)) for (e in names(ep$targets)) {
    d <- data.table::copy(ep$targets[[e]])
    xy <- .frame_xy(d$X, d$Y, frame)
    d[, `:=`(X = xy$X, Y = xy$Y)]
    ep$targets[[e]] <- d
  }
  ep$frame <- frame
  ep
}

#' Remove target and mast points from both epochs
#'
#' @param ep An [epochs()] object.
#' @param targets Target table (X, Y, Z) in the coordinates of `ep`.
#' @param radius Horizontal radius (m) removed around each target.
#' @param below Height (m) below the target centre down to which points are
#'   removed (the mast).
#' @return The `epochs` object without those points.
#' @export
remove_target_points <- function(ep, targets, radius = 0.25, below = 3) {
  if (is.null(targets) || nrow(targets) == 0L) return(ep)
  tp <- as.matrix(targets[, c("X", "Y")])
  for (e in c("t0", "t1")) {
    d <- ep[[e]]
    nn <- RANN::nn2(tp, as.matrix(d[, .(X, Y)]), k = 1)
    tz <- targets$Z[nn$nn.idx[, 1]]
    drop <- nn$nn.dists[, 1] < radius & d$Z > tz - below & d$Z < tz + radius
    ep[[e]] <- d[!drop]
  }
  ep
}

.frame_xy <- function(X, Y, f) {
  dx <- X - f$origin[1]; dy <- Y - f$origin[2]
  list(X = dx * f$ex[1] + dy * f$ex[2], Y = dx * f$ey[1] + dy * f$ey[2])
}

.poly_area <- function(x, y) abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2

.buffer_poly <- function(x, y, buffer) {
  if (buffer == 0) return(list(X = x, Y = y))
  cx <- mean(x); cy <- mean(y); r <- sqrt((x - cx)^2 + (y - cy)^2)
  list(X = cx + (x - cx) * (1 + buffer / r), Y = cy + (y - cy) * (1 + buffer / r))
}

.pip <- function(x, y, px, py) {
  n <- length(px); inside <- rep(FALSE, length(x)); j <- n
  for (i in seq_len(n)) {
    if (py[i] != py[j]) {
      cond <- ((py[i] > y) != (py[j] > y)) &
        (x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i]) + px[i])
      inside <- xor(inside, cond)
    }
    j <- i
  }
  inside
}
