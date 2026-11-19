#' Read or coerce stem seeds
#'
#' Seeds are approximate stem positions given by the user (for instance
#' points picked in CloudCompare, or a stem map from the field inventory).
#'
#' @param x A file path, `data.frame` or matrix. The first two numeric
#'   columns are taken as X and Y unless columns named `X` and `Y` exist. A
#'   column named `tree_id` (or `id`) provides the tree identifiers;
#'   otherwise trees are numbered in row order.
#' @return A `data.table` with `tree_id`, `X`, `Y`.
#' @export
as_seeds <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- data.table::fread(x)
  if (is.matrix(x)) x <- data.table::as.data.table(x)
  x <- data.table::as.data.table(x)
  nm <- toupper(names(x))
  pos <- match(c("X", "Y"), nm)
  if (anyNA(pos)) {
    num <- which(vapply(x, is.numeric, logical(1)))
    idc <- which(nm %in% c("TREE_ID", "ID"))
    num <- setdiff(num, idc)
    if (length(num) < 2L) stop("Seeds need two numeric columns (X, Y).")
    pos <- num[1:2]
  }
  idc <- which(nm %in% c("TREE_ID", "ID"))
  id <- if (length(idc)) x[[idc[1]]] else seq_len(nrow(x))
  data.table::data.table(tree_id = id, X = as.numeric(x[[pos[1]]]),
                         Y = as.numeric(x[[pos[2]]]))
}

#' Detect stems and estimate DBH, with or without seeds
#'
#' Stems are located from their vertical continuity: the cloud between
#' `column` heights is projected on a fine horizontal grid and only cells
#' occupied along most of that height range are kept, which discards
#' branches and understorey. Each group of such cells is a stem candidate
#' and a circle is fitted by RANSAC to the points of the `slice` around
#' breast height.
#'
#' Without `seeds`, every valid candidate becomes a tree (fully automatic).
#' With `seeds`, each seed is linked to the nearest candidate within
#' `seed_radius`; if there is none, a circle is searched directly around
#' the seed; if that also fails, the seed position is kept with
#' `found = FALSE` so that the tree is not lost from the analysis. Seeds
#' therefore control *which* trees are analysed and their identifiers,
#' which keeps the link with field inventory data.
#'
#' @param cloud A height-normalized cloud (`data.table` with X, Y, Z = height).
#' @param seeds Optional seeds; anything accepted by [as_seeds()].
#' @param slice Height range (m) used for the circle fit (DBH).
#' @param column Height range (m) used for the vertical continuity test.
#' @param cell Horizontal cell size (m) of the continuity grid.
#' @param min_vertical Minimum fraction of 10-cm height bins of `column`
#'   that must be occupied in a cell.
#' @param eps Distance (m) used to group continuity cells into candidates.
#' @param dbh_range Accepted DBH range (m).
#' @param seed_radius Search radius (m) around each seed.
#' @param tol Inlier tolerance (m) of the circle fit; use about the ranging
#'   noise of the sensor (0.01-0.03 m for SLAM scanners).
#' @param min_coverage Minimum fraction of the circumference (12 sectors)
#'   covered by inliers.
#' @return A `data.table` with one row per tree: `tree_id`, `X`, `Y`, `dbh`
#'   (m), fit `rmse`, `coverage`, number of inliers `n` and `found`.
#'   The circle fit uses random sampling; call `set.seed()` for exact
#'   reproducibility.
#' @export
detect_stems <- function(cloud, seeds = NULL, slice = c(1.1, 1.5),
                         column = c(0.8, 2.4), cell = 0.05, min_vertical = 0.6,
                         eps = 0.12, dbh_range = c(0.05, 0.8), seed_radius = 0.5,
                         tol = 0.02, min_coverage = 0.4) {
  if (!data.table::is.data.table(cloud)) cloud <- as_cloud(cloud)
  r_range <- dbh_range / 2
  empty <- data.table::data.table(tree_id = integer(), X = numeric(), Y = numeric(),
                                  dbh = numeric(), rmse = numeric(),
                                  coverage = numeric(), n = integer(), found = logical())
  col <- cloud[Z >= column[1] & Z <= column[2], .(X, Y, Z)]
  sl <- col[Z >= slice[1] & Z <= slice[2]]
  cand <- empty
  if (nrow(col) > 0L) {
    nbin <- ceiling(diff(column) / 0.1)
    key <- unique(data.table::data.table(ci = floor(col$X / cell), cj = floor(col$Y / cell),
                                         zb = floor((col$Z - column[1]) / 0.1)))
    occ <- key[, .(nz = .N), by = .(ci, cj)][nz >= min_vertical * nbin]
    if (nrow(occ) > 0L && nrow(sl) > 0L) {
      occ[, `:=`(cx = (ci + 0.5) * cell, cy = (cj + 0.5) * cell)]
      occ[, cl := dbscan::dbscan(cbind(cx, cy), eps = eps, minPts = 1)$cluster]
      cen <- occ[, .(X = mean(cx), Y = mean(cy)), by = cl]
      rr <- r_range[2] + 0.1
      nn <- RANN::nn2(as.matrix(cen[, .(X, Y)]), as.matrix(sl[, .(X, Y)]), k = 1)
      sl2 <- sl[, .(X, Y, id = nn$nn.idx[, 1], d = nn$nn.dists[, 1])][d <= rr]
      fits <- lapply(split(sl2, by = "id"), function(p)
        .circle_ransac(p$X, p$Y, tol = tol, r_range = r_range))
      fits <- fits[!vapply(fits, is.null, logical(1))]
      if (length(fits)) {
        cand <- data.table::rbindlist(lapply(fits, data.table::as.data.table))
        cand <- cand[coverage >= min_coverage & n >= 12L]
        if (nrow(cand)) {
          vf <- vapply(seq_len(nrow(cand)), function(q) {
            dd <- sqrt((col$X - cand$X[q])^2 + (col$Y - cand$Y[q])^2)
            zz <- col$Z[abs(dd - cand$r[q]) < 2 * tol]
            length(unique(floor((zz - column[1]) / 0.1))) / nbin
          }, numeric(1))
          hol <- vapply(seq_len(nrow(cand)), function(q) {
            dd <- sqrt((sl$X - cand$X[q])^2 + (sl$Y - cand$Y[q])^2)
            sum(dd < cand$r[q] - 3 * tol) / max(1, cand$n[q])
          }, numeric(1))
          cand <- cand[vf >= min_vertical & hol <= 0.2]
        }
        cand <- .dedup_stems(cand)
      }
    }
  }
  if (is.null(seeds)) {
    if (nrow(cand) == 0L) return(empty)
    data.table::setorder(cand, Y, X)
    return(cand[, .(tree_id = seq_len(.N), X, Y, dbh = 2 * r, rmse, coverage, n,
                    found = TRUE)])
  }
  seeds <- as_seeds(seeds)
  out <- vector("list", nrow(seeds))
  for (s in seq_len(nrow(seeds))) {
    sx <- seeds$X[s]; sy <- seeds$Y[s]; fit <- NULL
    if (nrow(cand) > 0L) {
      dd <- sqrt((cand$X - sx)^2 + (cand$Y - sy)^2)
      if (min(dd) <= seed_radius) fit <- as.list(cand[which.min(dd)])
    }
    if (is.null(fit) && nrow(sl) > 0L) {
      p <- sl[(X - sx)^2 + (Y - sy)^2 <= seed_radius^2]
      f <- .circle_ransac(p$X, p$Y, tol = tol, r_range = r_range)
      if (!is.null(f) && f$coverage >= min_coverage && f$n >= 12L &&
          sqrt((f$X - sx)^2 + (f$Y - sy)^2) <= seed_radius) fit <- f
    }
    out[[s]] <- if (is.null(fit))
      data.table::data.table(X = sx, Y = sy, dbh = NA_real_, rmse = NA_real_,
                             coverage = NA_real_, n = 0L, found = FALSE)
    else data.table::data.table(X = fit$X, Y = fit$Y, dbh = 2 * fit$r, rmse = fit$rmse,
                                coverage = fit$coverage, n = as.integer(fit$n),
                                found = TRUE)
  }
  out <- data.table::rbindlist(out)
  cbind(data.table::data.table(tree_id = seeds$tree_id), out)
}

.dedup_stems <- function(cand) {
  if (nrow(cand) < 2L) return(cand)
  data.table::setorder(cand, -n)
  keep <- rep(TRUE, nrow(cand))
  for (a in seq_len(nrow(cand))) {
    if (!keep[a]) next
    dd <- sqrt((cand$X - cand$X[a])^2 + (cand$Y - cand$Y[a])^2)
    dup <- dd < pmax(cand$r, cand$r[a]) * 1.5 & seq_len(nrow(cand)) > a
    keep[dup] <- FALSE
  }
  cand[keep]
}

.circle_kasa <- function(x, y) {
  mx <- mean(x); my <- mean(y); u <- x - mx; v <- y - my
  s <- tryCatch(qr.solve(cbind(2 * u, 2 * v, 1), u^2 + v^2), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  r2 <- s[3] + s[1]^2 + s[2]^2
  if (!is.finite(r2) || r2 <= 0) return(NULL)
  c(s[1] + mx, s[2] + my, sqrt(r2))
}

.circle_ransac <- function(x, y, tol, r_range, n_iter = 300L, max_n = 1500L) {
  n <- length(x)
  if (n < 12L) return(NULL)
  if (n > max_n) { s <- sample.int(n, max_n); x <- x[s]; y <- y[s]; n <- max_n }
  best <- NULL; bk <- 0L
  for (it in seq_len(n_iter)) {
    s <- sample.int(n, 3L)
    c3 <- .circle_kasa(x[s], y[s])
    if (is.null(c3) || c3[3] < r_range[1] || c3[3] > r_range[2]) next
    kk <- sum(abs(sqrt((x - c3[1])^2 + (y - c3[2])^2) - c3[3]) < tol)
    if (kk > bk) { bk <- kk; best <- c3 }
  }
  if (is.null(best)) return(NULL)
  for (it in 1:3) {
    inl <- abs(sqrt((x - best[1])^2 + (y - best[2])^2) - best[3]) < tol
    if (sum(inl) < 8L) break
    f <- .circle_kasa(x[inl], y[inl])
    if (is.null(f) || f[3] < r_range[1] || f[3] > r_range[2]) break
    best <- f
  }
  res <- sqrt((x - best[1])^2 + (y - best[2])^2) - best[3]
  inl <- abs(res) < tol
  if (!any(inl)) return(NULL)
  ang <- atan2(y[inl] - best[2], x[inl] - best[1])
  sect <- floor((ang + pi) / (2 * pi) * 12) %% 12
  list(X = best[1], Y = best[2], r = best[3], n = sum(inl),
       rmse = sqrt(mean(res[inl]^2)), coverage = length(unique(sect)) / 12)
}

#' Test whether each stem is still standing in a cloud
#'
#' Counts the points of a cloud inside a vertical cylinder around each stem
#' and measures how much of the height range is occupied. A standing stem
#' occupies almost all height bins; a felled tree leaves the cylinder empty.
#' This is the criterion used to flag trees removed by thinning.
#'
#' @param cloud A height-normalized cloud (usually `ep$t1`).
#' @param stems Stem table from [detect_stems()].
#' @param zrange Height range (m) tested; keep it above stump height and
#'   harvest residues.
#' @param buffer Distance (m) added to the stem radius.
#' @param bin Height bin (m).
#' @param min_occupancy Minimum fraction of occupied bins to call the stem
#'   present.
#' @return `stems` with the columns `occupancy`, `n_pts` and `present`.
#' @export
stem_presence <- function(cloud, stems, zrange = c(1, 2.5), buffer = 0.1,
                          bin = 0.1, min_occupancy = 0.5) {
  stems <- data.table::copy(stems)
  nbin <- ceiling(diff(zrange) / bin)
  stems[, `:=`(occupancy = 0, n_pts = 0L)]
  sl <- cloud[Z >= zrange[1] & Z < zrange[2], .(X, Y, Z)]
  if (nrow(sl) > 0L && nrow(stems) > 0L) {
    nn <- RANN::nn2(as.matrix(stems[, .(X, Y)]), as.matrix(sl[, .(X, Y)]), k = 1)
    rad <- ifelse(is.na(stems$dbh), 0.15, stems$dbh / 2) + buffer
    idx <- nn$nn.idx[, 1]
    keep <- nn$nn.dists[, 1] <= rad[idx]
    if (any(keep)) {
      a <- data.table::data.table(idx = idx[keep],
                                  zb = floor((sl$Z[keep] - zrange[1]) / bin))
      a <- a[, .(nb = data.table::uniqueN(zb), n = .N), by = idx]
      stems[a$idx, `:=`(occupancy = a$nb / nbin, n_pts = a$n)]
    }
  }
  stems[, present := occupancy >= min_occupancy]
  stems[]
}

#' Match stems between epochs
#'
#' One-to-one matching of independently detected stems of t0 and t1 by
#' position (closest pairs first). Stems of t0 without a partner are flagged
#' `"removed"`; stems of t1 without a partner are flagged `"new"` (usually
#' trees hidden by occlusion in t0).
#'
#' @param stems_t0,stems_t1 Stem tables from [detect_stems()].
#' @param max_dist Maximum distance (m) between matched stems.
#' @return A `data.table` with `tree_id` (t0), `id_t1`, `X`, `Y`, `dbh_t0`,
#'   `dbh_t1`, `dist` and `status` (`"kept"`, `"removed"` or `"new"`).
#' @export
match_trees <- function(stems_t0, stems_t1, max_dist = 0.3) {
  s0 <- stems_t0; s1 <- stems_t1
  out <- data.table::data.table(tree_id = s0$tree_id, id_t1 = NA_integer_,
                                X = s0$X, Y = s0$Y, dbh_t0 = s0$dbh,
                                dbh_t1 = NA_real_, dist = NA_real_, status = "removed")
  used <- rep(FALSE, nrow(s1))
  if (nrow(s0) > 0L && nrow(s1) > 0L) {
    kk <- min(3L, nrow(s1))
    nn <- RANN::nn2(as.matrix(s1[, .(X, Y)]), as.matrix(s0[, .(X, Y)]), k = kk)
    pr <- data.table::data.table(a = rep(seq_len(nrow(s0)), kk),
                                 b = as.vector(nn$nn.idx), d = as.vector(nn$nn.dists))
    pr <- pr[d <= max_dist][order(d)]
    done <- rep(FALSE, nrow(s0))
    for (r in seq_len(nrow(pr))) {
      a <- pr$a[r]; b <- pr$b[r]
      if (done[a] || used[b]) next
      done[a] <- TRUE; used[b] <- TRUE
      data.table::set(out, a, c("id_t1", "dbh_t1", "dist", "status"),
                      list(as.integer(b), s1$dbh[b], pr$d[r], "kept"))
    }
  }
  if (any(!used)) {
    nw <- s1[!used]
    out <- rbind(out, data.table::data.table(
      tree_id = NA_integer_, id_t1 = which(!used), X = nw$X, Y = nw$Y,
      dbh_t0 = NA_real_, dbh_t1 = nw$dbh, dist = NA_real_, status = "new"))
  }
  out[]
}
