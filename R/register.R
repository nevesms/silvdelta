#' Fine co-registration of t1 onto t0 using stable structures
#'
#' Management interventions change the crown, so a registration on the whole
#' cloud is biased by the change itself. This function aligns the epochs on
#' structures expected to be stable: (1) a trimmed ICP on the lower stem
#' layer fixes the horizontal position, rotation and tilt; (2) the vertical
#' offset is then set from the difference between the two terrain surfaces.
#' The clouds must already be roughly aligned (offsets smaller than
#' `max_dist`), e.g. by georeferencing, a manual alignment or
#' [register_targets()]. When the epochs were already registered, the new
#' transformation is composed with the previous one.
#'
#' @param ep An [epochs()] object, not yet normalized.
#' @param dof `"rigid"` (rotation + translation) or `"translation"`.
#' @param layer Height range above ground (m) of the layer used by the ICP.
#' @param thin Voxel size (m) used to thin the layer before matching.
#' @param max_dist Maximum pairing distance (m).
#' @param trim Fraction of the closest pairs kept at each iteration
#'   (trimmed ICP; discards pairs falling on removed branches or trees).
#' @param max_iter,tol Iteration limit and convergence tolerance on the RMSE.
#' @param fix_z Adjust the vertical offset from the terrain surfaces.
#' @param ground_res Cell size (m) of the terrain models used internally.
#' @param max_points Maximum number of points per epoch used in the ICP.
#' @return The `epochs` object with t1 transformed and `ep$registration`
#'   holding the 4x4 `matrix`, the residual `rmse` (m) of the stable layer,
#'   the vertical correction `dz` and the number of `iterations`.
#' @export
register_epochs <- function(ep, dof = c("rigid", "translation"), layer = c(0.5, 3),
                            thin = 0.05, max_dist = 0.5, trim = 0.8,
                            max_iter = 60, tol = 1e-6, fix_z = TRUE,
                            ground_res = 0.5, max_points = 2e5) {
  .check_epochs(ep)
  dof <- match.arg(dof)
  if (isTRUE(ep$normalized))
    stop("register_epochs() must be run before normalize_epochs().")
  g0 <- ground_model(ep$t0, res = ground_res)
  g1 <- ground_model(ep$t1, res = ground_res)
  pick <- function(d, g) {
    hh <- d$Z - ground_height(g, d$X, d$Y)
    m <- .thin_xyz(d[hh >= layer[1] & hh <= layer[2]], thin)
    if (nrow(m) > max_points) m <- m[sample.int(nrow(m), max_points), , drop = FALSE]
    m
  }
  Fx <- pick(ep$t0, g0); Mv <- pick(ep$t1, g1)
  if (nrow(Fx) < 50L || nrow(Mv) < 50L)
    stop("Too few points in the registration layer; check 'layer'.")
  Tm <- diag(4); prev <- Inf; it <- 0L
  for (it in seq_len(max_iter)) {
    nn <- RANN::nn2(Fx, Mv, k = 1)
    dd <- nn$nn.dists[, 1]
    keep <- dd <= min(max_dist, quantile(dd, trim, names = FALSE))
    if (sum(keep) < 20L) stop("Registration failed: epochs are too far apart.")
    P <- Mv[keep, , drop = FALSE]; Q <- Fx[nn$nn.idx[keep, 1], , drop = FALSE]
    if (dof == "translation") {
      R <- diag(3); tv <- colMeans(Q - P)
    } else {
      kb <- .kabsch(P, Q); R <- kb$R; tv <- kb$t
    }
    Mv <- sweep(Mv %*% t(R), 2, tv, "+")
    step <- diag(4); step[1:3, 1:3] <- R; step[1:3, 4] <- tv
    Tm <- step %*% Tm
    rmse <- sqrt(mean(dd[keep]^2))
    if (abs(prev - rmse) < tol) break
    prev <- rmse
  }
  t1 <- data.table::copy(ep$t1)
  xyz <- sweep(as.matrix(t1[, .(X, Y, Z)]) %*% t(Tm[1:3, 1:3]), 2, Tm[1:3, 4], "+")
  dz <- 0
  if (fix_z) {
    g1b <- ground_model(data.table::data.table(X = xyz[, 1], Y = xyz[, 2], Z = xyz[, 3]),
                        res = ground_res)
    gx <- g0$ox + (seq_len(nrow(g0$z)) - 0.5) * g0$res
    gy <- g0$oy + (seq_len(ncol(g0$z)) - 0.5) * g0$res
    gg <- expand.grid(X = gx, Y = gy)
    inb <- gg$X >= g1b$ox & gg$X <= g1b$ox + nrow(g1b$z) * g1b$res &
      gg$Y >= g1b$oy & gg$Y <= g1b$oy + ncol(g1b$z) * g1b$res
    if (any(inb)) {
      dz <- median(ground_height(g0, gg$X[inb], gg$Y[inb]) -
                     ground_height(g1b, gg$X[inb], gg$Y[inb]))
      xyz[, 3] <- xyz[, 3] + dz
      Tm[3, 4] <- Tm[3, 4] + dz
      Mv[, 3] <- Mv[, 3] + dz
    }
  }
  t1[, `:=`(X = xyz[, 1], Y = xyz[, 2], Z = xyz[, 3])]
  dd <- RANN::nn2(Fx, Mv, k = 1)$nn.dists[, 1]
  keep <- dd <= min(max_dist, quantile(dd, trim, names = FALSE))
  ep$t1 <- t1
  reg <- list(method = "icp", matrix = Tm, rmse = sqrt(mean(dd[keep]^2)),
              dz = dz, iterations = it, dof = dof)
  prev <- ep$registration
  if (!is.null(prev)) {
    reg$matrix <- Tm %*% prev$matrix
    reg$method <- paste(prev$method, "+ icp")
    reg$rmse_targets <- prev$rmse_targets
    reg$n_targets <- prev$n_targets
    reg$ambiguous <- prev$ambiguous
  }
  if (!is.null(ep$targets$t1)) ep$targets$t1 <- .transform_xyz(ep$targets$t1, Tm)
  ep$registration <- reg
  ep
}

.thin_xyz <- function(d, res) {
  if (nrow(d) == 0L) return(matrix(numeric(0), 0, 3))
  key <- data.table::data.table(i = floor(d$X / res), j = floor(d$Y / res),
                                k = floor(d$Z / res))
  as.matrix(d[!duplicated(key), .(X, Y, Z)])
}

.kabsch <- function(P, Q) {
  cp <- colMeans(P); cq <- colMeans(Q)
  H <- crossprod(sweep(P, 2, cp), sweep(Q, 2, cq))
  s <- svd(H)
  dsign <- sign(det(s$v %*% t(s$u)))
  R <- s$v %*% diag(c(1, 1, dsign)) %*% t(s$u)
  list(R = R, t = as.vector(cq - R %*% cp))
}

.transform_xyz <- function(d, Tm) {
  d <- data.table::copy(d)
  xyz <- sweep(as.matrix(d[, .(X, Y, Z)]) %*% t(Tm[1:3, 1:3]), 2, Tm[1:3, 4], "+")
  d[, `:=`(X = xyz[, 1], Y = xyz[, 2], Z = xyz[, 3])]
  d[]
}
