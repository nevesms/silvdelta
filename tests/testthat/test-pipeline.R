test_that("registration recovers the simulated misalignment", {
  set.seed(10)
  sim <- simulate_stand(nx = 3, ny = 2, seed = 10)
  ep <- register_epochs(sim)
  expect_lt(abs(ep$registration$matrix[3, 4] + 0.06), 0.02)
  s0 <- detect_stems(normalize_epochs(ep)$t0)
  s1 <- detect_stems(normalize_epochs(ep)$t1)
  m <- match_trees(s0, s1, max_dist = 0.1)
  expect_true(all(m$status == "kept"))
  expect_lt(mean(m$dist), 0.02)
})

test_that("stems and DBH are found without seeds", {
  set.seed(2)
  sim <- simulate_stand(nx = 3, ny = 2, seed = 2)
  ep <- normalize_epochs(register_epochs(sim))
  st <- detect_stems(ep$t0)
  expect_equal(nrow(st), 6)
  nn <- RANN::nn2(as.matrix(sim$truth[, .(X, Y)]), as.matrix(st[, .(X, Y)]), k = 1)
  expect_lt(max(nn$nn.dists), 0.05)
  expect_lt(max(abs(st$dbh - sim$truth$dbh[nn$nn.idx[, 1]])), 0.02)
})

test_that("pruning is detected and quantified per tree (seeded)", {
  set.seed(3)
  sim <- simulate_stand(nx = 4, ny = 3, intervention = "pruning", seed = 3)
  res <- silvdelta(sim, seeds = sim$seeds, level = "tree", verbose = FALSE)
  tr <- res$trees
  expect_equal(tr$tree_id, sim$truth$tree_id)
  expect_true(all(tr$change == "pruned"))
  expect_lt(mean(abs(tr$prune_height - sim$truth$prune_h)), 0.6)
  expect_gt(cor(tr$lost_volume, sim$truth$removed_length), 0.7)
  expect_equal(res$stand$n_removed, 0)
})

test_that("thinning is detected at tree and area level", {
  set.seed(4)
  sim <- simulate_stand(nx = 4, ny = 3, intervention = "thinning", seed = 4)
  res <- silvdelta(sim, seeds = sim$seeds, level = "tree", verbose = FALSE)
  expect_equal(res$trees$change == "removed", sim$truth$removed)
  expect_true(all(res$trees$change[!sim$truth$removed] == "unchanged"))

  area <- silvdelta(sim, level = "area", verbose = FALSE)
  expect_equal(area$stand$n_t0, 12)
  expect_equal(area$stand$n_removed, sum(sim$truth$removed))
  expect_equal(area$stand$intensity_N, 100 * mean(sim$truth$removed))
  expect_lt(sum(area$profile$d_volume), 0)
})

test_that("mixed intervention separates pruned and removed trees", {
  set.seed(5)
  sim <- simulate_stand(nx = 4, ny = 3, intervention = "both", seed = 5)
  res <- silvdelta(sim, level = "tree", verbose = FALSE)
  nn <- RANN::nn2(as.matrix(sim$truth[, .(X, Y)]),
                  as.matrix(res$trees[, .(X, Y)]), k = 1)$nn.idx[, 1]
  truth <- ifelse(sim$truth$removed, "removed", "pruned")[nn]
  expect_equal(res$trees$change, truth)
})

test_that("biomass calibration reports leave-one-out statistics", {
  set.seed(6)
  v <- runif(20, 0.05, 0.4)
  cal <- calibrate_biomass(v, 60 * v + rnorm(20, 0, 1))
  expect_s3_class(cal, "biomass_calibration")
  expect_lt(abs(coef(cal$model)[[1]] - 60), 5)
  expect_equal(length(predict(cal, c(0.1, 0.2))), 2)
})
