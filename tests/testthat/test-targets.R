test_that("intensity targets register strongly misaligned epochs", {
  set.seed(7)
  sim <- simulate_stand(nx = 4, ny = 3, intervention = "both", rotation = 35,
                        offset = c(4, -3, 0.8), seed = 7)
  d0 <- detect_targets(sim$t0)
  expect_equal(nrow(d0), 5)
  ep <- register_targets(sim)
  expect_false(ep$registration$ambiguous)
  expect_equal(ep$registration$n_targets, 5)
  nn <- RANN::nn2(as.matrix(sim$targets_truth[, .(X, Y, Z)]),
                  as.matrix(ep$targets$t1[, .(X, Y, Z)]), k = 1)
  expect_lt(max(nn$nn.dists), 0.03)
  g <- target_plot(ep$targets$t0)
  expect_equal(nrow(g$polygon), 4)
  expect_false(is.na(g$midpoint))
  expect_false(is.null(g$frame))
})

test_that("a symmetric arrangement (no middle target) is flagged", {
  sim <- simulate_stand(nx = 3, ny = 3, intervention = "pruning", midpoint = FALSE,
                        rotation = 150, offset = c(2, 1, 0.3), seed = 8)
  expect_warning(ep <- register_targets(sim, refine = FALSE), "symmetric")
  expect_true(ep$registration$ambiguous)
})

test_that("symmetric targets on pre-aligned epochs resolve to the closest solution", {
  sim <- simulate_stand(nx = 3, ny = 3, intervention = "pruning", midpoint = FALSE,
                        rotation = 0.5, offset = c(0.1, -0.05, 0.08), seed = 8)
  expect_warning(ep <- register_targets(sim), "least")
  nn <- RANN::nn2(as.matrix(sim$targets_truth[, .(X, Y, Z)]),
                  as.matrix(ep$targets$t1[, .(X, Y, Z)]), k = 1)
  expect_lt(max(nn$nn.dists), 0.03)
})

test_that("the operation and thinning type are identified", {
  ops <- list(
    list(int = "pruning", mode = "random", op = "pruning", tt = NA_character_),
    list(int = "thinning", mode = "systematic", op = "thinning", tt = "systematic"),
    list(int = "thinning", mode = "below", op = "thinning", tt = "from below"),
    list(int = "both", mode = "systematic", op = "thinning + pruning", tt = "systematic"))
  for (o in ops) {
    sim <- simulate_stand(nx = 5, ny = 4, intervention = o$int, thin_mode = o$mode,
                          thin_fraction = 0.3, rotation = 25, offset = c(3, 2, 0.5),
                          seed = 11)
    for (lv in c("tree", "area")) {
      res <- silvdelta(sim, registration = "targets", level = lv, verbose = FALSE)
      expect_equal(res$operation$operation, o$op, info = paste(o$int, o$mode, lv))
      expect_equal(res$operation$thinning_type, o$tt, info = paste(o$int, o$mode, lv))
      expect_equal(res$stand$n_removed, sum(sim$truth$removed))
      expect_lt(abs(res$area_ha - res$plot$area_ha), 1e-9)
    }
  }
})

test_that("results can be expressed in the plot frame", {
  sim <- simulate_stand(nx = 3, ny = 3, intervention = "pruning", rotation = 70,
                        offset = c(10, 5, 1), seed = 9)
  res <- silvdelta(sim, registration = "targets", plot_frame = TRUE, verbose = FALSE)
  expect_gte(min(res$plot$polygon$X), -0.05)
  expect_gte(min(res$plot$polygon$Y), -0.05)
  expect_true(all(res$trees$change == "pruned"))
})
