fw_problem <- function(seed, P = 120, n = 60, lam = 0.1) {
  set.seed(seed)
  X <- matrix(rnorm(P * n), P, n)
  w_true <- rep(0, n); w_true[sample(n, 5)] <- rep(0.2, 5)
  y <- as.numeric(X %*% w_true + rnorm(P, sd = 0.5)) + 0.3   # treated slightly outside the hull
  list(X = X, y = y, lam = lam)
}
objective <- function(p, w) sum((p$y - p$X %*% w)^2) + p$lam * sum(w^2)

test_that("solution satisfies the simplex KKT conditions", {
  for (seed in 1:4) {
    p <- fw_problem(seed)
    w <- fit_weights_fw(p$X, p$y, lambda = p$lam)
    expect_true(attr(w, "fw_converged"))
    expect_lte(attr(w, "fw_rel_gap"), 1e-8)
    expect_equal(sum(w), 1, tolerance = 1e-12)
    expect_true(all(w >= 0))
    g <- as.numeric(2 * (crossprod(p$X) %*% w - crossprod(p$X, p$y)) + 2 * p$lam * w)
    act <- w > 1e-10
    scale <- max(1, abs(objective(p, w)))
    expect_lte(max(g[act]) - min(g), 1e-5 * scale)          # active gradients equal the minimum
  }
})

test_that("objective matches an exact QP solver", {
  skip_if_not_installed("quadprog")
  for (seed in 1:3) {
    p <- fw_problem(seed)
    n <- ncol(p$X)
    D <- 2 * (crossprod(p$X) + p$lam * diag(n)) + 1e-10 * diag(n)
    qp <- quadprog::solve.QP(D, 2 * crossprod(p$X, p$y), cbind(1, diag(n)), c(1, rep(0, n)), meq = 1)
    w <- fit_weights_fw(p$X, p$y, lambda = p$lam)
    expect_lte(objective(p, w) - objective(p, qp$solution), 1e-8 * objective(p, qp$solution))
    expect_lte(sum(abs(w - qp$solution)), 1e-4)
  }
})

test_that("warm start and tolerance behave", {
  p <- fw_problem(5)
  w1 <- fit_weights_fw(p$X, p$y, lambda = p$lam)
  w2 <- fit_weights_fw(p$X, p$y, lambda = p$lam, w0 = as.numeric(w1))
  expect_lte(attr(w2, "fw_iter"), 2)
  expect_equal(as.numeric(w2), as.numeric(w1), tolerance = 1e-8)
  w3 <- fit_weights_fw(p$X, p$y, lambda = p$lam, max_iter = 3)
  expect_false(attr(w3, "fw_converged"))
  expect_error(fit_weights_fw(p$X, p$y, w0 = 1:3), "w0")
})

test_that("multiout_synth reports solver convergence", {
  set.seed(3)
  sim <- simulate_multi_outcome_panel(N = 120, T = 40, M = 2, K = 3, L = 10)
  fit <- multiout_synth(sim$Yobs, sim$treat_time, L = 10, K = 3, max_donors = 30, lambda = 0.1)
  expect_true(all(fit$fw_diagnostics$converged))
  expect_equal(nrow(fit$fw_diagnostics), length(fit$treated_units))
  expect_warning(multiout_synth(sim$Yobs, sim$treat_time, L = 10, K = 3, max_donors = 30,
                                lambda = 0.1, fw_max_iter = 2), "duality-gap")
})
