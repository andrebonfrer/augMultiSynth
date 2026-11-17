aug_problem <- function(seed, P = 40, n = 25) {
  set.seed(seed)
  A <- matrix(rnorm(P * n), P, n)
  y <- as.numeric(A %*% c(rep(0.25, 4), rep(0, n - 4))) + 0.8 + rnorm(P, sd = 0.2)  # outside the hull
  w <- as.numeric(fit_weights_fw(A, y, lambda = 0.1))
  list(A = A, y = y, w = w)
}

test_that("augmented weights follow Lemma 1 and sum to one", {
  p <- aug_problem(1)
  lam <- c(0.01, 1, 100, Inf)
  G <- augMultiSynth:::.ridge_augment(p$A, p$y, p$w, lam)
  Ac <- p$A - rowMeans(p$A)
  sc <- sum(Ac^2) / nrow(Ac)
  for (k in 1:3) {
    direct <- p$w + as.numeric(t(Ac) %*% solve(tcrossprod(Ac) + lam[k] * sc * diag(nrow(Ac)),
                                               p$y - p$A %*% p$w))
    expect_equal(G[, k], direct, tolerance = 1e-10)
  }
  expect_equal(colSums(G), rep(1, length(lam)), tolerance = 1e-12)
  expect_equal(G[, 4], p$w)
  # P > n branch gives the same answer
  pt <- aug_problem(2, P = 30, n = 50)
  A2 <- cbind(pt$A, pt$A[, 1:5] + 0.01)                # n > P
  w2 <- as.numeric(fit_weights_fw(A2, pt$y, lambda = 0.1))
  G2 <- augMultiSynth:::.ridge_augment(A2, pt$y, w2, 0.5)
  Ac2 <- A2 - rowMeans(A2); sc2 <- sum(Ac2^2) / nrow(Ac2)
  d2 <- w2 + as.numeric(t(Ac2) %*% solve(tcrossprod(Ac2) + 0.5 * sc2 * diag(nrow(Ac2)), pt$y - A2 %*% w2))
  expect_equal(G2[, 1], d2, tolerance = 1e-9)
})

test_that("the weighting form equals SC plus a ridge outcome-model correction", {
  p <- aug_problem(3)
  lam <- 0.3
  g <- augMultiSynth:::.ridge_augment(p$A, p$y, p$w, lam)[, 1]
  Ac <- p$A - rowMeans(p$A); sc <- sum(Ac^2) / nrow(Ac)
  set.seed(9); tgt <- rnorm(ncol(p$A))                   # donors' post-period outcome
  beta <- solve(tcrossprod(Ac) + lam * sc * diag(nrow(Ac)), Ac %*% (tgt - mean(tgt)))
  m_hat <- function(x) mean(tgt) + sum((x - rowMeans(p$A)) * beta)
  pred_aug <- sum(p$w * tgt) + (m_hat(p$y) - sum(p$w * apply(p$A, 2, m_hat)))
  expect_equal(sum(g * tgt), pred_aug, tolerance = 1e-10)
})

test_that("multiout_synth(augment = 'ridge') is consistent with inference and diagnostics", {
  set.seed(5)
  sim <- simulate_multi_outcome_panel(N = 200, T = 50, M = 3, K = 4, L = 12, treated_frac = 0.3,
                                      treated_loading_shift = 2)
  tr <- sim$treated[sim$treat_time[sim$treated] - 12 >= 1 & sim$treat_time[sim$treated] + 4 <= 50][1:8]
  args <- list(Y_list = sim$Yobs, treat_time = sim$treat_time, treated_units = tr, L = 12, K = 4,
               max_donors = 40, lambda = 0.1, intercept = "outcome", standardize_outcomes = TRUE)
  f0 <- do.call(multiout_synth, args)
  fa <- do.call(multiout_synth, c(args, list(augment = "ridge", augment_lambda = 0.5)))
  fi <- do.call(multiout_synth, c(args, list(augment = "ridge", augment_lambda = Inf)))
  expect_equal(sapply(fa$weights, sum), rep(1, length(tr)), tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(any(unlist(fa$weights) < 0))
  expect_equal(fi$tau, f0$tau, tolerance = 1e-12)
  expect_equal(fa$augment$weights_scm, f0$weights, tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(all(fa$augment$diagnostics$pre_rmse_aug <= fa$augment$diagnostics$pre_rmse_scm + 1e-12))
  ti <- expect_silent(tau_inference(fa, sim$Yobs, sim$treat_time))
  expect_equal(nrow(ti), length(tr) * 3)
  fp <- do.call(multiout_synth, c(args, list(augment = "ridge", augment_lambda = c(0.5, Inf))))
  expect_equal(fp$augment$path[[1]][, 1], fa$weights[[1]], tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(fp$weights, f0$weights, tolerance = 1e-12)
})

test_that("tune_augment chooses a penalty from the grid", {
  set.seed(6)
  sim <- simulate_multi_outcome_panel(N = 200, T = 50, M = 2, K = 4, L = 12, treated_frac = 0.3,
                                      treated_loading_shift = 2)
  tu <- tune_augment(sim$Yobs, sim$treat_time, L = 12, K = 4, grid = c(0.1, 1, 10), h = 3,
                     tune_units = 8, max_donors = 40, lambda = 0.1, intercept = "outcome",
                     standardize_outcomes = TRUE)
  expect_true(tu$lambda %in% c(0.1, 1, 10, Inf))
  expect_gte(tu$lambda, tu$lambda_min)
  expect_equal(nrow(tu$summary), 4)
  expect_error(tune_augment(sim$Yobs, sim$treat_time, L = 4, K = 4, h = 3), "L - h")
})

test_that("treated_loading_shift = 0 leaves the simulation unchanged", {
  set.seed(9); a <- simulate_multi_outcome_panel(N = 60, T = 30, M = 2, K = 3, L = 8)
  set.seed(9); b <- simulate_multi_outcome_panel(N = 60, T = 30, M = 2, K = 3, L = 8, treated_loading_shift = 0)
  expect_identical(a$Yobs, b$Yobs)
  set.seed(9); d <- simulate_multi_outcome_panel(N = 60, T = 30, M = 2, K = 3, L = 8, treated_loading_shift = 1.5)
  expect_equal(d$lambda_i[d$treated, 1] - a$lambda_i[a$treated, 1], rep(1.5, length(d$treated)))
  expect_equal(d$Yobs[[1]][-d$treated, ], a$Yobs[[1]][-a$treated, ])
})

test_that("augmented weights pass through scmBayesPost's conditional likelihood", {
  skip_if_not_installed("scmBayesPost")
  skip_if(utils::packageVersion("scmBayesPost") < "0.4.9")
  set.seed(5)
  sim <- simulate_multi_outcome_panel(N = 200, T = 50, M = 2, K = 4, L = 12, treated_frac = 0.3,
                                      treated_loading_shift = 2)
  tr <- sim$treated[sim$treat_time[sim$treated] - 12 >= 1 & sim$treat_time[sim$treated] + 4 <= 50][1:6]
  fa <- multiout_synth(sim$Yobs, sim$treat_time, treated_units = tr, L = 12, K = 4, max_donors = 40,
                       lambda = 0.1, intercept = "outcome", standardize_outcomes = TRUE,
                       augment = "ridge", augment_lambda = 0.5)
  ids <- as.character(fa$unit_ids); N <- length(ids); TT <- 50; tt <- sim$treat_time
  wid <- rep(seq_len(TT), each = N); ttl <- rep(tt, times = TT)
  dta <- data.table::data.table(id = rep(ids, times = TT), wID = wid, y = as.numeric(sim$Yobs[[1]]),
                                tvg.dummy = as.integer(is.finite(ttl) & wid >= ttl))
  W <- scmBayesPost::build_W_from_augMultiSynth(fa, id_universe = ids)
  expect_true(any(W < 0))
  g <- scmBayesPost::prepare_data_general(
    dta = dta, W = W, y_name = "y", f.X = y ~ 1 + tvg.dummy, id_col = "id", time_col = "wID",
    tr_col = "tvg.dummy", treat_type = "binary", second_stage = "none", first_stage = "none",
    verbose = FALSE, likelihood = "conditional", L = 12, post = 1:4)
  ti <- tau_inference(fa, sim$Yobs, tt, outcomes = 1)
  ord <- match(as.character(g$cov$treated_ids), as.character(fa$treated_unit_ids))
  b <- vapply(seq_along(g$X_list), function(k)
    stats::lm.fit(g$X_list[[k]], g$y_list[[k]])$coefficients[["tvg.dummy"]], numeric(1))
  expect_equal(b, ti$tau[ord], tolerance = 1e-10)
})
