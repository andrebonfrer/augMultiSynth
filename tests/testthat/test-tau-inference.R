make_fit <- function(intercept = "outcome", pooled = TRUE, seed = 11) {
  set.seed(seed)
  sim <- simulate_multi_outcome_panel(N = 150, T = 50, M = 2, K = 5, L = 12)
  fit <- multiout_synth(sim$Yobs, sim$treat_time, L = 12, K = 5, max_donors = 40,
                        lambda = 0.1, intercept = intercept,
                        standardize_outcomes = TRUE, pooled_adjustment = pooled, nu = 1)
  list(sim = sim, fit = fit)
}

test_that("tau equals the average of fit$tau over post event times", {
  for (ic in c("none", "outcome")) {
    s <- make_fit(intercept = ic)
    ti <- expect_silent(tau_inference(s$fit, s$sim$Yobs, s$sim$treat_time))
    J <- length(s$fit$treated_units)
    ref <- c(rowMeans(s$fit$tau[, 1, 2:6]), rowMeans(s$fit$tau[, 2, 2:6]))
    expect_equal(ti$tau, ref, tolerance = 1e-10)
    expect_equal(nrow(ti), 2 * J)

    ti0 <- tau_inference(s$fit, s$sim$Yobs, s$sim$treat_time, post = 0:5, outcomes = 2)
    expect_equal(ti0$tau, rowMeans(s$fit$tau[, 2, 1:6]), tolerance = 1e-10)
  }
})

test_that("estimate and standard error equal OLS on the gap series", {
  s <- make_fit()
  gp <- scm_gap_panel(s$fit, s$sim$Yobs, s$sim$treat_time)
  ti <- tau_inference(s$fit, s$sim$Yobs, s$sim$treat_time, level = 0.9)
  for (r in c(1, 4, nrow(ti))) {
    d <- gp[gp$treated_id == ti$treated_id[r] & gp$outcome == ti$outcome[r], ]
    o <- stats::lm(gap ~ post, data = d)
    cf <- summary(o)$coefficients["post", ]
    expect_equal(ti$tau[r], unname(cf["Estimate"]), tolerance = 1e-10)
    expect_equal(ti$se[r], unname(cf["Std. Error"]), tolerance = 1e-10)
    expect_equal(ti$p_value[r], unname(cf["Pr(>|t|)"]), tolerance = 1e-8)
    expect_equal(ti$df[r], o$df.residual)
    ci <- stats::confint(o, "post", level = 0.9)
    expect_equal(c(ti$lo[r], ti$hi[r]), unname(as.numeric(ci)), tolerance = 1e-10)
  }
})

test_that("weighted two-way FE on the pseudo-panel reproduces the point estimate", {
  s <- make_fit()
  ti <- tau_inference(s$fit, s$sim$Yobs, s$sim$treat_time, outcomes = 1)
  Y <- s$sim$Yobs[[1]]; L <- 12; ev <- c(-L:-1, 1:5)
  for (jj in 1:3) {
    j <- s$fit$treated_units[jj]; w <- s$fit$weights[[jj]]; d <- s$fit$donors[[jj]][w > 0]
    units <- c(j, d); per <- s$sim$treat_time[j] + ev
    pp <- data.frame(u = factor(rep(units, times = length(per))),
                     e = factor(rep(ev, each = length(units))),
                     y = as.numeric(Y[units, per]),
                     D = as.integer(rep(units == j, times = length(per)) & rep(ev > 0, each = length(units))),
                     wt = rep(c(1, w[w > 0]), times = length(per)))
    b <- stats::coef(stats::lm(y ~ u + e + D, data = pp, weights = wt))[["D"]]
    expect_equal(b, ti$tau[jj], tolerance = 1e-8)
    # without fixed effects the pooled regression does not recover tau
    b0 <- stats::coef(stats::lm(y ~ D, data = pp, weights = wt))[["D"]]
    expect_false(isTRUE(all.equal(b0, ti$tau[jj], tolerance = 1e-6)))
  }
})

test_that("scmBayesPost's conditional likelihood uses the same gap series", {
  skip_if_not_installed("scmBayesPost")
  skip_if(utils::packageVersion("scmBayesPost") < "0.4.8")
  s <- make_fit(); m <- 2L; L <- 12; post <- 1:5
  Y <- s$sim$Yobs[[m]]; N <- nrow(Y); TT <- ncol(Y)
  tt <- s$sim$treat_time
  wid <- rep(seq_len(TT), each = N); ttl <- rep(tt, times = TT)
  dta <- data.table::data.table(id = rep(as.character(s$fit$unit_ids), times = TT), wID = wid,
                                y = as.numeric(Y), tvg.dummy = as.integer(is.finite(ttl) & wid >= ttl))
  W <- scmBayesPost::build_W_from_augMultiSynth(s$fit, id_universe = as.character(s$fit$unit_ids),
                                                self_weight = 1)
  g <- scmBayesPost::prepare_data_general(
    dta = dta, W = W, y_name = "y", f.X = y ~ 1 + tvg.dummy, id_col = "id", time_col = "wID",
    tr_col = "tvg.dummy", treat_type = "binary", second_stage = "none", first_stage = "none",
    verbose = FALSE, likelihood = "conditional", L = L, post = post)
  gp <- scm_gap_panel(s$fit, s$sim$Yobs, tt, outcomes = m)
  ord <- match(as.character(g$cov$treated_ids), as.character(s$fit$treated_unit_ids))
  expect_equal(unlist(g$y_list), unname(unlist(split(gp$gap, gp$treated_unit)[as.character(s$fit$treated_units[ord])])),
               tolerance = 1e-10)
  ti <- tau_inference(s$fit, s$sim$Yobs, tt, outcomes = m)
  b <- vapply(seq_along(g$X_list), function(k)
    stats::lm.fit(g$X_list[[k]], g$y_list[[k]])$coefficients[["tvg.dummy"]], numeric(1))
  expect_equal(b, ti$tau[ord], tolerance = 1e-10)
})

test_that("fits without stored L need L explicitly", {
  s <- make_fit()
  f <- s$fit
  f$L <- NULL
  expect_error(tau_inference(f, s$sim$Yobs, s$sim$treat_time), "pass L explicitly")
  expect_silent(tau_inference(f, s$sim$Yobs, s$sim$treat_time, L = 12))
  expect_warning(tau_inference(f, s$sim$Yobs, s$sim$treat_time, L = 11), "differ")
})
