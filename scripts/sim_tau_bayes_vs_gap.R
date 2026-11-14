# scripts/sim_tau_bayes_vs_gap.R
#
# Unit-level treatment effects from the augMultiSynth simulation, three ways:
#   tau_gap   : augMultiSynth gap estimator, mean(fit$tau[j, m, 2:(K+1)])
#   tau_inf   : tau_inference() - same estimate with SE / CI from the gap-series
#               likelihood (conditional on donor paths and weights)
#   tau_bayes : scmBayesPost posterior mean under the same likelihood,
#               prepare_data_general(likelihood = "conditional"), plus the
#               sampler's hierarchical shrinkage (needs scmBayesPost >= 0.4.8)
#
# Optional: scmBayesPost's default pseudo-panel likelihood (treated + donors
# over the full panel, y ~ tvg.dummy, no fixed effects). It is slow (about 5 min
# per outcome on one core) and biased; off by default.
#
# Run from the package root:
#   source("scripts/sim_tau_bayes_vs_gap.R")

library(augMultiSynth)
library(scmBayesPost)
library(data.table)

set.seed(1)

# ---- Settings ---------------------------------------------------------------
parallel_backend <- "fork"   # "none" | "fork" | "psock" (use "psock" in RStudio)
n_cores <- max(1L, parallel::detectCores() - 1L)
n_iter  <- 1000
burn_in <- 500
run_pseudo_panel <- FALSE
out_dir <- "scripts"
L <- 30; K <- 10

# ---- 1. Simulate and fit (as in sim_example.R) -------------------------------
res <- run_demo(
  N = 1000, T = 150, M = 8, L = L, K = K,
  max_donors = 500, treated_eval = 300,
  intercept = "outcome", standardize_outcomes = TRUE,
  pooled_adjustment = TRUE, nu = 1,
  parallel_backend = parallel_backend, n_cores = n_cores,
  tune_lambda = TRUE, tune_nu = TRUE, test_ids = TRUE, verbose = TRUE
)
fit <- res$fit; sim <- res$sim
M <- length(sim$Yobs)
tr_ids <- as.character(fit$treated_unit_ids)

# ---- 2. Frequentist inference on the gap-series likelihood -------------------
ti <- as.data.table(tau_inference(fit, sim$Yobs, sim$treat_time))
ti[, treated_id := as.character(treated_id)]

tau_gap <- res$avg_1K$tau_hat_avg_1K
tau_true <- res$avg_1K$tau_true_avg_1K
cmp <- rbindlist(lapply(seq_len(M), function(m)
  data.table(treated_id = tr_ids, outcome = m,
             tau_gap = tau_gap[, m], tau_true = tau_true[, m])))
cmp <- merge(cmp, ti[, .(treated_id, outcome, tau_inf = tau, se_inf = se,
                         lo_inf = lo, hi_inf = hi)],
             by = c("treated_id", "outcome"))
cat(sprintf("max |tau_inference - gap| = %.2e\n", cmp[, max(abs(tau_inf - tau_gap))]))

# ---- 3. scmBayesPost under the conditional (gap-series) likelihood -----------
bayes_tau <- function(gdata, m, spec) {
  gdata$cov$intX <- match("tvg.dummy", gdata$cov$Xcols)
  set.seed(100 + m)
  invisible(utils::capture.output(
    post <- gibbs_postscm(gdata, n_iter = n_iter, burn_in = burn_in)
  ))
  Kx  <- length(gdata$cov$Xcols)
  idx <- (0:(gdata$cov$J0 - 1)) * Kx + gdata$cov$intX
  draws <- post$beta_samples[, idx, drop = FALSE]
  data.table(treated_id = as.character(gdata$cov$treated_ids), outcome = m, spec = spec,
             mean = colMeans(draws),
             lo = apply(draws, 2, stats::quantile, 0.025),
             hi = apply(draws, 2, stats::quantile, 0.975))
}

W_all <- build_W_from_augMultiSynth(fit, id_universe = as.character(res$unit_ids), self_weight = 1)

panel_m <- function(m) {
  Y <- sim$Yobs[[m]]; N <- nrow(Y); TT <- ncol(Y); tt <- sim$treat_time
  wid <- rep(seq_len(TT), each = N); ttl <- rep(tt, times = TT)
  data.table(id = rep(as.character(res$unit_ids), times = TT), wID = wid, y = as.numeric(Y),
             tvg.dummy = as.integer(is.finite(ttl) & wid >= ttl))
}

fit_bayes <- function(m, likelihood) {
  gdata <- prepare_data_general(
    dta = panel_m(m), W = W_all, y_name = "y", f.X = y ~ 1 + tvg.dummy,
    id_col = "id", time_col = "wID", tr_col = "tvg.dummy",
    treat_type = "binary", second_stage = "none", first_stage = "none", verbose = FALSE,
    likelihood = likelihood, L = L, post = 1:K)
  bayes_tau(gdata, m, if (likelihood == "conditional") "gap" else "pseudo_panel")
}

bayes <- rbindlist(lapply(seq_len(M), function(m) {
  message(sprintf("scmBayesPost: outcome %d of %d", m, M))
  rbind(fit_bayes(m, "conditional"), if (run_pseudo_panel) fit_bayes(m, "pseudo_panel"))
}))
cmp <- merge(cmp, bayes, by = c("treated_id", "outcome"), allow.cartesian = TRUE)

summ <- cmp[, .(
  max_abs_inf_gap = max(abs(tau_inf - tau_gap)),
  cor_bayes_gap   = cor(mean, tau_gap),
  rmse_gap_true   = sqrt(mean((tau_gap - tau_true)^2)),
  rmse_bayes_true = sqrt(mean((mean - tau_true)^2)),
  cover95_inf     = mean(lo_inf <= tau_true & tau_true <= hi_inf),
  cover95_bayes   = mean(lo <= tau_true & tau_true <= hi)
), by = .(spec, outcome)]
print(summ, digits = 3)
fwrite(cmp, file.path(out_dir, "sim_tau_bayes_vs_gap.csv"))

# ---- 4. Plot -----------------------------------------------------------------
plot_cmp <- function() {
  ncolp <- min(M, 4L); nrowp <- ceiling(M / ncolp)
  op <- par(mfrow = c(nrowp, ncolp), mar = c(4, 4, 3, 1), oma = c(0, 0, 2.5, 0))
  on.exit(par(op))
  lim <- range(cmp[, c(tau_gap, lo, hi)], finite = TRUE)
  cols <- c(gap = "#1f5fa8", pseudo_panel = "#d9822b")
  for (m in seq_len(M)) {
    plot(NA, xlim = lim, ylim = lim,
         xlab = "Gap tau (augMultiSynth, avg k=1..K)",
         ylab = "scmBayesPost tau (posterior mean)", main = paste0("Outcome ", m))
    for (s in intersect(c("pseudo_panel", "gap"), unique(cmp$spec))) {
      d <- cmp[outcome == m & spec == s]
      segments(d$tau_gap, d$lo, d$tau_gap, d$hi, col = grDevices::adjustcolor(cols[s], 0.25))
      points(d$tau_gap, d$mean, pch = 16, cex = 0.55, col = grDevices::adjustcolor(cols[s], 0.8))
    }
    abline(0, 1, lty = 2)
    s <- summ[outcome == m & spec == "gap"]
    mtext(sprintf("cor %.4f | coverage: Bayes %.2f, freq %.2f",
                  s$cor_bayes_gap, s$cover95_bayes, s$cover95_inf), side = 3, line = 0.2, cex = 0.65)
    if (m == 1 && run_pseudo_panel) {
      legend("topleft", c("gap-series likelihood", "naive pseudo-panel"),
             col = cols, pch = 16, bty = "n", cex = 0.8)
    }
  }
  mtext("Unit-level tau: scmBayesPost (gap-series likelihood) vs augMultiSynth gap; bars = 95% credible intervals",
        outer = TRUE, cex = 1)
}

out_pdf <- file.path(out_dir, "sim_tau_bayes_vs_gap.pdf")
pdf(out_pdf, width = 12, height = 7); plot_cmp(); dev.off()
png(sub("\\.pdf$", ".png", out_pdf), width = 1800, height = 1050, res = 150); plot_cmp(); dev.off()
cat("\nWrote plot to: ", out_pdf, "\n", sep = "")
