# scripts/augment_experiment.R
#
# Ridge augmentation vs plain synthetic control when treated units sit
# increasingly far outside the donors' convex hull (treated_loading_shift on
# the first factor). For each shift and replication: fit plain SC and
# ridge-augmented SC (penalty from tune_augment(), one-SE rule and minimum),
# then compare against the true effects:
#   post_rmse : RMSE of tau_hat_k - tau_k over k = 1..K (= counterfactual error)
#   avg_rmse  : RMSE of the unit-level average effect over k = 1..K
#   coverage  : share of 95% tau_inference() intervals covering the true average
#   neg_mass  : mean negative weight mass (extrapolation)
#
# Run from the package root:  source("scripts/augment_experiment.R")

library(augMultiSynth)

parallel_backend <- "fork"           # "psock" in RStudio
n_cores <- max(1L, parallel::detectCores() - 1L)
N <- 600; TT <- 100; M <- 6; L <- 30; K <- 10
n_eval <- 60; reps <- 5
shifts <- c(0, 1, 2, 3)
out_dir <- "scripts"

base_args <- list(L = L, K = K, max_donors = 150, lambda = 0.1, intercept = "outcome",
                  standardize_outcomes = TRUE, parallel_backend = parallel_backend,
                  n_cores = n_cores)

evaluate <- function(fit, sim, label, shift, rep) {
  tt <- sim$tau_true_event[fit$treated_units, , 2:(K + 1), drop = FALSE]
  th <- fit$tau[, , 2:(K + 1), drop = FALSE]
  ti <- tau_inference(fit, sim$Yobs, sim$treat_time)
  true_avg <- as.numeric(apply(tt, c(1, 2), mean))           # J x M, column-major = ti order
  neg <- if (is.null(fit$augment)) 0 else mean(fit$augment$diagnostics$neg_mass)
  data.frame(shift = shift, rep = rep, method = label,
             lambda = if (is.null(fit$augment)) NA else fit$augment$lambda,
             post_rmse = sqrt(mean((th - tt)^2)),
             avg_rmse = sqrt(mean((ti$tau - true_avg)^2)),
             coverage = mean(ti$lo <= true_avg & true_avg <= ti$hi),
             neg_mass = neg)
}

res <- list()
for (s in shifts) for (r in seq_len(reps)) {
  set.seed(1000 * s + r)
  sim <- simulate_multi_outcome_panel(N = N, T = TT, M = M, L = L, K = K, treated_frac = 0.3,
                                      treated_loading_shift = s)
  feas <- sim$treated[sim$treat_time[sim$treated] - L >= 1 & sim$treat_time[sim$treated] + K <= TT]
  tr <- sort(sample(feas, min(n_eval, length(feas))))
  args <- c(list(Y_list = sim$Yobs, treat_time = sim$treat_time, treated_units = tr), base_args)
  f0 <- do.call(multiout_synth, args)
  tu <- do.call(tune_augment, c(list(Y_list = sim$Yobs, treat_time = sim$treat_time, L = L, K = K,
                                     treated_units = tr, tune_units = 40), base_args[-(1:2)]))
  f1 <- do.call(multiout_synth, c(args, list(augment = "ridge", augment_lambda = tu$lambda)))
  f2 <- do.call(multiout_synth, c(args, list(augment = "ridge", augment_lambda = tu$lambda_min)))
  res[[length(res) + 1L]] <- rbind(evaluate(f0, sim, "SC", s, r),
                                   evaluate(f1, sim, "ridge ASCM (one-SE)", s, r),
                                   evaluate(f2, sim, "ridge ASCM (min CV)", s, r))
  message(sprintf("shift %g rep %d done", s, r))
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_dir, "augment_experiment.csv"), row.names = FALSE)

summ <- aggregate(cbind(post_rmse, avg_rmse, coverage, neg_mass) ~ shift + method, data = res, FUN = mean)
summ <- summ[order(summ$shift, summ$method), ]
print(summ, digits = 3, row.names = FALSE)

# ---- Plot ---------------------------------------------------------------------
plot_aug <- function() {
  op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1)); on.exit(par(op))
  cols <- c("SC" = "#222222", "ridge ASCM (one-SE)" = "#1f5fa8", "ridge ASCM (min CV)" = "#d9822b")
  for (v in c("post_rmse", "coverage", "neg_mass")) {
    yl <- range(res[[v]], if (v == "coverage") c(0.95, 1) else 0)
    plot(NA, xlim = range(shifts) + c(-0.2, 0.2), ylim = yl, xlab = "Treated loading shift (factor 1)",
         ylab = c(post_rmse = "Post-adoption RMSE (tau_k)", coverage = "95% interval coverage",
                  neg_mass = "Negative weight mass")[v],
         main = c(post_rmse = "Twin error after adoption", coverage = "tau_inference() coverage",
                  neg_mass = "Extrapolation")[v])
    for (k in seq_along(cols)) {
      m <- names(cols)[k]; d <- res[res$method == m, ]
      off <- (k - 2) * 0.08
      points(d$shift + off, d[[v]], col = grDevices::adjustcolor(cols[k], 0.4), pch = 16, cex = 0.8)
      a <- tapply(d[[v]], d$shift, mean)
      lines(as.numeric(names(a)) + off, a, col = cols[k], lwd = 2, type = "b", pch = 19)
    }
    if (v == "coverage") abline(h = 0.95, lty = 3)
    if (v == "post_rmse") legend("topleft", names(cols), col = cols, lwd = 2, bty = "n", cex = 0.85)
  }
}
pdf(file.path(out_dir, "augment_experiment.pdf"), width = 13, height = 4.6); plot_aug(); dev.off()
png(file.path(out_dir, "augment_experiment.png"), width = 1950, height = 690, res = 150); plot_aug(); dev.off()
