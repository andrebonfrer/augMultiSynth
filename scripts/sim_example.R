# sim_example.R
#
# Standalone simulation script for augMultiSynth
#
# This script demonstrates a full simulation–estimation–evaluation pipeline
# using multi-outcome augmented synthetic control. It is intended for
# reproducible experiments and benchmarking, not as part of the package API.
#
# To run:
#   source("scripts/sim_example.R")
#   or
#   Rscript scripts/sim_example.R

library(augMultiSynth)

set.seed(1)

res <- run_demo(
  N = 800,
  T = 120,
  M = 3,
  treated_eval = 150,
  L = 25,
  K = 8,
  max_donors = 400,
  intercept = "outcome",
  standardize_outcomes = TRUE,
  pooled_adjustment = TRUE,
  nu = 1
)

# point-at-event-time metrics (k=0)
res$cor
res$rmse

# avg over event times 1..K
res$avg_1K$cor_avg_1K
res$avg_1K$rmse_avg_1K

# Estimated vs True unit-level avg treatment effects (k=1..K), by outcome
tau_hat <- res$avg_1K$tau_hat_avg_1K
tau_true <- res$avg_1K$tau_true_avg_1K

stopifnot(all(dim(tau_hat) == dim(tau_true)))
M <- ncol(tau_hat)

op <- par(mfrow = c(1, M), mar = c(4,4,3,1))
for (m in 1:M) {
  x <- tau_true[, m]
  y <- tau_hat[, m]

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]

  plot(x, y,
       xlab = "True avg tau (k=1..K)",
       ylab = "Estimated avg tau (k=1..K)",
       main = paste0("Outcome ", m))
  abline(0, 1, lty = 2)

  # annotate with correlation and RMSE
  r <- if (length(x) >= 3) cor(x, y) else NA_real_
  rmse <- sqrt(mean((y - x)^2))
  mtext(sprintf("cor=%.3f  rmse=%.3f", r, rmse), side = 3, line = 0.2, cex = 0.85)
}
par(op)

