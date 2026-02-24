# scripts/sim_example.R
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

# Choose parallel settings:
# - In RStudio: prefer parallel_backend="psock"
# - In macOS terminal: "fork" is fast
# - For no parallel: parallel=FALSE
parallel_flag <- FALSE
parallel_backend <- "none"   # "none" | "fork" | "psock"
parallel_backend <- "psock"

n_cores <- 8

res <- run_demo(
  N = 1000,
  T = 150,
  M = 8,
  L = 30,
  K = 10,
  max_donors = 500,
  treated_eval = 300,
  intercept = "outcome",
  standardize_outcomes = TRUE,
  pooled_adjustment = TRUE,
  nu = 1,
  parallel_backend = parallel_backend,
  n_cores = n_cores,
  tune_lambda = TRUE,
  tune_nu = TRUE,
  test_ids = TRUE,
  verbose = TRUE
)

# If it returns without error, the ID mappings are correct.
cat("Treated unit IDs:\n")
print(head(res$fit$treated_unit_ids))
cat("\nDonor IDs for first treated unit:\n")
print(head(res$fit$donor_ids[[1]]))

# point-at-event-time metrics (k=0)
cat("\nCor (k=0):\n"); print(res$cor)
cat("\nRMSE (k=0):\n"); print(res$rmse)

# avg over event times 1..K
cat("\nCor avg_1K:\n"); print(res$avg_1K$cor_avg_1K)
cat("\nRMSE avg_1K:\n"); print(res$avg_1K$rmse_avg_1K)

# Estimated vs True unit-level avg treatment effects (k=1..K), by outcome
tau_hat  <- res$avg_1K$tau_hat_avg_1K
tau_true <- res$avg_1K$tau_true_avg_1K

stopifnot(all(dim(tau_hat) == dim(tau_true)))
M <- ncol(tau_hat)

# ---- Plotting: use a device big enough (PDF recommended) ----
out_pdf <- file.path("scripts", "sim_example_scatter.pdf")
pdf(out_pdf, width = 12, height = 7)

# 2x4 layout works well for M=8; generalize via rows/cols
ncolp <- if (M <= 4) M else 4
nrowp <- ceiling(M / ncolp)

op <- par(mfrow = c(nrowp, ncolp), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (m in seq_len(M)) {
  x <- tau_true[, m]
  y <- tau_hat[, m]

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]

  plot(
    x, y,
    xlab = "True avg tau (k=1..K)",
    ylab = "Estimated avg tau (k=1..K)",
    main = paste0("Outcome ", m)
  )
  abline(0, 1, lty = 2)

  r <- if (length(x) >= 3) stats::cor(x, y) else NA_real_
  rmse <- if (length(x) >= 1) sqrt(mean((y - x)^2)) else NA_real_
  mtext(sprintf("cor=%.3f  rmse=%.3f", r, rmse), side = 3, line = 0.2, cex = 0.85)
}
mtext("Unit-level avg(1:K): Estimated vs True", outer = TRUE, cex = 1.1)
par(op)
dev.off()

cat("\nWrote plot to: ", out_pdf, "\n", sep = "")
