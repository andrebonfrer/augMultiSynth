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

res <- run_demo(
  N = 2000, T = 200, M = 3,
  treated_eval = 300,
  L = 50, K = 10,
  max_donors = 800,
  outcome_specific_intercepts = TRUE,
  standardize_outcomes = TRUE
)

res$cor
res$avg_1K$cor_avg_1K

m <- 1
plot(
  res$avg_1K$tau_true_avg_1K[, m],
  res$avg_1K$tau_hat_avg_1K[, m],
  xlab = "True avg tau (k = 1..K)",
  ylab = "Estimated avg tau (k = 1..K)"
)
abline(0, 1)
