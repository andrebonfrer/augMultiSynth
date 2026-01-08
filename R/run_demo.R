#' Run a complete end-to-end demo (simulation + estimation + diagnostics)
#'
#' Simulates a staggered-adoption multi-outcome panel with latent factor structure,
#' fits multi-outcome synthetic control with unit-level treatment effects (optionally
#' with pooled intercept adjustment), and reports accuracy diagnostics.
#'
#' @param N Integer. Number of units.
#' @param T Integer. Number of periods.
#' @param M Integer. Number of outcomes.
#' @param treated_eval Integer. Number of treated units to evaluate (subset of treated).
#' @param L Integer. Number of pre-treatment periods used for weights and intercept handling.
#' @param K Integer. Maximum event time (computes k = 0..K).
#' @param max_donors Integer. Maximum donors after screening.
#' @param screen_outcome Integer. Outcome index used for donor screening.
#' @param screen_method Character. "cor" or "mse" donor screening.
#' @param lambda Numeric. Ridge penalty on weights (solver regularization).
#' @param solver Character. Currently "fw".
#' @param pooled_adjustment Logical. If TRUE, fit pooled intercept adjustment (multisynth-like).
#' @param nu Numeric. Ridge penalty for pooled intercept adjustment (per outcome).
#' @param use_shrinkage Logical. If TRUE, compute Option-B post-τ shrinkage (multisynth_shrink).
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress.
#' @param intercept Character. Intercept handling in the weight-fitting objective passed to
#'   \code{multiout_synth}. One of \code{"none"}, \code{"global"}, or \code{"outcome"}.
#' @param standardize_outcomes Logical. Passed to \code{multiout_synth}.
#' @param outcome_weights Optional numeric vector length M. Passed to \code{multiout_synth}.
#' @param time_weights Optional numeric vector length L. Passed to \code{multiout_synth}.
#' @param eps_sd Small jitter used in standardization to avoid division by near-zero SDs.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{fit}: result of \code{multiout_synth()}
#'     \item \code{sim}: simulated data object
#'     \item \code{cor}: length-M correlations for k=0 unit effects (treated subset)
#'     \item \code{rmse}: length-M RMSE for k=0 unit effects (treated subset)
#'     \item \code{avg_1K}: list from \code{eval_avg_1K()} (unit-level avg(1:K) accuracy)
#'     \item \code{shrink}: list from \code{multisynth_shrink()} (if enabled), plus accuracy vs truth
#'   }
#'
#' @export
run_demo <- function(N = 2000, T = 200, M = 3,
                     treated_eval = 300,
                     L = 30, K = 10,
                     max_donors = 800,
                     screen_outcome = 1,
                     screen_method = c("cor", "mse"),
                     lambda = 1e-3,
                     solver = c("fw"),
                     pooled_adjustment = FALSE,
                     nu = 0,
                     use_shrinkage = TRUE,
                     seed = 1,
                     verbose = FALSE,
                     intercept = c("none","global","outcome"),
                     standardize_outcomes = FALSE,
                     outcome_weights = NULL,
                     time_weights = NULL,
                     eps_sd = 1e-8
) {
  screen_method <- match.arg(screen_method)
  solver <- match.arg(solver)
  intercept <- match.arg(intercept)

  if (!is.null(seed)) set.seed(seed)

  # 1) Simulate data
  # NOTE: simulate_multi_outcome_panel() (as previously defined) does NOT take `intercept`.
  sim <- simulate_multi_outcome_panel(N = N, T = T, M = M, K = K, L = L)

  treated_units_all <- which(is.finite(sim$treat_time))
  if (length(treated_units_all) < 1) stop("Simulation produced no treated units.")
  if (treated_eval > length(treated_units_all)) treated_eval <- length(treated_units_all)

  treated_units <- sample(treated_units_all, treated_eval)

  # 2) Fit multi-outcome SC (optionally pooled adjustment)
  fit <- multiout_synth(
    Y_list = sim$Yobs,
    treat_time = sim$treat_time,
    treated_units = treated_units,
    L = L, K = K,
    max_donors = max_donors,
    screen_outcome = screen_outcome,
    screen_method = screen_method,
    lambda = lambda,
    solver = solver,
    pooled_adjustment = pooled_adjustment,
    nu = nu,
    verbose = verbose,
    standardize_outcomes = standardize_outcomes,
    outcome_weights = outcome_weights,
    time_weights = time_weights,
    intercept = intercept,
    eps_sd = eps_sd
  )

  # 3) Accuracy for k = 0
  J <- length(treated_units)
  stopifnot(all(dim(fit$tau)[1:3] == c(J, M, K + 1)))

  tau_true_k0 <- matrix(NA_real_, nrow = J, ncol = M)

  # Prefer: sim$tau_true_event (N x M x T) or sim$tau_it (N x M x (K+1))
  if (!is.null(sim$tau_it) && length(dim(sim$tau_it)) == 3 && dim(sim$tau_it)[2] == M) {
    if (dim(sim$tau_it)[3] >= (K + 1)) {
      tau_true_k0 <- sim$tau_it[treated_units, , 1, drop = FALSE]
      tau_true_k0 <- matrix(tau_true_k0, nrow = J, ncol = M)
    }
  } else if (!is.null(sim$tau_true_event) && length(dim(sim$tau_true_event)) == 3 && dim(sim$tau_true_event)[2] == M) {
    # tau_true_event is by event time k; k=0 corresponds to multiplier at event time 0
    # For constant effects, this equals tau_unit. For dynamic, it is tau_unit * shape[1].
    tau_true_k0 <- sim$tau_true_event[treated_units, , 1, drop = FALSE]
    tau_true_k0 <- matrix(tau_true_k0, nrow = J, ncol = M)
  } else if (!is.null(sim$tau_unit)) {
    tau_true_k0 <- sim$tau_unit[treated_units, , drop = FALSE]
  } else {
    stop("Simulation object does not contain tau truth (tau_it, tau_true_event, or tau_unit).")
  }

  tau_hat_k0 <- fit$tau[, , 1, drop = FALSE]
  tau_hat_k0 <- matrix(tau_hat_k0, nrow = J, ncol = M)

  cor_out <- rmse_out <- rep(NA_real_, M)
  for (m in seq_len(M)) {
    ok <- is.finite(tau_true_k0[, m]) & is.finite(tau_hat_k0[, m])
    cor_out[m]  <- if (sum(ok) >= 3) stats::cor(tau_true_k0[ok, m], tau_hat_k0[ok, m]) else NA_real_
    rmse_out[m] <- if (sum(ok) >= 1) sqrt(mean((tau_hat_k0[ok, m] - tau_true_k0[ok, m])^2)) else NA_real_
  }

  # 4) Average over event times 1..K
  avg_1K <- eval_avg_1K(fit = fit, sim = sim, K = K)

  # 5) Optional post-τ shrinkage layer (Option B)
  shrink <- NULL
  if (isTRUE(use_shrinkage)) {
    shrink <- multisynth_shrink(tau_hat = fit$tau, K = K, use_k = 1:K, nu = NULL)

    if (!is.null(avg_1K$tau_true_avg_1K)) {
      cor_shr <- rmse_shr <- rep(NA_real_, M)
      for (m in seq_len(M)) {
        ok <- is.finite(avg_1K$tau_true_avg_1K[, m]) & is.finite(shrink$theta_shrunk[, m])
        cor_shr[m]  <- if (sum(ok) >= 3) stats::cor(avg_1K$tau_true_avg_1K[ok, m], shrink$theta_shrunk[ok, m]) else NA_real_
        rmse_shr[m] <- if (sum(ok) >= 1) sqrt(mean((shrink$theta_shrunk[ok, m] - avg_1K$tau_true_avg_1K[ok, m])^2)) else NA_real_
      }
      shrink$cor_vs_true_avg_1K <- cor_shr
      shrink$rmse_vs_true_avg_1K <- rmse_shr
    }
  }

  list(
    fit = fit,
    sim = sim,
    cor = cor_out,
    rmse = rmse_out,
    avg_1K = avg_1K,
    shrink = shrink
  )
}
