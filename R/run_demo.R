#' Run a complete end-to-end demo (simulation + estimation + diagnostics)
#'
#' Simulates a staggered-adoption multi-outcome panel with latent factor structure,
#' fits multi-outcome synthetic control with unit-level treatment effects (optionally
#' with pooled intercept adjustment), and reports accuracy diagnostics.
#'
#' @param N Integer. Number of units.
#' @param T Integer. Number of periods.
#' @param M Integer. Number of outcomes.
#' @param treated_eval Integer. Number of treated units to estimate (subset of treated).
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
#' @param test_ids Logical. If TRUE, runs internal mapping checks for treated and donor IDs.
#'
#' @param tune_lambda Logical. If TRUE, runs \code{tune_lambda()} on a subset and uses the
#'   selected lambda. If FALSE, uses the supplied \code{lambda}.
#' @param lambda_grid Numeric vector of candidate lambdas (only if tune_lambda=TRUE).
#' @param holdout Integer. Holdout length for CV (only if tune_lambda=TRUE).
#' @param tune_units Integer. Number of treated units used for CV (only if tune_lambda=TRUE).
#'
#' @param parallel_backend One of "none", "psock", or "fork". Passed to \code{multiout_synth}
#'   and (if used) \code{tune_lambda}.
#' @param n_cores Integer. Number of workers/cores for the chosen backend.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{fit}: result of \code{multiout_synth()}
#'     \item \code{sim}: simulated data object
#'     \item \code{unit_ids}: unit identifiers used
#'     \item \code{lambda_used}: lambda used in estimation (tuned or supplied)
#'     \item \code{cor}: length-M correlations for k=0 unit effects (treated subset)
#'     \item \code{rmse}: length-M RMSE for k=0 unit effects (treated subset)
#'     \item \code{avg_1K}: list from \code{eval_avg_1K()} (unit-level avg(1:K) accuracy)
#'     \item \code{shrink}: list from \code{multisynth_shrink()} (if enabled), plus accuracy vs truth
#'     \item \code{tune}: NULL or tune_lambda results (if tune_lambda=TRUE)
#'   }
#'
#' @export
run_demo <- function(
    N = 2000, T = 200, M = 3,
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
    intercept = c("none", "global", "outcome"),
    standardize_outcomes = FALSE,
    outcome_weights = NULL,
    time_weights = NULL,
    eps_sd = 1e-8,
    test_ids = TRUE,
    tune_lambda = FALSE,
    tune_nu = FALSE,
    lambda_grid = 10^seq(-8, -1, by = 1),
    holdout = 5,
    tune_units = 120,
    parallel_backend = c("none", "fork", "psock"),
    n_cores = max(1L, parallel::detectCores() - 1L)
) {
  screen_method <- match.arg(screen_method)
  solver <- match.arg(solver)
  intercept <- match.arg(intercept)
  parallel_backend <- match.arg(parallel_backend)
  n_cores <- as.integer(n_cores)

  if (!is.null(seed)) set.seed(seed)

  # 1) Simulate data (simulator does NOT take intercept)
  sim <- simulate_multi_outcome_panel(N = N, T = T, M = M, K = K, L = L)

  # Stable unit IDs for downstream joins (and to test ID plumbing)
  unit_ids <- sprintf("U%06d", seq_len(N))

  treated_all <- which(is.finite(sim$treat_time))
  if (length(treated_all) < 1L) stop("Simulation produced no treated units.")
  treated_eval <- min(as.integer(treated_eval), length(treated_all))

  treated_units <- sample(treated_all, treated_eval)


  # 2) Optional tune_lambda() (uses same backend controls)
  tune <- NULL
  lambda_used <- lambda

  if (isTRUE(tune_lambda)) {
    if(verbose) message("Running lambda tuning step.\n")
    tune <- tune_lambda(
      Y_list = sim$Yobs,
      treat_time = sim$treat_time,
      treated_units = treated_units,
      L = L, K = K,
      lambda_grid = lambda_grid,
      holdout = holdout,
      intercept = intercept,
      max_donors = min(max_donors, 1000L),
      tune_units = tune_units,
      parallel_backend = parallel_backend,
      n_cores = n_cores
    )
    lambda_used <- tune$lambda_opt
    if (verbose) {
      message(sprintf("tune_lambda(): lambda_opt = %f\n", lambda_used))
    }
  }

  if(isTRUE(tune_nu)) {
    message("Running nu tuning step.\n")
    tm <- tune_nu(
      Y_list = sim$Yobs,
      treat_time = sim$treat_time,
      treated_units = treated_units,
      L = L, K = K,
      lambda = lambda_used,          # fixed lambda from tune_lambda()
      nu_grid = c(0, 1e-3, 1e-2, 1e-1, 1, 10),
      holdout = 5,
      intercept = intercept,
      max_donors = max_donors,
      tune_units = 120,
      parallel_backend = parallel_backend,
      n_cores = n_cores,
      verbose = verbose
    )

  nu_opt <- tm$nu_opt
  }

  # 3) Fit multiout_synth
  fit <- multiout_synth(
    Y_list = sim$Yobs,
    treat_time = sim$treat_time,
    treated_units = treated_units,
    unit_ids = unit_ids,
    L = L, K = K,
    max_donors = max_donors,
    screen_outcome = screen_outcome,
    screen_method = screen_method,
    lambda = lambda_used,
    solver = solver,
    pooled_adjustment = pooled_adjustment,
    nu = nu,
    verbose = verbose,
    standardize_outcomes = standardize_outcomes,
    outcome_weights = outcome_weights,
    time_weights = time_weights,
    intercept = intercept,
    eps_sd = eps_sd,
    parallel_backend = parallel_backend,
    n_cores = n_cores
  )

  # IMPORTANT: multiout_synth can filter treated_units for feasibility; use fit$treated_units
  treated_units_fit <- fit$treated_units
  J <- length(treated_units_fit)

  # 4) ID mapping self-test (robust to feasibility filtering)
  if (isTRUE(test_ids)) {
    stopifnot(identical(fit$unit_ids, unit_ids))
    stopifnot(identical(fit$treated_unit_ids, unit_ids[treated_units_fit]))

    check_jj <- seq_len(min(5L, J))
    for (jj in check_jj) {
      stopifnot(identical(fit$donor_ids[[jj]], unit_ids[fit$donors[[jj]]]))
      stopifnot(length(fit$donor_ids[[jj]]) == length(fit$weights[[jj]]))
    }
  }

  # 5) Accuracy for k = 0
  stopifnot(all(dim(fit$tau)[1:3] == c(J, M, K + 1)))

  # truth at k=0 for the same units used in fit
  tau_true_k0 <- matrix(NA_real_, nrow = J, ncol = M)

  if (!is.null(sim$tau_it) && length(dim(sim$tau_it)) == 3 && dim(sim$tau_it)[2] == M && dim(sim$tau_it)[3] >= (K + 1)) {
    tau_true_k0 <- sim$tau_it[treated_units_fit, , 1, drop = FALSE]
    tau_true_k0 <- matrix(tau_true_k0, nrow = J, ncol = M)
  } else if (!is.null(sim$tau_true_event) && length(dim(sim$tau_true_event)) == 3 && dim(sim$tau_true_event)[2] == M) {
    tau_true_k0 <- sim$tau_true_event[treated_units_fit, , 1, drop = FALSE]
    tau_true_k0 <- matrix(tau_true_k0, nrow = J, ncol = M)
  } else if (!is.null(sim$tau_unit)) {
    tau_true_k0 <- sim$tau_unit[treated_units_fit, , drop = FALSE]
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

  # 6) Average over event times 1..K evaluation (expects eval_avg_1K(fit, sim, K))
  avg_1K <- eval_avg_1K(fit = fit, sim = sim, K = K)

  # 7) Optional post-τ shrinkage layer (Option B)
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
    unit_ids = unit_ids,
    lambda_used = lambda_used,
    cor = cor_out,
    rmse = rmse_out,
    avg_1K = avg_1K,
    shrink = shrink,
    tune = tune
  )
}
