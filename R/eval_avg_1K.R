#' Evaluate avg treatment effects over event times 1..K
#'
#' @param fit Fit object returned by multiout_synth()
#' @param sim Simulation object returned by simulate_multi_outcome_panel()
#' @param K Integer. Max event time used.
#'
#' @return List with tau_hat_avg_1K, tau_true_avg_1K, cor_avg_1K, rmse_avg_1K
#' @export
eval_avg_1K <- function(fit, sim, K) {
  tau_hat <- fit$tau                 # J x M x (K+1)
  if (is.null(dim(tau_hat)) || length(dim(tau_hat)) != 3) {
    stop("fit$tau must be a 3D array J x M x (K+1).")
  }

  J <- dim(tau_hat)[1]
  M <- dim(tau_hat)[2]
  Khat <- dim(tau_hat)[3] - 1L

  treated_units <- fit$treated_units
  if (is.null(treated_units)) stop("fit$treated_units is NULL.")
  treated_units <- as.integer(treated_units)

  if (length(treated_units) != J) stop("fit$treated_units length mismatch with fit$tau.")
  if (!is.matrix(sim$tau_unit)) stop("sim$tau_unit not found or not a matrix; cannot compute truth.")
  if (ncol(sim$tau_unit) != M) stop("sim$tau_unit column count does not match number of outcomes M.")
  if (K < 1) stop("K must be >= 1 for avg_1K.")
  if (K > Khat) stop("Requested K exceeds the K used to compute fit$tau.")
  if (anyNA(treated_units)) stop("fit$treated_units contains NA.")
  if (any(treated_units < 1L)) stop("fit$treated_units contains indices < 1.")
  if (any(treated_units > nrow(sim$tau_unit))) stop("fit$treated_units contains indices > N in sim$tau_unit.")

  # Estimated avg over event times 1..K (exclude k=0)
  tau_hat_avg_1K <- apply(tau_hat[, , 2:(K + 1), drop = FALSE], c(1, 2), mean)
  tau_hat_avg_1K <- matrix(tau_hat_avg_1K, nrow = J, ncol = M)

  # True avg over event times 1..K
  tau_true_avg_1K <- sim$tau_unit[treated_units, , drop = FALSE]

  if (isTRUE(sim$dynamic_te)) {
    # Prefer stored path if you add it to simulator; otherwise reconstruct as in simulator
    if (!is.null(sim$tau_path_shape)) {
      tau_path_shape <- sim$tau_path_shape
    } else {
      ramp <- function(k) if (k < 5) (k + 1) / 5 else 1
      TT <- nrow(sim$f_t)
      tau_path_shape <- sapply(0:(TT - 1), ramp)
    }

    mult <- mean(tau_path_shape[2:(K + 1)])
    tau_true_avg_1K <- tau_true_avg_1K * mult
  }

  # Simple diagnostic: if columns are identical, stop with a short message
  if (M >= 2 && isTRUE(all.equal(tau_true_avg_1K[, 1], tau_true_avg_1K[, 2]))) {
    stop("tau_true_avg_1K columns are identical. This indicates mis-indexing or stale code being used.")
  }

  cor_avg_1K <- rmse_avg_1K <- rep(NA_real_, M)
  for (m in seq_len(M)) {
    ok <- is.finite(tau_true_avg_1K[, m]) & is.finite(tau_hat_avg_1K[, m])
    cor_avg_1K[m] <- if (sum(ok) >= 3) stats::cor(tau_true_avg_1K[ok, m], tau_hat_avg_1K[ok, m]) else NA_real_
    rmse_avg_1K[m] <- if (sum(ok) >= 1) sqrt(mean((tau_hat_avg_1K[ok, m] - tau_true_avg_1K[ok, m])^2)) else NA_real_
  }

  list(
    tau_hat_avg_1K = tau_hat_avg_1K,
    tau_true_avg_1K = tau_true_avg_1K,
    cor_avg_1K = cor_avg_1K,
    rmse_avg_1K = rmse_avg_1K,
    treated_units_used = treated_units
  )
}
