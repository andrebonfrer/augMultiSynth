#' Multi-outcome synthetic control with unit-level treatment effects
#'
#' Estimates donor weights and unit-level event-time treatment effects for multiple
#' outcomes in staggered adoption settings. Designed to be scalable via donor screening
#' and simplex solvers. Optionally includes a pooled regression-adjustment intercept
#' (multisynth-like) estimated from pooled pre-treatment residuals across treated units.
#'
#' @param Y_list List of outcome matrices, each \code{N x T}. All must share \code{N} and \code{T}.
#' @param treat_time Numeric vector length \code{N} giving first treated period for each unit
#'   (or \code{Inf} for never-treated).
#' @param treated_units Optional integer vector of treated unit indices to estimate.
#'   Defaults to all units with finite \code{treat_time}.
#' @param L Integer. Number of pre-treatment periods used to fit weights and (if enabled)
#'   compute intercept shifts.
#' @param K Integer. Maximum event time; returns effects for \code{k=0,...,K}.
#' @param max_donors Integer. Maximum donor set size after screening. Defaults to \code{1000}.
#' @param screen_outcome Integer. Index in \code{Y_list} to use for donor screening.
#'   Defaults to \code{1}.
#' @param screen_method Screening method for donors. \code{"cor"} or \code{"mse"}.
#' @param demean Logical. If \code{TRUE}, demeans pre-treatment stacked design before fitting.
#'   Recommended when using intercept-shift interpretation.
#' @param lambda Nonnegative scalar ridge penalty on weights in the solver.
#' @param solver Solver choice. Currently \code{"fw"} (Frank–Wolfe).
#' @param pooled_adjustment Logical. If \code{TRUE}, estimates a pooled intercept correction
#'   per outcome from pre-treatment residuals across treated units (multisynth-like).
#' @param nu Nonnegative ridge penalty for the pooled intercept correction (per outcome).
#'   Larger \code{nu} shrinks the pooled intercept toward zero.
#' @param verbose Logical. If \code{TRUE}, prints progress.
#'
#' @return A list with components:
#' \describe{
#'   \item{weights}{A named list of donor weight vectors, one per treated unit.}
#'   \item{donors}{A named list of donor index vectors, one per treated unit.}
#'   \item{tau}{A 3D array of dimension \code{J x M x (K+1)} containing treatment effects,
#'     where \code{J = length(treated_units)} and \code{M = length(Y_list)}.}
#'   \item{treated_units}{Integer vector of treated units used in estimation.}
#'   \item{pooled_beta0}{NULL or numeric vector length \code{M} of pooled intercept adjustments.}
#'   \item{pooled_adjustment}{Logical flag as passed in.}
#'   \item{nu}{Penalty parameter used for pooled adjustment.}
#' }
#'
#' @export
multiout_synth <- function(Y_list, treat_time, treated_units = NULL,
                           L, K,
                           max_donors = 1000,
                           screen_outcome = 1,
                           screen_method = c("cor", "mse"),
                           demean = TRUE,
                           lambda = 1e-3,
                           solver = c("fw"),
                           pooled_adjustment = FALSE,
                           nu = 0,
                           verbose = FALSE,
                           standardize_outcomes = FALSE,
                           eps_sd = 1e-8) {

  screen_method <- match.arg(screen_method)
  solver <- match.arg(solver)

  M <- length(Y_list)
  if (M < 1) stop("Y_list must contain at least one outcome matrix.")
  if (!is.matrix(Y_list[[1]])) stop("Each element of Y_list must be a matrix.")
  N <- nrow(Y_list[[1]])
  TT <- ncol(Y_list[[1]])

  # Ensure all outcomes share dimensions
  for (m in seq_len(M)) {
    if (!is.matrix(Y_list[[m]])) stop("Each element of Y_list must be a matrix.")
    if (nrow(Y_list[[m]]) != N || ncol(Y_list[[m]]) != TT) {
      stop("All outcome matrices in Y_list must have identical dimensions N x T.")
    }
  }

  if (length(treat_time) != N) stop("treat_time must have length N.")

  if (is.null(treated_units)) {
    treated_units <- which(is.finite(treat_time))
  } else {
    treated_units <- as.integer(treated_units)
  }
  J <- length(treated_units)
  if (J < 1) stop("No treated units supplied or found.")

  if (!is.numeric(L) || length(L) != 1 || L < 1) stop("L must be a positive integer.")
  if (!is.numeric(K) || length(K) != 1 || K < 0) stop("K must be a nonnegative integer.")
  L <- as.integer(L)
  K <- as.integer(K)

  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) stop("lambda must be >= 0.")
  if (!is.numeric(nu) || length(nu) != 1 || nu < 0) stop("nu must be >= 0.")

  tau <- array(NA_real_, dim = c(J, M, K + 1))
  weights <- vector("list", J)
  donors_list <- vector("list", J)
  names(weights) <- names(donors_list) <- as.character(treated_units)

  # -----------------------------
  # PASS 1: Fit weights per treated unit
  # -----------------------------
  for (jj in seq_len(J)) {
    j <- treated_units[jj]
    if (verbose && (jj %% 50 == 0)) {
      message(sprintf("Fitting weights for treated unit %d of %d", jj, J))
    }

    donors <- donor_screen(
      Y_ref = Y_list[[screen_outcome]],
      treat_time = treat_time,
      j = j, K = K, L = L,
      max_donors = max_donors,
      method = screen_method
    )

    XY <- build_Xy_for_unit(Y_list, treat_time, j, donors, L,
                            standardize_outcomes, eps_sd)
    X <- XY$X
    y <- XY$y

    if (demean) {
      X <- demean_rows(X)
      y <- as.numeric(demean_rows(matrix(y, nrow = 1)))
    }

    w <- switch(
      solver,
      fw = fit_weights_fw(X, y, lambda = lambda, max_iter = 2000, tol = 1e-6, verbose = FALSE)
    )

    weights[[jj]] <- w
    donors_list[[jj]] <- donors
  }

  # -----------------------------
  # PASS 2: Optional pooled intercept adjustment (per outcome)
  # -----------------------------
  pooled_beta0 <- NULL
  if (isTRUE(pooled_adjustment)) {
    pooled_beta0 <- numeric(M)

    for (m in seq_len(M)) {
      y_pre_all <- numeric(0)
      ysc_pre_all <- numeric(0)

      Y <- Y_list[[m]]

      for (jj in seq_len(J)) {
        i <- treated_units[jj]
        Ti <- treat_time[i]
        if (!is.finite(Ti)) next
        Ti <- as.integer(Ti)

        pre_times <- (Ti - L):(Ti - 1)
        if (min(pre_times) < 1) stop("Not enough pre periods for pooled adjustment given L.")

        donors <- donors_list[[jj]]
        w <- weights[[jj]]

        y_i_pre <- Y[i, pre_times]
        y_sc_pre <- as.numeric(crossprod(w, Y[donors, pre_times, drop = FALSE]))

        y_pre_all <- c(y_pre_all, y_i_pre)
        ysc_pre_all <- c(ysc_pre_all, y_sc_pre)
      }

      pooled_beta0[m] <- fit_pooled_intercept(y_pre_all, ysc_pre_all, nu = nu)
    }
  }

  # -----------------------------
  # PASS 3: Compute unit effects using fitted weights (and optional pooled_beta0)
  # -----------------------------
  for (jj in seq_len(J)) {
    j <- treated_units[jj]

    tau_j <- unit_effects(
      Y_list = Y_list,
      treat_time = treat_time,
      j = j,
      donors = donors_list[[jj]],
      w = weights[[jj]],
      K = K,
      use_intercept = TRUE,
      L = L,
      pooled_beta0 = pooled_beta0
    )

    tau[jj, , ] <- tau_j
  }

  list(
    weights = weights,
    donors = donors_list,
    tau = tau,
    treated_units = treated_units,
    pooled_beta0 = pooled_beta0,
    pooled_adjustment = pooled_adjustment,
    nu = nu
  )
}
