#' Tune pooled-adjustment ridge penalty nu by pre-period cross-validation
#'
#' Tunes \code{nu} (a.k.a. \code{nu}) for the multisynth-like pooled intercept adjustment.
#' For each candidate \code{nu}, we:
#' (i) fit pooled intercept adjustments \code{beta0[m]} using pooled pre-treatment residuals
#' across a subset of treated units, then (ii) evaluate validation error on held-out
#' pre-treatment periods after applying the pooled adjustment.
#'
#' @param Y_list List of outcome matrices, each \code{N x T}.
#' @param treat_time Numeric vector length \code{N}: first treated period (or \code{Inf}).
#' @param treated_units Integer vector of treated unit indices (or IDs if \code{unit_ids} provided).
#' @param L Integer. Total number of pre-treatment periods in the design.
#' @param K Integer. Post window size used only for donor eligibility in screening.
#' @param lambda Fixed ridge penalty used to fit donor weights during tuning.
#' @param nu_grid Numeric vector of candidate \code{nu} values (>=0).
#' @param holdout Integer. Number of pre periods to hold out (must be < \code{L}).
#' @param method Character. \code{"blocked"} holds out last \code{holdout} pre periods;
#'   \code{"random"} holds out \code{holdout} randomly (per unit).
#' @param screen_outcome Integer. Outcome index used for donor screening.
#' @param screen_method Character. \code{"cor"} or \code{"mse"}.
#' @param max_donors Integer. Maximum donors after screening.
#' @param standardize_outcomes Logical. Passed to \code{build_Xy_for_unit()}.
#' @param eps_sd Numeric. Passed to \code{build_Xy_for_unit()} to stabilize standardization.
#' @param outcome_weights Optional numeric length \code{M}. Outcome weights for objective.
#' @param time_weights Optional numeric length \code{L}. Time weights for objective.
#' @param intercept Character. One of \code{"none","global","outcome"} intercept handling.
#' @param solver Character. Currently \code{"fw"}.
#' @param seed Optional integer for reproducibility when \code{method="random"}.
#' @param verbose Logical. Print progress.
#' @param parallel_backend One of \code{"none","fork","psock"}.
#' @param n_cores Integer. Number of cores for parallel backends.
#' @param unit_ids Optional vector length \code{N} of unit IDs for mapping.
#' @param tune_units Integer. If provided, sample this many treated units for tuning.
#'
#' @return A list with:
#' \describe{
#'   \item{nu_opt}{Chosen nu (scalar) minimizing mean validation MSE across tuned units.}
#'   \item{cv}{Data.frame with columns \code{nu}, \code{mse_mean}, \code{mse_median}, \code{n_units}.}
#'   \item{nu_grid}{The searched grid.}
#'   \item{holdout}{Holdout size.}
#'   \item{method}{CV method.}
#'   \item{units_used}{Treated units used for tuning.}
#'   \item{beta0_by_nu}{Optional list of pooled intercept vectors (length M) for each nu.}
#' }
#'
#' @export
tune_nu <- function(
    Y_list,
    treat_time,
    treated_units,
    L,
    K,
    lambda = 1e-3,
    nu_grid = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
    holdout = 5,
    method = c("blocked", "random"),
    screen_outcome = 1,
    screen_method = c("cor", "mse"),
    max_donors = 1000,
    standardize_outcomes = FALSE,
    eps_sd = 1e-8,
    outcome_weights = NULL,
    time_weights = NULL,
    intercept = c("none", "global", "outcome"),
    solver = c("fw"),
    seed = 1,
    verbose = FALSE,
    parallel_backend = c("none","fork","psock"),
    n_cores = 2L,
    unit_ids = NULL,
    tune_units = 200,
    return_beta0 = FALSE
) {
  method <- match.arg(method)
  screen_method <- match.arg(screen_method)
  intercept <- match.arg(intercept)
  solver <- match.arg(solver)
  parallel_backend <- match.arg(parallel_backend)

  M <- length(Y_list)
  if (M < 1) stop("Y_list must contain at least one outcome matrix.")
  N <- nrow(Y_list[[1]])
  TT <- ncol(Y_list[[1]])

  for (m in seq_len(M)) {
    Y_list[[m]] <- as.matrix(Y_list[[m]])
    storage.mode(Y_list[[m]]) <- "double"
    if (nrow(Y_list[[m]]) != N || ncol(Y_list[[m]]) != TT) {
      stop("All outcome matrices in Y_list must have identical dimensions N x T.")
    }
  }
  if (length(treat_time) != N) stop("treat_time must have length N.")

  if (!is.numeric(L) || length(L) != 1 || L < 2) stop("L must be integer >= 2.")
  if (!is.numeric(K) || length(K) != 1 || K < 0) stop("K must be integer >= 0.")
  L <- as.integer(L); K <- as.integer(K)

  holdout <- as.integer(holdout)
  if (holdout < 1L || holdout >= L) stop("holdout must be >=1 and < L.")

  nu_grid <- as.numeric(nu_grid)
  if (length(nu_grid) < 1L || any(!is.finite(nu_grid)) || any(nu_grid < 0)) {
    stop("nu_grid must be non-empty, finite, and >= 0.")
  }

  # Map IDs -> indices if needed
  if (!is.null(unit_ids) && (is.character(treated_units) || is.factor(treated_units))) {
    treated_units <- match(as.character(treated_units), unit_ids)
    if (anyNA(treated_units)) stop("Some treated_units IDs were not found in unit_ids.")
  }
  treated_units <- as.integer(treated_units)
  treated_units <- treated_units[is.finite(treat_time[treated_units])]
  if (length(treated_units) < 5L) stop("Need at least a few treated units to tune nu.")

  # optional sub-sample for speed
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(tune_units) && is.finite(tune_units) && tune_units < length(treated_units)) {
    units_used <- sample(treated_units, as.integer(tune_units))
  } else {
    units_used <- treated_units
  }
  J <- length(units_used)

  # objective weights
  if (is.null(outcome_weights)) outcome_weights <- rep(1, M)
  if (length(outcome_weights) != M) stop("outcome_weights must have length M.")
  if (is.null(time_weights)) time_weights <- rep(1, L)
  if (length(time_weights) != L) stop("time_weights must have length L.")

  if (verbose) {
    message(sprintf("Tuning nu on %d treated units, grid size %d (backend=%s, cores=%d)",
                    J, length(nu_grid), parallel_backend, n_cores))
  }

  # ---- Stage A: per-unit prep: donors, weights (fixed lambda), and pre/val slices
  prep_one_unit <- function(j) {
    tryCatch({
      Tj <- treat_time[j]
      if (!is.finite(Tj)) stop("Unit j must be treated (finite treat_time[j]).")
      Tj <- as.integer(Tj)

      pre_times_full <- (Tj - L):(Tj - 1)
      if (min(pre_times_full) < 1L) stop("Not enough pre periods for unit j given L.")

      idx_full <- seq_len(L)
      idx_val <- if (method == "blocked") (L - holdout + 1L):L else sort(sample(idx_full, holdout))
      idx_tr  <- setdiff(idx_full, idx_val)

      pre_tr <- pre_times_full[idx_tr]
      pre_va <- pre_times_full[idx_val]

      donors <- donor_screen(
        Y_ref = Y_list[[screen_outcome]],
        treat_time = treat_time,
        j = j, K = K, L = L,
        max_donors = max_donors,
        method = screen_method
      )
      donors <- as.integer(donors)

      # build full design on L pre periods
      XY <- build_Xy_for_unit(
        Y_list, treat_time, j, donors, L,
        standardize_outcomes = standardize_outcomes,
        eps_sd = eps_sd
      )
      X_full <- XY$X  # (M*L) x Nd
      y_full <- XY$y  # length M*L

      # stacked row weights
      row_w_full <- rep(outcome_weights, each = L) * rep(time_weights, times = M)
      row_w_full <- as.numeric(row_w_full)

      # train rows are outcome-block replicated
      tr_rows <- unlist(lapply(seq_len(M), function(m) (m - 1L) * L + idx_tr))

      # train pieces for weight fitting
      X <- X_full[tr_rows, , drop = FALSE]
      y <- y_full[tr_rows]
      row_w <- row_w_full[tr_rows]

      # intercept handling inside objective (train only)
      if (intercept == "global") {
        sw <- sum(row_w)
        y <- y - sum(row_w * y) / sw
        Xbar <- as.numeric(crossprod(row_w, X) / sw)
        X <- X - matrix(Xbar, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
      } else if (intercept == "outcome") {
        tmp <- demean_within_outcome_blocks(X, y, L = length(idx_tr), time_weights = time_weights[idx_tr])
        X <- tmp$X; y <- tmp$y
      }

      # pre-whiten
      sqrtw <- sqrt(row_w)
      Xw <- X * sqrtw
      yw <- y * sqrtw

      # fit weights once (fixed lambda)
      w_hat <- switch(
        solver,
        fw = fit_weights_fw(Xw, yw, lambda = lambda, max_iter = 800, tol = 1e-5, verbose = FALSE)
      )
      w_hat <- as.numeric(w_hat)
      if (length(w_hat) != length(donors)) stop("Weight length mismatch with donors.")

      list(ok = TRUE, j = j, donors = donors, w = w_hat, pre_tr = pre_tr, pre_va = pre_va)
    }, error = function(e) {
      list(ok = FALSE, j = j, msg = conditionMessage(e))
    })
  }

  export_prep <- c(
    "prep_one_unit",
    "Y_list","treat_time","M","N","TT","L","K",
    "max_donors","screen_outcome","screen_method",
    "standardize_outcomes","eps_sd",
    "outcome_weights","time_weights","intercept",
    "lambda","solver","method","holdout",
    "donor_screen","build_Xy_for_unit",
    "demean_within_outcome_blocks","fit_weights_fw"
  )

  prep_list <- run_parallel(
    X = units_used,
    FUN = prep_one_unit,
    parallel_backend = parallel_backend,
    n_cores = n_cores,
    export = export_prep
  )

  prep_ok <- vapply(prep_list, function(z) is.list(z) && isTRUE(z$ok), logical(1))
  if (!any(prep_ok)) {
    msgs <- unique(vapply(prep_list, function(z) if (!isTRUE(z$ok)) as.character(z$msg)[1] else NA_character_, character(1)))
    msgs <- msgs[!is.na(msgs)]
    stop("No successful prep evaluations. Example errors: ", paste(head(msgs, 5), collapse = " | "))
  }

  prep_list <- prep_list[prep_ok]
  units_used_ok <- vapply(prep_list, `[[`, integer(1), "j")
  J_ok <- length(units_used_ok)

  # ---- Stage B: for each mu compute pooled beta0 (per outcome), then validation MSE per unit
  # Helper: compute pooled beta0 for a given mu
  compute_beta0_nu <- function(nu) {
    beta0 <- numeric(M)

    for (m in seq_len(M)) {
      y_pre_all <- numeric(0)
      ysc_pre_all <- numeric(0)
      Y <- Y_list[[m]]

      for (zz in prep_list) {
        j <- zz$j
        donors <- zz$donors
        w <- zz$w
        pre_tr <- zz$pre_tr

        y_pre <- Y[j, pre_tr]
        y_sc  <- as.numeric(crossprod(w, Y[donors, pre_tr, drop = FALSE]))

        y_pre_all <- c(y_pre_all, y_pre)
        ysc_pre_all <- c(ysc_pre_all, y_sc)
      }

      beta0[m] <- fit_pooled_intercept(y_pre_all, ysc_pre_all, nu = nu)
    }
    beta0
  }

  # Helper: compute validation MSE per unit for a given nu and beta0
  val_mse_one_unit <- function(zz, beta0) {
    j <- zz$j
    donors <- zz$donors
    w <- zz$w
    pre_va <- zz$pre_va

    mse_m <- numeric(M)
    for (m in seq_len(M)) {
      Y <- Y_list[[m]]
      yv <- Y[j, pre_va]
      ysc <- as.numeric(crossprod(w, Y[donors, pre_va, drop = FALSE]))
      # pooled adjustment: subtract beta0[m] (level shift)
      resid <- (yv - ysc - beta0[m])
      mse_m[m] <- mean(resid^2)
    }
    mean(mse_m)  # aggregate across outcomes
  }

  mse_by_nu <- numeric(length(nu_grid))
  beta0_by_nu <- if (isTRUE(return_beta0)) vector("list", length(nu_grid)) else NULL

  for (g in seq_along(nu_grid)) {
    nu <- nu_grid[g]
    beta0 <- compute_beta0_nu(nu)
    if (isTRUE(return_beta0)) beta0_by_nu[[g]] <- beta0

    mse_unit <- vapply(prep_list, function(zz) val_mse_one_unit(zz, beta0), numeric(1))
    mse_by_nu[g] <- mean(mse_unit, na.rm = TRUE)

    if (verbose) message(sprintf("nu=%g  mean_val_mse=%g", nu, mse_by_nu[g]))
  }

  nu_opt <- nu_grid[which.min(mse_by_nu)]

  cv <- data.frame(
    nu = nu_grid,
    mse_mean = as.numeric(mse_by_nu),
    mse_median = as.numeric(mse_by_nu),  # optional: could store unit-level medians if you want
    n_units = J_ok
  )

  out <- list(
    nu_opt = nu_opt,
    cv = cv,
    nu_grid = nu_grid,
    holdout = holdout,
    method = method,
    units_used = units_used_ok
  )
  if (isTRUE(return_beta0)) out$beta0_by_nu <- beta0_by_nu
  out
}
