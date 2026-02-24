#' Tune weight ridge penalty lambda by pre-period cross-validation
#'
#' Selects \code{lambda} for the donor weight solver by minimizing pre-treatment
#' validation error on held-out pre periods. This tunes the ridge penalty in the
#' simplex weight problem:
#' \deqn{\min_{w\in\Delta}\ \|y - Xw\|_W^2 + \lambda \|w\|_2^2}
#' where \code{y, X} are stacked multi-outcome pre-treatment data (after any
#' intercept projection and row-weight pre-whitening).
#'
#' @param Y_list List of outcome matrices, each \code{N x T}.
#' @param treat_time Numeric vector length \code{N}: first treated period (or \code{Inf}).
#' @param treated_units Integer vector of treated unit indices (or IDs if \code{unit_ids} provided).
#' @param L Integer. Total number of pre-treatment periods in the design.
#' @param K Integer. Post window size used only for donor eligibility in screening.
#' @param lambda_grid Numeric vector of candidate \code{lambda} values (>=0).
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
#' @param max_iter Integer. Max iterations for solver inside tuning (keep lower than main fit).
#' @param tol Numeric. Tolerance for solver.
#' @param seed Optional integer for reproducibility when \code{method="random"}.
#' @param verbose Logical. Print progress.
#' @param parallel_backend One of \code{"none","fork","psock"}.
#' @param n_cores Integer. Number of workers/cores.
#' @param unit_ids Optional vector length \code{N} of unit IDs for mapping.
#' @param tune_units Integer. If provided, sample this many treated units for tuning.
#'
#' @return A list with:
#' \describe{
#'   \item{lambda_opt}{Chosen lambda minimizing mean validation MSE across tuned units.}
#'   \item{cv}{Data.frame with columns \code{lambda}, \code{mse_mean}, \code{mse_median}, \code{n_units}.}
#'   \item{lambda_grid}{The searched grid.}
#'   \item{holdout}{Holdout size.}
#'   \item{method}{CV method.}
#'   \item{units_used}{Treated units used for tuning (indices).}
#' }
#'
#' @export
tune_lambda <- function(
    Y_list,
    treat_time,
    treated_units,
    L,
    K,
    lambda_grid = 10^seq(-8, -2, by = 1),
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
    max_iter = 800,
    tol = 1e-5,
    seed = 1,
    verbose = FALSE,
    parallel_backend = c("none", "fork", "psock"),
    n_cores = 2L,
    unit_ids = NULL,
    tune_units = 200
) {
  method <- match.arg(method)
  screen_method <- match.arg(screen_method)
  intercept <- match.arg(intercept)
  solver <- match.arg(solver)
  parallel_backend <- match.arg(parallel_backend)
  n_cores <- as.integer(n_cores)

  # ---- Dimensions & coercions ----
  M <- length(Y_list)
  if (M < 1) stop("Y_list must contain at least one outcome matrix.")
  for (m in seq_len(M)) {
    Y_list[[m]] <- as.matrix(Y_list[[m]])
    storage.mode(Y_list[[m]]) <- "double"
  }
  N <- nrow(Y_list[[1]])
  TT <- ncol(Y_list[[1]])
  for (m in seq_len(M)) {
    if (nrow(Y_list[[m]]) != N || ncol(Y_list[[m]]) != TT) {
      stop("All outcome matrices in Y_list must have identical dimensions N x T.")
    }
  }
  if (length(treat_time) != N) stop("treat_time must have length N.")

  # ---- Validate scalars ----
  if (!is.numeric(L) || length(L) != 1 || L < 2) stop("L must be integer >= 2.")
  if (!is.numeric(K) || length(K) != 1 || K < 0) stop("K must be integer >= 0.")
  L <- as.integer(L); K <- as.integer(K)

  if (!is.numeric(holdout) || length(holdout) != 1 || holdout < 1) stop("holdout must be >= 1.")
  holdout <- as.integer(holdout)
  if (holdout >= L) stop("holdout must be < L.")

  lambda_grid <- as.numeric(lambda_grid)
  if (length(lambda_grid) < 1L) stop("lambda_grid must be non-empty.")
  if (any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) stop("lambda_grid must be finite and >= 0.")

  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 10) stop("max_iter must be a reasonable integer.")
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) stop("tol must be > 0.")

  if (!is.null(seed)) set.seed(seed)

  # ---- Map treated IDs -> indices if needed ----
  if (!is.null(unit_ids) && (is.character(treated_units) || is.factor(treated_units))) {
    treated_units <- match(as.character(treated_units), unit_ids)
    if (anyNA(treated_units)) stop("Some treated_units IDs were not found in unit_ids.")
  }
  treated_units <- as.integer(treated_units)
  treated_units <- treated_units[is.finite(treat_time[treated_units])]
  treated_units <- treated_units[treated_units >= 1L & treated_units <= N]
  if (length(treated_units) < 5L) stop("Need at least a few treated units to tune lambda.")

  # ---- Optional subsample for speed ----
  if (!is.null(tune_units) && is.finite(tune_units) && tune_units < length(treated_units)) {
    units_used <- sample(treated_units, as.integer(tune_units))
  } else {
    units_used <- treated_units
  }
  units_used <- unique(units_used)

  # ---- Objective weights ----
  ow <- outcome_weights
  if (is.null(ow)) ow <- rep(1, M)
  if (length(ow) != M) stop("outcome_weights must have length M.")
  if (any(!is.finite(ow)) || any(ow < 0) || sum(ow) <= 0) stop("outcome_weights must be finite, nonnegative, sum>0.")

  tw <- time_weights
  if (is.null(tw)) tw <- rep(1, L)
  if (length(tw) != L) stop("time_weights must have length L.")
  if (any(!is.finite(tw)) || any(tw < 0) || sum(tw) <= 0) stop("time_weights must be finite, nonnegative, sum>0.")

  row_w_full <- as.numeric(rep(ow, each = L) * rep(tw, times = M))  # length M*L

  if (verbose) {
    message(sprintf("tune_lambda(): units=%d, grid=%d, backend=%s, cores=%d",
                    length(units_used), length(lambda_grid), parallel_backend, n_cores))
  }

  # ---- Per-unit CV worker ----
  one_unit_cv <- function(j) {
    tryCatch({
      Tj <- treat_time[j]
      if (!is.finite(Tj)) stop("Unit j must be treated (finite treat_time[j]).")
      Tj <- as.integer(Tj)

      pre_times_full <- (Tj - L):(Tj - 1)
      if (min(pre_times_full) < 1L) stop("Not enough pre periods for unit j given L.")

      idx_full <- seq_len(L)
      idx_val <- if (method == "blocked") {
        (L - holdout + 1L):L
      } else {
        sort(sample(idx_full, holdout))
      }
      idx_tr <- setdiff(idx_full, idx_val)
      if (length(idx_tr) < 1L) stop("Training set empty after holdout; reduce holdout.")

      # donors screening (eligibility enforced inside donor_screen)
      donors <- donor_screen(
        Y_ref = Y_list[[screen_outcome]],
        treat_time = treat_time,
        j = j, K = K, L = L,
        max_donors = max_donors,
        method = screen_method
      )
      donors <- as.integer(donors)
      donors <- donors[is.finite(donors) & donors >= 1L & donors <= N]
      if (length(donors) < 2L) stop("Too few donors after screening.")

      # Build full stacked design on full pre window
      XY <- build_Xy_for_unit(
        Y_list, treat_time, j, donors, L,
        standardize_outcomes = standardize_outcomes,
        eps_sd = eps_sd
      )
      X_full <- XY$X   # (M*L) x Nd
      y_full <- XY$y   # length M*L

      # stacked row indices for train/valid
      tr_rows <- unlist(lapply(seq_len(M), function(m) (m - 1L) * L + idx_tr))
      va_rows <- unlist(lapply(seq_len(M), function(m) (m - 1L) * L + idx_val))

      # training slices
      X_tr <- X_full[tr_rows, , drop = FALSE]
      y_tr <- y_full[tr_rows]
      w_tr <- row_w_full[tr_rows]

      # apply intercept projection on training
      if (intercept == "global") {
        sw <- sum(w_tr)
        ybar <- sum(w_tr * y_tr) / sw
        Xbar <- as.numeric(crossprod(w_tr, X_tr) / sw)  # 1 x Nd

        y_tr <- y_tr - ybar
        X_tr <- X_tr - matrix(Xbar, nrow = nrow(X_tr), ncol = ncol(X_tr), byrow = TRUE)
      } else if (intercept == "outcome") {
        # training has M blocks of length length(idx_tr)
        tmp <- demean_within_outcome_blocks(X_tr, y_tr, L = length(idx_tr), time_weights = tw[idx_tr])
        X_tr <- tmp$X
        y_tr <- tmp$y
      }

      # pre-whiten training
      sqrtw_tr <- sqrt(w_tr)
      Xw_tr <- X_tr * sqrtw_tr
      yw_tr <- y_tr * sqrtw_tr

      # validation slices
      X_va <- X_full[va_rows, , drop = FALSE]
      y_va <- y_full[va_rows]
      w_va <- row_w_full[va_rows]

      # apply *compatible* intercept transform to validation
      if (intercept == "global") {
        # use TRAINING ybar/Xbar
        y_va <- y_va - ybar
        X_va <- X_va - matrix(Xbar, nrow = nrow(X_va), ncol = ncol(X_va), byrow = TRUE)
      } else if (intercept == "outcome") {
        # validation has M blocks of length length(idx_val)
        tmpv <- demean_within_outcome_blocks(X_va, y_va, L = length(idx_val), time_weights = tw[idx_val])
        X_va <- tmpv$X
        y_va <- tmpv$y
      }

      # pre-whiten validation
      sqrtw_va <- sqrt(w_va)
      Xw_va <- X_va * sqrtw_va
      yw_va <- y_va * sqrtw_va

      # evaluate lambdas
      mse <- numeric(length(lambda_grid))
      for (g in seq_along(lambda_grid)) {
        lam <- lambda_grid[g]

        w_hat <- switch(
          solver,
          fw = fit_weights_fw(Xw_tr, yw_tr, lambda = lam, max_iter = max_iter, tol = tol, verbose = FALSE)
        )

        # defensive coercion (in case solver returns list)
        if (is.list(w_hat)) {
          if (!is.null(w_hat$w)) w_hat <- w_hat$w
          else if (!is.null(w_hat$weights)) w_hat <- w_hat$weights
          else stop("fit_weights_fw() returned a list but no $w/$weights.")
        }
        w_hat <- as.numeric(w_hat)
        if (length(w_hat) != ncol(Xw_tr)) stop("Weight length mismatch in tuning.")

        resid <- yw_va - as.numeric(Xw_va %*% w_hat)
        mse[g] <- mean(resid^2)
      }

      list(ok = TRUE, j = j, mse = mse)
    }, error = function(e) {
      list(ok = FALSE, j = j, msg = conditionMessage(e))
    })
  }

  # ---- Run CV in chosen backend via run_parallel() ----
  # run_parallel(X, FUN, backend, n_cores) should return a list of length(X)
  out_list <- run_parallel(
    X = units_used,
    FUN = one_unit_cv,
    parallel_backend = parallel_backend,
    n_cores = n_cores,
    export = c(
      "one_unit_cv",
      "Y_list","treat_time","M","N","TT","L","K",
      "max_donors","screen_outcome","screen_method",
      "standardize_outcomes","eps_sd",
      "outcome_weights","time_weights","intercept",
      "lambda_grid","solver",
      "donor_screen","build_Xy_for_unit",
      "demean_within_outcome_blocks","fit_weights_fw"
    ),
    envir = environment()
  )

  # Normalize outputs (handles NULLs / weird returns)
  out_list <- lapply(out_list, function(z) {
    if (is.null(z)) return(list(ok = FALSE, j = NA_integer_, msg = "worker returned NULL"))
    if (!is.list(z)) return(list(ok = FALSE, j = NA_integer_, msg = paste("unexpected return:", paste(class(z), collapse = "/"))))
    if (is.null(z$ok) || length(z$ok) == 0) z$ok <- FALSE
    if (is.null(z$msg) || length(z$msg) == 0) z$msg <- NA_character_
    z$msg <- as.character(z$msg)[1]
    z
  })

  ok <- vapply(out_list, function(z) isTRUE(z$ok), logical(1))

  if (!any(ok)) {
    msgs <- unique(vapply(out_list, function(z) if (!isTRUE(z$ok)) z$msg else NA_character_, character(1)))
    msgs <- msgs[!is.na(msgs)]
    stop("No successful CV evaluations. Example errors: ", paste(head(msgs, 6), collapse = " | "))
  }

  mse_mat <- do.call(rbind, lapply(out_list[ok], function(z) z$mse))
  if (!is.matrix(mse_mat) || ncol(mse_mat) != length(lambda_grid)) {
    stop("Internal error: mse_mat shape mismatch. Check one_unit_cv returns mse with length(lambda_grid).")
  }

  mse_mean <- colMeans(mse_mat, na.rm = TRUE)
  mse_median <- apply(mse_mat, 2, stats::median, na.rm = TRUE)

  lambda_opt <- lambda_grid[which.min(mse_mean)]

  cv <- data.frame(
    lambda = lambda_grid,
    mse_mean = as.numeric(mse_mean),
    mse_median = as.numeric(mse_median),
    n_units = sum(ok)
  )

  list(
    lambda_opt = lambda_opt,
    cv = cv,
    lambda_grid = lambda_grid,
    holdout = holdout,
    method = method,
    units_used = units_used
  )
}
