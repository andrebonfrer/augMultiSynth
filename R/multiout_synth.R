#' Multi-outcome synthetic control with unit-level treatment effects
#'
#' Fits simplex-constrained donor weights for each treated unit by minimizing a
#' stacked, weighted pre-treatment objective across outcomes, then computes
#' unit-level event-time treatment effects for k = 0..K. Optionally estimates a
#' pooled intercept adjustment per outcome (multisynth-like) from pooled
#' pre-treatment residuals across treated units.
#'
#' @param Y_list List of outcome matrices (N x T), all sharing the same N and T.
#' @param treat_time Numeric vector length N; first treated period per unit (Inf = never-treated).
#' @param treated_units Optional integer vector of treated unit indices to estimate.
#'   Defaults to all finite treat_time.
#' @param unit_ids Optional vector length N of unit identifiers (character/integer/etc.).
#'   If NULL, uses 1:N.
#' @param L Integer; number of pre-treatment periods used for weight fitting.
#' @param K Integer; maximum event time (k = 0..K).
#' @param max_donors Integer; cap on donor count after screening.
#' @param screen_outcome Integer; which outcome to use for donor screening.
#' @param screen_method "cor" or "mse".
#' @param lambda Nonnegative ridge penalty on weights in the solver.
#' @param solver Currently "fw".
#' @param pooled_adjustment Logical; if TRUE fit pooled intercept adjustment per outcome.
#' @param nu Nonnegative ridge penalty for pooled adjustment (per outcome).
#' @param verbose Logical.
#' @param standardize_outcomes Logical; standardize each outcome block in X/y for weight fitting only.
#' @param outcome_weights Optional length M weights on outcome blocks in the objective.
#' @param time_weights Optional length L weights on pre-periods (shared across outcomes).
#' @param intercept One of "none","global","outcome" controlling intercept projection in objective.
#' @param eps_sd Small jitter for standardization stability.
#' @param parallel_backend One of "none","fork","psock". Recommended on macOS: "psock".
#' @param n_cores Integer; number of cores/workers when parallel_backend != "none".
#'
#' @return A list with weights, donors, donor_ids, tau, treated_units, treated_unit_ids,
#'   unit_ids, pooled_beta0, pooled_adjustment, nu.
#'
#' @export
multiout_synth <- function(
    Y_list, treat_time,
    treated_units = NULL,
    unit_ids = NULL,
    L, K,
    max_donors = 1000,
    screen_outcome = 1,
    screen_method = c("cor", "mse"),
    lambda = 1e-3,
    solver = c("fw"),
    pooled_adjustment = FALSE,
    nu = 0,
    verbose = FALSE,
    standardize_outcomes = FALSE,
    outcome_weights = NULL,
    time_weights = NULL,
    intercept = c("none", "global", "outcome"),
    eps_sd = 1e-8,
    parallel_backend = c("none", "fork", "psock"),
    n_cores = max(1L, parallel::detectCores() - 1L)
) {
  # -----------------------------
  # Parse args / basic checks
  # -----------------------------
  screen_method <- match.arg(screen_method)
  solver <- match.arg(solver)
  intercept <- match.arg(intercept)
  parallel_backend <- match.arg(parallel_backend)

  if (!is.list(Y_list) || length(Y_list) < 1L) stop("Y_list must be a non-empty list of matrices.")
  if (!is.numeric(treat_time)) stop("treat_time must be numeric.")
  if (!is.numeric(L) || length(L) != 1L || L < 1) stop("L must be a positive integer.")
  if (!is.numeric(K) || length(K) != 1L || K < 0) stop("K must be a nonnegative integer.")
  L <- as.integer(L); K <- as.integer(K)

  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0) stop("lambda must be >= 0.")
  if (!is.numeric(nu) || length(nu) != 1L || nu < 0) stop("nu must be >= 0.")

  # -----------------------------
  # Coerce outcomes to double matrices; check dimensions
  # -----------------------------
  M <- length(Y_list)
  for (m in seq_len(M)) {
    Y_list[[m]] <- as.matrix(Y_list[[m]])
    storage.mode(Y_list[[m]]) <- "double"
  }

  N <- nrow(Y_list[[1]])
  TT <- ncol(Y_list[[1]])

  for (m in seq_len(M)) {
    if (!is.matrix(Y_list[[m]])) stop("Each element of Y_list must be a matrix.")
    if (nrow(Y_list[[m]]) != N || ncol(Y_list[[m]]) != TT) {
      stop("All outcome matrices in Y_list must share identical N and T.")
    }
  }

  if (length(treat_time) != N) stop("treat_time must have length N (nrow of outcomes).")

  # -----------------------------
  # Unit identifiers
  # -----------------------------
  if (is.null(unit_ids)) {
    unit_ids <- seq_len(N)
  } else {
    if (length(unit_ids) != N) stop("unit_ids must have length N.")
  }

  # -----------------------------
  # Determine treated units and apply feasibility filter
  # -----------------------------
  if (is.null(treated_units)) {
    treated_units <- which(is.finite(treat_time))
  } else {
    treated_units <- as.integer(treated_units)
  }

  # keep only truly treated
  treated_units <- treated_units[is.finite(treat_time[treated_units])]
  if (length(treated_units) < 1L) stop("No treated units with finite treat_time.")

  # feasibility: L pre periods and K post periods exist
  Ti <- as.integer(treat_time[treated_units])
  ok_feas <- (Ti - L >= 1L) & (Ti + K <= TT)
  treated_units <- treated_units[ok_feas]

  if (length(treated_units) < 1L) {
    stop("No treated units satisfy feasibility constraints given L and K.")
  }

  treated_unit_ids <- unit_ids[treated_units]
  J <- length(treated_units)

  # -----------------------------
  # Outcome/time weights for stacked objective
  # -----------------------------
  ow <- outcome_weights
  if (is.null(ow)) ow <- rep(1, M)
  if (length(ow) != M) stop("outcome_weights must have length M.")
  if (any(!is.finite(ow)) || any(ow < 0) || sum(ow) <= 0) {
    stop("outcome_weights must be finite, nonnegative, and sum to > 0.")
  }

  tw <- time_weights
  if (is.null(tw)) tw <- rep(1, L)
  if (length(tw) != L) stop("time_weights must have length L.")
  if (any(!is.finite(tw)) || any(tw < 0) || sum(tw) <= 0) {
    stop("time_weights must be finite, nonnegative, and sum to > 0.")
  }

  # precompute row weights for stacked (M blocks of length L)
  row_w_base <- rep(ow, each = L) * rep(tw, times = M)
  row_w_base <- as.numeric(row_w_base)

  # -----------------------------
  # Allocate outputs (named by treated IDs)
  # -----------------------------
  nm <- as.character(treated_unit_ids)
  weights <- vector("list", J); names(weights) <- nm
  donors_list <- vector("list", J); names(donors_list) <- nm
  donor_ids_list <- vector("list", J); names(donor_ids_list) <- nm
  tau <- array(NA_real_, dim = c(J, M, K + 1))

  # -----------------------------
  # Worker: fit weights for one treated unit
  # -----------------------------
  fit_one_unit <- function(j) {
    tryCatch({
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

      XY <- build_Xy_for_unit(
        Y_list = Y_list,
        treat_time = treat_time,
        j = j,
        donors = donors,
        L = L,
        standardize_outcomes = standardize_outcomes,
        eps_sd = eps_sd
      )
      X <- XY$X
      y <- XY$y

      # row weights (copy base)
      row_w <- row_w_base
      if (length(y) != length(row_w) || nrow(X) != length(row_w)) {
        stop("Internal error: stacked dimension mismatch for X/y/row weights.")
      }

      # intercept projection inside objective
      if (intercept == "global") {
        sw <- sum(row_w)
        y <- y - sum(row_w * y) / sw
        Xbar <- as.numeric(crossprod(row_w, X) / sw)  # 1 x Nd
        X <- X - matrix(Xbar, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
      } else if (intercept == "outcome") {
        tmp <- demean_within_outcome_blocks(X, y, L = L, time_weights = tw)
        X <- tmp$X
        y <- tmp$y
      }

      # apply W via pre-whitening
      sqrtw <- sqrt(row_w)
      Xw <- X * sqrtw
      yw <- y * sqrtw

      w <- switch(
        solver,
        fw = fit_weights_fw(Xw, yw, lambda = lambda, max_iter = 1000, tol = 1e-5, verbose = FALSE)
      )

      # normalize possible solver return types
      if (is.list(w)) {
        if (!is.null(w$w)) w <- w$w
        else if (!is.null(w$weights)) w <- w$weights
        else stop("Solver returned list without $w/$weights.")
      }
      w <- as.numeric(w)

      if (length(w) != length(donors)) stop("Weight length mismatch with donors.")
      if (any(!is.finite(w))) stop("Non-finite weights returned by solver.")

      list(ok = TRUE, j = j, donors = donors, w = w)
    }, error = function(e) {
      list(ok = FALSE, j = j, msg = conditionMessage(e))
    })
  }

  # -----------------------------
  # PASS 1: Fit weights (parallel via run_parallel)
  # -----------------------------
  if (verbose) {
    message(sprintf("PASS 1: fitting weights for J=%d treated units (backend=%s, cores=%d)",
                    J, parallel_backend, as.integer(n_cores)))
  }

  out_list <- run_parallel(
    X = treated_units,
    FUN = fit_one_unit,
    parallel_backend = parallel_backend,
    n_cores = as.integer(n_cores),
    # PSOCK needs exports; fork/none ignore export safely
    export = c(
      "fit_one_unit",
      "Y_list", "treat_time", "N", "TT", "M", "L", "K",
      "max_donors", "screen_outcome", "screen_method",
      "standardize_outcomes", "eps_sd",
      "ow", "tw", "row_w_base", "intercept",
      "lambda", "solver",
      "donor_screen", "build_Xy_for_unit",
      "demean_within_outcome_blocks", "fit_weights_fw"
    ),
    envir = environment()
  )

  # normalize NULLs / missing fields (fork crashes may yield NULL)
  out_list <- lapply(out_list, function(z) {
    if (is.null(z)) return(list(ok = FALSE, j = NA_integer_, msg = "Worker returned NULL (no result)."))
    if (!is.list(z)) return(list(ok = FALSE, j = NA_integer_, msg = "Worker returned non-list."))
    if (is.null(z$ok) || length(z$ok) == 0) z$ok <- FALSE
    if (is.null(z$msg) || length(z$msg) == 0) z$msg <- NA_character_
    z$msg <- as.character(z$msg)[1]
    z
  })

  ok_jobs <- vapply(out_list, function(z) isTRUE(z$ok), logical(1))
  if (!all(ok_jobs)) {
    bad <- which(!ok_jobs)
    msgs <- unique(vapply(out_list[bad], function(z) z$msg, character(1)))
    msgs <- msgs[!is.na(msgs)]
    stop(
      sprintf("Weight fitting failed for %d/%d treated units. Example errors: %s",
              length(bad), J, paste(head(msgs, 5), collapse = " | "))
    )
  }

  # unpack
  for (jj in seq_len(J)) {
    weights[[jj]] <- out_list[[jj]]$w
    donors_list[[jj]] <- out_list[[jj]]$donors
    donor_ids_list[[jj]] <- unit_ids[out_list[[jj]]$donors]
  }

  # -----------------------------
  # PASS 2: Optional pooled intercept adjustment (per outcome)
  # -----------------------------
  pooled_beta0 <- NULL
  if (isTRUE(pooled_adjustment)) {
    if (verbose) message("PASS 2: fitting pooled intercept adjustment per outcome.")
    pooled_beta0 <- numeric(M)

    for (m in seq_len(M)) {
      y_pre_all <- numeric(0)
      ysc_pre_all <- numeric(0)
      Y <- Y_list[[m]]

      for (jj in seq_len(J)) {
        i <- treated_units[jj]
        Ti <- as.integer(treat_time[i])
        pre_times <- (Ti - L):(Ti - 1)
        if (min(pre_times) < 1L) stop("Not enough pre periods for pooled adjustment given L.")

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
  # PASS 3: Compute unit-level effects
  # -----------------------------
  if (verbose) message("PASS 3: computing unit-level event-time effects.")
  for (jj in seq_len(J)) {
    j <- treated_units[jj]
    tau[jj, , ] <- unit_effects(
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
  }

  list(
    weights = weights,
    donors = donors_list,
    donor_ids = donor_ids_list,
    tau = tau,
    treated_units = treated_units,
    treated_unit_ids = treated_unit_ids,
    unit_ids = unit_ids,
    pooled_beta0 = pooled_beta0,
    pooled_adjustment = pooled_adjustment,
    nu = nu
  )
}
