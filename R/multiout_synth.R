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
#' @param donors_list Optional list of candidate donor sets, one per treated unit: each element an
#'   integer vector of unit indices (integers in 1..N) or a character vector of unit ids (matched against
#'   \code{unit_ids}; pass ids as character to avoid ambiguity with indices). Named by
#'   treated unit id (as character) or, if unnamed, aligned with the treated units in the order of
#'   \code{treated_units}. For a unit with an entry, screening runs inside that set (the eligibility
#'   rule \code{treat_time > treat_time[j] + K} still applies); units without an entry use the full
#'   eligible pool as before. \code{NULL} (default) reproduces the previous behaviour exactly.
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
#' @param augment `"none"` (default) or `"ridge"`: ridge-augmented synthetic control
#'   (Ben-Michael, Feller and Rothstein, Lemma 1). The donors' pre-treatment features
#'   used to fit the weights (all outcomes stacked over the L pre periods, standardised and
#'   de-meaned as for the weights) are centred across donors, and the weights are
#'   corrected by `A_c' (A_c A_c' + lambda I)^{-1} (y - A w)`. The corrected weights sum
#'   to one, can be negative (extrapolation outside the donors' convex hull), and are the
#'   same for all outcomes. They replace `weights` (the simplex weights are kept in
#'   `augment$weights_scm`), so `tau`, [tau_inference()] and [scm_gap_panel()] use them.
#' @param augment_lambda Ridge penalty relative to the mean eigenvalue of `A_c A_c'`
#'   (`Inf` = no augmentation). `NULL` (default) chooses it with [tune_augment()].
#'   A vector returns the weights for every value in `augment$path` (used for tuning)
#'   and leaves `weights` unaugmented.
#' @param augment_tune List of further arguments for [tune_augment()] (e.g. `grid`,
#'   `h`, `tune_units`, `one_se`).
#' @param fw_max_iter,fw_tol Solver controls passed to [fit_weights_fw()]: maximum
#'   iterations and tolerance on the relative Frank–Wolfe duality gap. Per-unit
#'   convergence is reported in `fw_diagnostics`.
#' @param standardize_outcomes Logical; standardize each outcome block in X/y for weight fitting only.
#' @param outcome_weights Optional length M weights on outcome blocks in the objective.
#' @param time_weights Optional length L weights on pre-periods (shared across outcomes).
#' @param intercept One of "none","global","outcome" controlling intercept projection in objective.
#' @param eps_sd Small jitter for standardization stability.
#' @param parallel_backend One of "none","fork","psock". Recommended on macOS: "psock".
#' @param n_cores Integer; number of cores/workers when parallel_backend != "none".
#'
#' @return A list with weights, donors, donor_ids, tau, treated_units, treated_unit_ids,
#'   unit_ids, pooled_beta0, pooled_adjustment, nu, L, K, fw_diagnostics, augment
#'   (`NULL` unless `augment = "ridge"`: method, lambda, tuning, weights_scm, and per-unit
#'   diagnostics `pre_rmse_scm`, `pre_rmse_aug`, `neg_mass`, `corr_l1`). `L` and `K` let [tau_inference()]
#'   and the panel helpers default to the window used in fitting.
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
    n_cores = max(1L, parallel::detectCores() - 1L),
    donors_list = NULL,
    fw_max_iter = 50000L,
    fw_tol = 1e-8,
    augment = c("none", "ridge"),
    augment_lambda = NULL,
    augment_tune = list()
) {
  # -----------------------------
  # Parse args / basic checks
  # -----------------------------
  screen_method <- match.arg(screen_method)
  solver <- match.arg(solver)
  intercept <- match.arg(intercept)
  parallel_backend <- match.arg(parallel_backend)
  augment <- match.arg(augment)

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
  # Optional per-unit candidate donor sets (resolved to unit indices, keyed by treated unit index)
  # -----------------------------
  cand_by_unit <- resolve_donors_list(donors_list, treated_units, unit_ids, N)

  # Ridge augmentation: choose the penalty by in-time placebo forecasts if not given
  augment_tuning <- NULL
  if (augment == "ridge") {
    if (is.null(augment_lambda)) {
      if (verbose) message("Tuning the ridge augmentation penalty (tune_augment).")
      augment_tuning <- do.call(tune_augment, c(list(
        Y_list = Y_list, treat_time = treat_time, L = L, K = K, treated_units = treated_units,
        max_donors = max_donors, screen_outcome = screen_outcome, screen_method = screen_method,
        lambda = lambda, standardize_outcomes = standardize_outcomes,
        outcome_weights = outcome_weights, time_weights = time_weights, intercept = intercept,
        eps_sd = eps_sd, parallel_backend = parallel_backend, n_cores = n_cores,
        fw_max_iter = fw_max_iter, fw_tol = fw_tol), augment_tune))
      augment_lambda <- augment_tuning$lambda
    }
    if (!is.numeric(augment_lambda) || !length(augment_lambda) || anyNA(augment_lambda) ||
        any(augment_lambda < 0)) {
      stop("augment_lambda must be nonnegative numeric (Inf = no augmentation).")
    }
  }

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
        method = screen_method,
        candidates = cand_by_unit[[as.character(j)]]
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
        fw = fit_weights_fw(Xw, yw, lambda = lambda, max_iter = fw_max_iter, tol = fw_tol,
                            verbose = FALSE)
      )
      a_or_na <- function(nm) { v <- attr(w, nm); if (is.null(v)) NA_real_ else as.numeric(v) }
      fw_diag <- c(rel_gap = a_or_na("fw_rel_gap"), iter = a_or_na("fw_iter"),
                   converged = as.numeric(isTRUE(attr(w, "fw_converged"))))
      w <- as.numeric(w)

      aug <- NULL
      if (augment == "ridge") {
        G <- .ridge_augment(Xw, yw, w, augment_lambda)
        if (length(augment_lambda) == 1L) {
          aug <- list(diag = .augment_diag(Xw, yw, w, G[, 1]), w_scm = w)
          w <- as.numeric(G[, 1])
        } else {
          aug <- list(path = G, w_scm = w)
        }
      }

      # normalize possible solver return types
      if (is.list(w)) {
        if (!is.null(w$w)) w <- w$w
        else if (!is.null(w$weights)) w <- w$weights
        else stop("Solver returned list without $w/$weights.")
      }
      w <- as.numeric(w)

      if (length(w) != length(donors)) stop("Weight length mismatch with donors.")
      if (any(!is.finite(w))) stop("Non-finite weights returned by solver.")

      list(ok = TRUE, j = j, donors = donors, w = w, fw_diag = fw_diag, aug = aug)
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
      "lambda", "solver", "cand_by_unit", "fw_max_iter", "fw_tol",
      "augment", "augment_lambda", ".ridge_augment", ".augment_diag",
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

  fw_diagnostics <- data.frame(
    treated_unit = treated_units,
    rel_gap = vapply(out_list, function(z) z$fw_diag[["rel_gap"]], numeric(1)),
    iter = vapply(out_list, function(z) z$fw_diag[["iter"]], numeric(1)),
    converged = vapply(out_list, function(z) z$fw_diag[["converged"]] == 1, logical(1))
  )
  if (!all(fw_diagnostics$converged)) {
    warning(sprintf(paste0("fit_weights_fw() did not reach the duality-gap tolerance for %d of %d ",
                           "treated units (max relative gap %.2g); see fit$fw_diagnostics and fw_max_iter."),
                    sum(!fw_diagnostics$converged), J, max(fw_diagnostics$rel_gap)), call. = FALSE)
  }

  # unpack
  for (jj in seq_len(J)) {
    weights[[jj]] <- out_list[[jj]]$w
    donors_list[[jj]] <- out_list[[jj]]$donors
    donor_ids_list[[jj]] <- unit_ids[out_list[[jj]]$donors]
  }

  augment_out <- NULL
  if (augment == "ridge") {
    augment_out <- list(method = "ridge", lambda = augment_lambda, tuning = augment_tuning,
                        weights_scm = stats::setNames(lapply(out_list, function(z) z$aug$w_scm), nm))
    if (length(augment_lambda) == 1L) {
      dg <- do.call(rbind, lapply(out_list, function(z) z$aug$diag))
      augment_out$diagnostics <- data.frame(treated_unit = treated_units, dg, row.names = NULL)
    } else {
      augment_out$path <- lapply(out_list, function(z) z$aug$path)
    }
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
    nu = nu,
    L = L,
    K = K,
    fw_diagnostics = fw_diagnostics,
    augment = augment_out
  )
}

#' Resolve a user-supplied list of candidate donor sets to unit indices
#'
#' @param donors_list As in [multiout_synth()] (may be NULL).
#' @param treated_units Integer indices of the treated units being estimated.
#' @param unit_ids Vector of unit identifiers, length N.
#' @param N Number of units.
#' @return A list keyed by treated unit index (as character) with integer donor indices, or an
#'   empty list when \code{donors_list} is NULL. Entries that resolve to no valid unit are dropped
#'   (the unit then uses the full eligible pool).
#' @keywords internal
resolve_donors_list <- function(donors_list, treated_units, unit_ids, N) {
  if (is.null(donors_list)) return(list())
  if (!is.list(donors_list)) stop("donors_list must be a list (one element per treated unit) or NULL.")
  ids_chr <- as.character(unit_ids)
  keys <- names(donors_list)
  if (is.null(keys) || any(!nzchar(keys))) {
    if (length(donors_list) != length(treated_units))
      stop("An unnamed donors_list must have one element per treated unit (length(treated_units)).")
    keys <- as.character(unit_ids[treated_units])
  }
  out <- list()
  for (i in seq_along(donors_list)) {
    j <- match(keys[i], ids_chr)
    if (is.na(j)) next                                        # key is not a unit id: ignore silently
    d <- donors_list[[i]]
    if (is.null(d) || !length(d)) next
    # character entries are unit ids; integer-valued numbers inside 1..N are unit indices; other numbers are ids
    idx <- if (is.character(d) || is.factor(d)) match(as.character(d), ids_chr)
           else if (is.numeric(d) && all(d == round(d), na.rm = TRUE) && all(d >= 1 & d <= N, na.rm = TRUE)) as.integer(d)
           else match(as.character(d), ids_chr)
    idx <- unique(idx[is.finite(idx) & idx >= 1L & idx <= N])
    if (length(idx)) out[[as.character(j)]] <- idx
  }
  out
}
