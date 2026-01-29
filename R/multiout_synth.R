#' Multi-outcome synthetic control with unit-level treatment effects
#'
#' Estimates donor weights and unit-level event-time treatment effects for multiple
#' outcomes in staggered adoption settings, using a stacked pre-treatment objective
#' across outcomes. Designed for scalability via donor screening and simplex solvers
#' (Frank--Wolfe). Optionally includes a pooled intercept adjustment (multisynth-like)
#' estimated from pooled pre-treatment residuals across treated units.
#'
#' @param Y_list List of outcome matrices, each of dimension \code{N x T}. All outcomes
#'   must share the same \code{N} and \code{T}. Units are in rows; time periods are in columns.
#' @param treat_time Numeric vector of length \code{N}. \code{treat_time[i]} is the first treated
#'   period for unit \code{i}. Use \code{Inf} for never-treated units.
#' @param treated_units Optional integer vector of treated unit indices to estimate. Defaults
#'   to all units with finite \code{treat_time}.
#' @param L Integer. Number of pre-treatment periods used to fit donor weights. For treated
#'   unit \code{j}, pre-periods are \code{(treat_time[j]-L):(treat_time[j]-1)}.
#' @param K Integer. Maximum event time. Returns unit-level effects for \code{k = 0,...,K}.
#' @param max_donors Integer. Maximum donor set size after screening. Defaults to \code{1000}.
#' @param screen_outcome Integer. Index of the outcome in \code{Y_list} used for donor screening.
#'   Defaults to \code{1}.
#' @param screen_method Donor screening method. One of \code{"cor"} or \code{"mse"}.
#' @param lambda Nonnegative scalar ridge penalty applied to donor weights in the simplex solver.
#' @param solver Solver choice. Currently \code{"fw"} (Frank--Wolfe) for simplex-constrained
#'   ridge regression.
#' @param verbose Logical. If \code{TRUE}, prints progress.
#'
#' @param standardize_outcomes Logical. If \code{TRUE}, standardizes each outcome block in the
#'   stacked pre-treatment design using donor pre-treatment moments (to improve comparability
#'   across outcomes). Standardization affects the weight-fitting objective only; treatment
#'   effects are computed on the original outcome scale.
#' @param outcome_weights Optional numeric vector of length \code{M = length(Y_list)} giving
#'   relative weights for each outcome block in the stacked pre-treatment objective. Larger values
#'   put more emphasis on matching that outcome in pre-treatment periods. Defaults to equal weights.
#' @param time_weights Optional numeric vector of length \code{L} giving relative weights for the
#'   \code{L} pre-treatment periods within each outcome block. Defaults to equal weights. A common
#'   use is to emphasize recent pre-treatment periods.
#'
#' @param intercept Character. Intercept handling in the weight-fitting objective. One of:
#'   \describe{
#'     \item{\code{"none"}}{No intercept projection; weights attempt to match both levels and trends.}
#'     \item{\code{"global"}}{Project out a single intercept shared across all outcomes and pre-periods
#'       in the stacked objective (equivalent to minimizing over \eqn{b} in \eqn{\|y - Xw - b\mathbf{1}\|_W^2}).}
#'     \item{\code{"outcome"}}{Project out a separate intercept within each outcome block (recommended for
#'       multi-outcome settings). This is equivalent to minimizing over \eqn{b_m} for each outcome block
#'       in \eqn{\sum_m \|y^{(m)} - X^{(m)}w - b_m\mathbf{1}\|_{W_m}^2}.}
#'   }
#'
#' @param pooled_adjustment Logical. If \code{TRUE}, estimates a pooled intercept correction per outcome
#'   from pooled pre-treatment residuals across treated units (multisynth-like). This produces
#'   \code{pooled_beta0} which is applied when computing unit-level effects.
#' @param nu Nonnegative ridge penalty for the pooled intercept correction (per outcome). Larger \code{nu}
#'   shrinks the pooled intercept toward zero.
#'
#' @param eps_sd Small positive constant used as numerical jitter in standardization to avoid division by
#'   near-zero standard deviations (only relevant when \code{standardize_outcomes = TRUE}).
#'
#' @details
#' For each treated unit \code{j}, the estimator:
#' \enumerate{
#'   \item Constructs a donor set \code{donors_j} by screening eligible donors and restricting to at most
#'         \code{max_donors}.
#'   \item Builds a stacked pre-treatment design \code{X_j} and target \code{y_j} by stacking the \code{L}
#'         pre-treatment periods across the \code{M} outcomes.
#'   \item Optionally projects out intercepts in the stacked objective via \code{intercept}.
#'   \item Fits simplex-constrained weights \code{w_j} by minimizing a weighted ridge objective
#'         \eqn{(y_j - X_j w)^\top W (y_j - X_j w) + \lambda\|w\|_2^2}, implemented by row pre-whitening.
#'   \item Computes outcome-by-event-time treatment effects \code{tau[j,m,k]} using the fitted weights
#'         (and optional pooled adjustment).
#' }
#'
#' The returned unit-level effects can be exported to second-stage analyses (e.g., regressing
#' \code{tau[,m,1]} or average effects over \code{k=1..K} on unit covariates).
#'
#' @return A list with components:
#' \describe{
#'   \item{weights}{Named list of donor weight vectors, one per treated unit.}
#'   \item{donors}{Named list of donor index vectors, one per treated unit.}
#'   \item{tau}{Numeric array of dimension \code{J x M x (K+1)} containing unit-level event-time
#'     treatment effects, where \code{J = length(treated_units)} and \code{M = length(Y_list)}.}
#'   \item{treated_units}{Integer vector of treated units used in estimation.}
#'   \item{pooled_beta0}{NULL or numeric vector length \code{M} of pooled intercept adjustments.}
#'   \item{pooled_adjustment}{Logical flag as passed in.}
#'   \item{nu}{Penalty parameter used for pooled adjustment.}
#' }
#'
#' @examples
#' set.seed(1)
#' N <- 300; T <- 120
#' Y1 <- matrix(rnorm(N*T), N, T)
#' Y2 <- matrix(rnorm(N*T), N, T)
#' treat_time <- rep(Inf, N); treat_time[1:60] <- sample(60:90, 60, replace = TRUE)
#'
#' fit <- multiout_synth(
#'   Y_list = list(Y1, Y2),
#'   treat_time = treat_time,
#'   L = 20, K = 5,
#'   max_donors = 150,
#'   intercept = "outcome",
#'   outcome_weights = c(1, 1),
#'   time_weights = rep(1, 20),
#'   pooled_adjustment = FALSE,
#'   verbose = FALSE
#' )
#' dim(fit$tau)
#'
#' @export
multiout_synth <- function(Y_list, treat_time, treated_units = NULL,
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
                           outcome_weights = NULL,      # length M, default rep(1, M)
                           time_weights = NULL,         # length L, default rep(1, L)
                           intercept = c("none", "global", "outcome"),
                           eps_sd = 1e-8,
                           parallel = FALSE,
                           n_cores = max(1L, parallel::detectCores() - 1L),
                           backend = c("auto","fork","psock")
                           ) {

  screen_method <- match.arg(screen_method)
  solver <- match.arg(solver)
  intercept <- match.arg(intercept)

  # choose backend
  backend <- match.arg(backend)

  in_rstudio <- identical(Sys.getenv("RSTUDIO"), "1")
  if (backend == "auto") {
    if (.Platform$OS.type == "windows" || in_rstudio) backend <- "psock" else backend <- "fork"
  }
  if (backend == "fork" && .Platform$OS.type == "windows") {
    stop("backend='fork' is not supported on Windows. Use backend='psock'.")
  }

  # find dimension of outcome vector
  M <- length(Y_list)

  # coerce doubles on Y_list
  for (m in seq_len(M)) {
    Y_list[[m]] <- as.matrix(Y_list[[m]])
    storage.mode(Y_list[[m]]) <- "double"
  }

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

  # -----------------------------
  # Determine treated units and apply validity/feasibility filters
  # -----------------------------
  if (is.null(treated_units)) {
    treated_units <- which(is.finite(treat_time))
  } else {
    treated_units <- as.integer(treated_units)
  }

  # Keep only truly treated units (finite treat_time)
  treated_units <- treated_units[is.finite(treat_time[treated_units])]

  if (length(treated_units) < 1L) {
    stop("No treated units remain after filtering to finite treat_time.")
  }

  # Feasibility: must have L pre periods and K post periods available
  Ti <- as.integer(treat_time[treated_units])
  ok <- (Ti - L >= 1L) & (Ti + K <= TT)

  treated_units <- treated_units[ok]

  if (length(treated_units) < 1L) {
    stop("No treated units satisfy feasibility (enough pre/post periods) given L and K.")
  }

  # Now define J *after* filtering
  J <- length(treated_units)


  if (!is.numeric(L) || length(L) != 1 || L < 1) stop("L must be a positive integer.")
  if (!is.numeric(K) || length(K) != 1 || K < 0) stop("K must be a nonnegative integer.")
  L <- as.integer(L)
  K <- as.integer(K)

  # ---- weights for the stacked objective ----
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

  # check lambda and nu
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) stop("lambda must be >= 0.")
  if (!is.numeric(nu) || length(nu) != 1 || nu < 0) stop("nu must be >= 0.")

  tau <- array(NA_real_, dim = c(J, M, K + 1))
  weights <- vector("list", J)
  donors_list <- vector("list", J)
  names(weights) <- names(donors_list) <- as.character(treated_units)


  # -----------------------------
  # PASS 1: Fit weights per treated unit (parallelizable)
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

    # defensive coercions
    donors <- as.integer(donors)
    donors <- donors[is.finite(donors) & donors >= 1L & donors <= N]
    if (length(donors) < 2L) stop("Too few donors after screening.")

    XY <- build_Xy_for_unit(
      Y_list, treat_time, j, donors, L,
      standardize_outcomes = standardize_outcomes,
      eps_sd = eps_sd
    )
    X <- XY$X
    y <- XY$y

    ow <- outcome_weights
    tw <- time_weights
    if (is.null(ow)) ow <- rep(1, M)
    if (length(ow) != M) stop("outcome_weights must have length M.")
    if (is.null(tw)) tw <- rep(1, L)
    if (length(tw) != L) stop("time_weights must have length L.")

    row_w <- rep(ow, each = L) * rep(tw, times = M)
    if (any(!is.finite(row_w)) || any(row_w < 0) || sum(row_w) <= 0) {
      stop("row weights must be finite, nonnegative, and sum to > 0.")
    }

    if (intercept == "global") {
      sw <- sum(row_w)
      y <- y - sum(row_w * y) / sw
      Xbar <- as.numeric(crossprod(row_w, X) / sw)
      X <- X - matrix(Xbar, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
    } else if (intercept == "outcome") {
      tmp <- demean_within_outcome_blocks(X, y, L = L, time_weights = tw)
      X <- tmp$X
      y <- tmp$y
    }

    sqrtw <- sqrt(row_w)
    X <- X * sqrtw
    y <- y * sqrtw

    w <- switch(
      solver,
      fw = fit_weights_fw(X, y, lambda = lambda, max_iter = 1000, tol = 1e-5, verbose = FALSE)
    )

    # ensure numeric weights (covers case fit_weights_fw returns list)
    if (is.list(w)) {
      if (!is.null(w$w)) w <- w$w
      else if (!is.null(w$weights)) w <- w$weights
      else stop("fit_weights_fw() returned a list but no $w/$weights element.")
    }
    w <- as.numeric(w)

    if (length(w) != length(donors)) stop("Weight length mismatch with donors.")
    if (any(!is.finite(w))) stop("Non-finite weights returned by solver.")

    list(ok = TRUE, j = j, donors = donors, w = w)

  }, error = function(e) {
    list(ok = FALSE, j = j, msg = conditionMessage(e))
  })
}


  # run in parallel or sequential
  treated_units_chr <- as.character(treated_units)

  if (isTRUE(parallel)) {

    if (backend == "fork") {
      out_list <- parallel::mclapply(
        treated_units,
        function(j) tryCatch(fit_one_unit(j), error = function(e) list(ok = FALSE, j = j, msg = conditionMessage(e))),
        mc.cores = n_cores,
        mc.preschedule = FALSE
      )

    } else if (backend == "psock") {

      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export everything the worker needs
      parallel::clusterExport(
        cl,
        varlist = c(
          "fit_one_unit",
          "Y_list","treat_time","M","N","TT","L","K",
          "max_donors","screen_outcome","screen_method",
          "standardize_outcomes","eps_sd",
          "outcome_weights","time_weights","intercept",
          "lambda","solver",
          "donor_screen","build_Xy_for_unit",
          "demean_within_outcome_blocks","fit_weights_fw"
        ),
        envir = environment()
      )

      out_list <- parallel::parLapply(
        cl,
        treated_units,
        function(j) tryCatch(fit_one_unit(j), error = function(e) list(ok = FALSE, j = j, msg = conditionMessage(e)))
      )
    }

  } else {
    out_list <- lapply(
      treated_units,
      function(j) tryCatch(fit_one_unit(j), error = function(e) list(ok = FALSE, j = j, msg = conditionMessage(e)))
    )
  }

  # robustness check
  if (length(out_list) != J) {
    stop("Internal error: out_list length mismatch with J. This should not happen.")
  }


  # unpack results into weights/donors_list
  for (jj in seq_len(J)) {
    weights[[jj]] <- out_list[[jj]]$w
    donors_list[[jj]] <- out_list[[jj]]$donors
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
