#' Compute unit-level treatment effects for multiple outcomes
#'
#' Produces event-time treatment effect estimates for a single treated unit
#' \code{j} across multiple outcomes, using previously estimated donor weights.
#'
#' Supports an optional pooled regression-adjustment intercept (per outcome) that
#' shifts synthetic predictions by a constant estimated from pooled pre-treatment
#' residuals across treated units (multisynth-like adjustment).
#'
#' @param Y_list A list of outcome matrices, each \code{N x T}.
#' @param treat_time Numeric vector length \code{N}, with first treated period for each unit
#'   (or \code{Inf} for never-treated).
#' @param j Integer index of the treated unit.
#' @param donors Integer vector of donor unit indices.
#' @param w Numeric vector of donor weights of length \code{length(donors)}.
#' @param K Integer. Maximum event time to compute (computes \code{k=0,...,K}).
#' @param use_intercept Logical. If \code{TRUE}, uses a weighted DiD-style
#'   intercept-shift adjustment by subtracting pre-period means for treated and synthetic.
#'   If \code{FALSE}, returns raw gaps \code{Y_treated - Y_synth}.
#' @param L Integer. Number of pre-treatment periods used to define the intercept shift
#'   when \code{use_intercept=TRUE}.
#' @param pooled_beta0 Optional numeric vector of length \code{M} (number of outcomes).
#'   If provided, shifts the synthetic prediction by \code{pooled_beta0[m]} for outcome \code{m}
#'   in both pre and post periods. This corresponds to an optional pooled regression-adjustment
#'   intercept (multisynth-like) applied after weights are estimated.
#'
#' @return A numeric matrix of dimension \code{M x (K+1)} where \code{M = length(Y_list)}.
#'   Rows correspond to outcomes; columns correspond to event time \code{k=0,...,K}.
#'
#' @examples
#' N <- 100; TT <- 80
#' Y1 <- matrix(rnorm(N*TT), N, TT)
#' Y2 <- matrix(rnorm(N*TT), N, TT)
#' treat_time <- rep(Inf, N); treat_time[1:20] <- 50
#' j <- 1; donors <- which(treat_time > 60)
#' w <- rep(1/length(donors), length(donors))
#'
#' # No pooled adjustment:
#' tau <- unit_effects(list(Y1, Y2), treat_time, j, donors, w, K = 5, L = 10)
#'
#' # With pooled adjustment (illustrative):
#' tau2 <- unit_effects(list(Y1, Y2), treat_time, j, donors, w, K = 5, L = 10,
#'                      pooled_beta0 = c(0.1, -0.2))
#'
#' @export
unit_effects <- function(Y_list, treat_time, j, donors, w, K = 10,
                         use_intercept = TRUE, L = 20,
                         pooled_beta0 = NULL) {

  M <- length(Y_list)
  if (M < 1) stop("Y_list must contain at least one outcome matrix.")

  # Basic dimension checks
  N <- nrow(Y_list[[1]])
  TT <- ncol(Y_list[[1]])
  for (m in seq_len(M)) {
    if (!is.matrix(Y_list[[m]])) stop("Each element of Y_list must be a matrix.")
    if (nrow(Y_list[[m]]) != N || ncol(Y_list[[m]]) != TT) {
      stop("All outcome matrices in Y_list must have identical dimensions N x T.")
    }
  }

  if (length(treat_time) != N) stop("treat_time must have length N.")
  if (length(donors) < 1) stop("donors must contain at least one donor unit index.")
  if (length(w) != length(donors)) stop("Length of w must equal length(donors).")
  if (any(!is.finite(w))) stop("Weights w must be finite.")
  if (abs(sum(w) - 1) > 1e-6) warning("Weights do not sum to 1 (within tolerance).")

  if (!is.null(pooled_beta0)) {
    if (!is.numeric(pooled_beta0)) stop("pooled_beta0 must be numeric or NULL.")
    if (length(pooled_beta0) != M) {
      stop("pooled_beta0 must have length equal to length(Y_list) (i.e., M).")
    }
    if (any(!is.finite(pooled_beta0))) stop("pooled_beta0 values must be finite.")
  }

  if (!is.numeric(j) || length(j) != 1) stop("j must be a single integer index.")
  j <- as.integer(j)
  if (j < 1 || j > N) stop("j is out of bounds for units (1..N).")

  # Treatment time for unit j (1-indexed period)
  Tj <- treat_time[j]
  if (!is.finite(Tj)) stop("Unit j must be treated (finite treat_time[j]).")
  Tj <- as.integer(Tj)

  # Post times: event times k=0..K map to periods Tj..Tj+K
  post_times <- Tj:(Tj + K)
  if (max(post_times) > TT) stop("Not enough post periods for requested K.")

  # Pre times for intercept shift
  if (use_intercept) {
    pre_times <- (Tj - L):(Tj - 1)
    if (min(pre_times) < 1) stop("Not enough pre periods for intercept shift given L.")
  } else {
    pre_times <- integer(0)
  }

  tau <- matrix(NA_real_, nrow = M, ncol = K + 1)

  # Precompute donor outcome slices per outcome are potentially large; keep simple and safe here.
  for (m in seq_len(M)) {
    Y <- Y_list[[m]]

    # Observed treated outcomes post
    yj_post <- Y[j, post_times]

    # Synthetic post: w'Y_donors,post
    ysc_post <- as.numeric(crossprod(w, Y[donors, post_times, drop = FALSE]))

    # Apply pooled intercept adjustment to synthetic prediction, if provided
    if (!is.null(pooled_beta0)) {
      ysc_post <- ysc_post + pooled_beta0[m]
    }

    if (use_intercept) {
      yj_pre_mean <- mean(Y[j, pre_times])

      ysc_pre <- as.numeric(crossprod(w, Y[donors, pre_times, drop = FALSE]))
      if (!is.null(pooled_beta0)) {
        ysc_pre <- ysc_pre + pooled_beta0[m]
      }
      ysc_pre_mean <- mean(ysc_pre)

      tau[m, ] <- (yj_post - yj_pre_mean) - (ysc_post - ysc_pre_mean)
    } else {
      tau[m, ] <- yj_post - ysc_post
    }
  }

  tau
}
