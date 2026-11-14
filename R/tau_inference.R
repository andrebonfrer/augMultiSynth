# Inference for unit-level synthetic-control effects
#
# The estimand is the average post-treatment effect for treated unit j and
# outcome m over a set of post event times P (default 1..K):
#
#   tau_jm = mean_{k in P} g_{j,m}(T_j + k)  -  mean_{k = -L..-1} g_{j,m}(T_j + k),
#   g_{j,m}(t) = Y_{j,m}(t) - sum_d w_{jd} Y_{d,m}(t)               ("gap series")
#
# which is exactly mean_{k in P} fit$tau[j, m, k + 1] (unit_effects() de-means
# treated and synthetic paths on the same L pre periods; a pooled intercept,
# if any, cancels).
#
# Three representations of the same point estimate:
#   (1) OLS of g on (1, post) over event times -L..-1 and P;
#   (2) weighted two-way FE regression on the pseudo-panel (treated unit plus
#       donors, event-time window -L..-1 and P), y ~ unit FE + event-time FE + D,
#       regression weight 1 for the treated unit and w_d for donor d;
#   (3) mean_{k in P} fit$tau.
#
# (2) gives the right point estimate but the wrong sampling variance: donor
# rows enter the likelihood as if they were independent observations. In the
# synthetic-control model the donor paths and weights are conditioned on, and
# the treated outcome is Y_j(t) = mu_j + sum_d w_d Y_d(t) + tau_j D(t) + e(t).
# The likelihood of the pseudo-panel therefore factorises as
#   p(Y_j | Y_D; tau, mu, sigma^2) * p(Y_D),
# with no parameter of interest in p(Y_D). The first factor is the gap-series
# regression (1). Equivalently: apply the within-period contrast c = (1, -w)
# to the pseudo-panel. Because sum(c) = 0 it removes the event-time effects,
# it maps the unit effects to a single intercept mu_j = alpha_j - sum_d w_d alpha_d,
# and because donors are untreated in the window (donor_screen() eligibility:
# treat_time > T_j + K) it maps the treatment indicator to D(t). Inference on
# the contrasted data is inference on the gap series, so the two agree exactly.

#' Gap-series panel for unit-level inference
#'
#' Builds, for every treated unit and outcome, the series
#' \eqn{g(t) = Y_j(t) - \sum_d w_{jd} Y_d(t)} over event times \eqn{-L, \dots, -1}
#' and the post event times in `post`. This is the synthetic-control pseudo-panel
#' after the within-period contrast \eqn{(1, -w)}, i.e. the data on which the
#' likelihood for \eqn{\tau_j} is defined once donor paths are conditioned on
#' (see [tau_inference()]). A regression of `gap` on `post` within each
#' `treated_id` x `outcome` reproduces `mean(fit$tau[j, m, post + 1])` exactly.
#'
#' @param fit Output of [multiout_synth()].
#' @param Y_list List of N x T outcome matrices used to fit `fit`.
#' @param treat_time Length-N adoption times used to fit `fit`.
#' @param L Number of pre periods. Defaults to `fit$L` (stored by
#'   [multiout_synth()] from version 0.3.7).
#' @param post Integer vector of post event times (>= 0) that define the
#'   estimand. Default `1:K` with `K = dim(fit$tau)[3] - 1`, matching
#'   [eval_avg_1K()].
#' @param outcomes Outcome indices to include. Default all.
#' @return A data.frame with one row per treated unit x outcome x event time:
#'   `treated_id`, `treated_unit` (row index), `outcome`, `event_time`,
#'   `period`, `y_treated`, `y_synth`, `gap`, `post` (0/1).
#' @seealso [tau_inference()]
#' @export
scm_gap_panel <- function(fit, Y_list, treat_time, L = fit$L, post = NULL,
                          outcomes = NULL) {
  chk <- check_tau_inputs(fit, Y_list, treat_time, L, post, outcomes)
  L <- chk$L; post <- chk$post; outcomes <- chk$outcomes
  ev <- c(-L:-1, post)
  n_ev <- length(ev)

  out <- vector("list", length(fit$treated_units) * length(outcomes))
  r <- 0L
  for (jj in seq_along(fit$treated_units)) {
    j <- fit$treated_units[jj]
    per <- as.integer(treat_time[j]) + ev
    d <- fit$donors[[jj]]
    w <- fit$weights[[jj]]
    for (m in outcomes) {
      Y <- Y_list[[m]]
      y_tr <- as.numeric(Y[j, per])
      y_sc <- as.numeric(crossprod(w, Y[d, per, drop = FALSE]))
      r <- r + 1L
      out[[r]] <- data.frame(
        treated_id   = rep(fit$treated_unit_ids[jj], n_ev),
        treated_unit = rep(j, n_ev),
        outcome      = rep(m, n_ev),
        event_time   = ev,
        period       = per,
        y_treated    = y_tr,
        y_synth      = y_sc,
        gap          = y_tr - y_sc,
        post         = as.integer(ev >= 0),
        stringsAsFactors = FALSE
      )
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "L") <- L
  attr(res, "post") <- post
  res
}

#' Inference for unit-level synthetic-control effects
#'
#' For each treated unit \eqn{j} and outcome \eqn{m}, estimates
#' \deqn{\tau_{jm} = \bar g_{post} - \bar g_{pre},}
#' the post-period mean of the gap series minus its mean over the \eqn{L}
#' pre periods. The estimate equals `mean(fit$tau[j, m, post + 1])` exactly.
#'
#' Inference conditions on the donor paths and weights, as in the
#' synthetic-control model \eqn{Y_j(t) = \mu_j + \sum_d w_d Y_d(t) + \tau_j D(t) + e(t)}.
#' Under that model the pseudo-panel likelihood factorises into
#' \eqn{p(Y_j \mid Y_D; \tau, \mu, \sigma^2)\, p(Y_D)}, and the first factor is the
#' regression of the gap series on an intercept and a post indicator. Standard
#' errors are therefore
#' \deqn{se(\hat\tau) = \hat\sigma \sqrt{1/L + 1/|P|},\quad
#'       \hat\sigma^2 = RSS / (L + |P| - 2),}
#' with \eqn{L + |P| - 2} degrees of freedom. A weighted regression on the
#' pseudo-panel (treated unit plus donors, weight 1 and \eqn{w_d}) reproduces the
#' point estimate when it includes unit and event-time fixed effects, but not this
#' variance, because donor rows enter as observations.
#'
#' The variance treats the gap residuals as serially uncorrelated and does not
#' account for the pre-period fit of the weights, which makes pre-period
#' residuals smaller than out-of-sample ones. Both make the interval narrower
#' than a design-based (placebo or conformal) interval would be.
#'
#' For Bayesian inference under the same likelihood, use
#' `scmBayesPost::prepare_data_general(likelihood = "conditional", L = L, post = post)`
#' with the weights from `scmBayesPost::build_W_from_augMultiSynth()`, or pass
#' [scm_gap_panel()] to any sampler with an identity weight matrix. See
#' `scripts/sim_tau_bayes_vs_gap.R`.
#'
#' @inheritParams scm_gap_panel
#' @param level Confidence level for the interval.
#' @param check If `TRUE`, verify that the estimates equal the average of
#'   `fit$tau` over `post` (to 1e-8) and warn if not.
#' @return A data.frame with one row per treated unit x outcome: `treated_id`,
#'   `treated_unit`, `outcome`, `tau`, `se`, `df`, `t`, `p_value`, `lo`, `hi`,
#'   `mu` (pre-period mean gap), `sigma`, `n_pre`, `n_post`.
#' @seealso [scm_gap_panel()]
#' @examples
#' set.seed(1)
#' sim <- simulate_multi_outcome_panel(N = 200, T = 60, M = 2, K = 5, L = 15)
#' fit <- multiout_synth(sim$Yobs, sim$treat_time, L = 15, K = 5,
#'                       max_donors = 50, lambda = 0.1)
#' ti <- tau_inference(fit, sim$Yobs, sim$treat_time)
#' head(ti)
#' @export
tau_inference <- function(fit, Y_list, treat_time, L = fit$L, post = NULL,
                          outcomes = NULL, level = 0.95, check = TRUE) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1).")
  }
  chk <- check_tau_inputs(fit, Y_list, treat_time, L, post, outcomes)
  L <- chk$L; post <- chk$post; outcomes <- chk$outcomes
  ev <- c(-L:-1, post)
  n_post <- length(post)
  n <- L + n_post
  df <- n - 2L
  if (df < 1L) stop("Need L + length(post) >= 3 for a residual variance.")
  is_post <- ev >= 0
  J <- length(fit$treated_units)

  rows <- vector("list", length(outcomes))
  for (mi in seq_along(outcomes)) {
    m <- outcomes[mi]
    Y <- Y_list[[m]]
    G <- matrix(NA_real_, nrow = J, ncol = n)   # J x event times
    for (jj in seq_len(J)) {
      j <- fit$treated_units[jj]
      per <- as.integer(treat_time[j]) + ev
      G[jj, ] <- Y[j, per] -
        as.numeric(crossprod(fit$weights[[jj]], Y[fit$donors[[jj]], per, drop = FALSE]))
    }
    mu <- rowMeans(G[, !is_post, drop = FALSE])
    mp <- rowMeans(G[, is_post, drop = FALSE])
    tau <- mp - mu
    fitted <- matrix(mu, J, n) + outer(tau, as.numeric(is_post))
    rss <- rowSums((G - fitted)^2)
    sigma <- sqrt(rss / df)
    se <- sigma * sqrt(1 / L + 1 / n_post)
    tstat <- tau / se
    q <- stats::qt(1 - (1 - level) / 2, df)
    rows[[mi]] <- data.frame(
      treated_id   = fit$treated_unit_ids,
      treated_unit = fit$treated_units,
      outcome      = m,
      tau          = tau,
      se           = se,
      df           = df,
      t            = tstat,
      p_value      = 2 * stats::pt(-abs(tstat), df),
      lo           = tau - q * se,
      hi           = tau + q * se,
      mu           = mu,
      sigma        = sigma,
      n_pre        = L,
      n_post       = n_post,
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL

  if (isTRUE(check) && !is.null(fit$tau) && length(dim(fit$tau)) == 3L &&
      max(post) <= dim(fit$tau)[3] - 1L) {
    ref <- unlist(lapply(outcomes, function(m) {
      rowMeans(matrix(fit$tau[, m, post + 1L, drop = FALSE], nrow = J))
    }))
    dev <- max(abs(res$tau - ref))
    if (!is.finite(dev) || dev > 1e-8) {
      warning(sprintf(paste0(
        "tau_inference(): estimates differ from the average of fit$tau by up to %.3g. ",
        "Check that L, Y_list and treat_time are those used to fit `fit`."), dev))
    }
  }
  res
}

# Shared argument checks for the inference helpers
check_tau_inputs <- function(fit, Y_list, treat_time, L, post, outcomes) {
  if (!is.list(fit) || is.null(fit$weights) || is.null(fit$donors) || is.null(fit$treated_units)) {
    stop("fit must be the output of multiout_synth().")
  }
  if (!is.list(Y_list) || !length(Y_list)) stop("Y_list must be a non-empty list of matrices.")
  N <- nrow(Y_list[[1]]); TT <- ncol(Y_list[[1]])
  if (length(treat_time) != N) stop("treat_time must have length nrow(Y_list[[1]]).")
  if (is.null(L)) {
    stop("L is not stored in this fit (fits from versions before 0.3.7); pass L explicitly.")
  }
  L <- as.integer(L)
  if (length(L) != 1L || L < 1L) stop("L must be a positive integer.")
  Kfit <- if (!is.null(fit$tau)) dim(fit$tau)[3] - 1L else NA_integer_
  if (is.null(post)) {
    if (!is.finite(Kfit) || Kfit < 1L) stop("Cannot infer post from fit$tau; pass post explicitly.")
    post <- seq_len(Kfit)
  }
  post <- sort(unique(as.integer(post)))
  if (!length(post) || any(post < 0L)) stop("post must contain event times >= 0.")
  if (is.finite(Kfit) && max(post) > Kfit) {
    warning("post extends beyond the K used to fit `fit`; donors are only guaranteed ",
            "untreated through T_j + K.")
  }
  if (is.null(outcomes)) outcomes <- seq_along(Y_list)
  outcomes <- as.integer(outcomes)
  if (any(outcomes < 1L | outcomes > length(Y_list))) stop("outcomes out of range.")
  Tj <- as.integer(treat_time[fit$treated_units])
  if (any(Tj - L < 1L) || any(Tj + max(post) > TT)) {
    stop("Event window -L..max(post) falls outside the panel for some treated units.")
  }
  list(L = L, post = post, outcomes = outcomes)
}
