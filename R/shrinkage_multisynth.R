#' Estimate shrinkage intensity nu (empirical Bayes) for post-averaged effects
#'
#' This estimates a per-outcome shrinkage intensity nu_m that pools unit-level
#' post-averaged treatment effects toward the treated mean. The estimate uses
#' within-unit dispersion across event times k=use_k to proxy measurement noise.
#'
#' @param tau_hat Array of estimated event-time effects, dimension J x M x (K+1),
#'   typically \code{fit$tau} from \code{multiout_synth()}.
#' @param K Integer. Maximum event time used in \code{tau_hat} (must match the third dim).
#' @param use_k Integer vector. Event times to use for the post-average (default 1:K).
#' @param min_signal_var Numeric. Floor for the estimated signal variance (default 0).
#'
#' @return List with \code{nu} (length M), \code{V_bar} (avg measurement variance),
#'   and \code{sigma2_theta} (estimated between-unit signal variance).
#'
#' @export
estimate_nu_eb <- function(tau_hat, K, use_k = 1:K, min_signal_var = 0) {
  stopifnot(length(dim(tau_hat)) == 3)
  J <- dim(tau_hat)[1]
  M <- dim(tau_hat)[2]
  stopifnot(dim(tau_hat)[3] >= (K + 1))

  idx <- use_k + 1
  Kuse <- length(use_k)

  nu <- numeric(M)
  V_bar <- numeric(M)
  sigma2_theta <- numeric(M)

  for (m in seq_len(M)) {
    Tk <- tau_hat[, m, idx, drop = FALSE]
    Tk <- matrix(Tk, nrow = J, ncol = Kuse)

    # unit-level post-average
    theta_hat <- rowMeans(Tk)

    # within-unit variance over k as noise proxy
    s2_within <- apply(Tk, 1, stats::var)

    # var(mean) approx s2_within / Kuse
    V_i <- s2_within / Kuse
    Vb <- mean(V_i, na.rm = TRUE)

    Vtheta <- stats::var(theta_hat, na.rm = TRUE)
    s2_theta <- max(min_signal_var, Vtheta - Vb)

    nu[m] <- if ((Vb + s2_theta) > 0) Vb / (Vb + s2_theta) else 1

    V_bar[m] <- Vb
    sigma2_theta[m] <- s2_theta
  }

  list(nu = nu, V_bar = V_bar, sigma2_theta = sigma2_theta)
}

#' Multisynth-like shrinkage of unit-level post-averaged effects (Option B)
#'
#' Implements an explicit pooling layer on top of unit-level augmented synthetic control.
#' For each outcome m, the unit-level post-average \eqn{\hat\theta_i^{(m)}} is shrunk
#' toward the treated mean \eqn{\bar\theta^{(m)}}:
#' \deqn{\tilde\theta_i^{(m)} = (1-\nu_m)\hat\theta_i^{(m)} + \nu_m \bar\theta^{(m)}.}
#'
#' By default \eqn{\nu_m} is estimated via empirical Bayes using event-time dispersion.
#'
#' @param tau_hat Array J x M x (K+1) of event-time effects (usually \code{fit$tau}).
#' @param K Integer. Maximum event time.
#' @param use_k Integer vector. Event times used to compute post-average (default 1:K).
#' @param nu Numeric vector length M, or NULL to estimate via \code{estimate_nu_eb()}.
#'
#' @return List with:
#'   \itemize{
#'     \item \code{theta_hat}: J x M matrix of unit post-averages over use_k
#'     \item \code{theta_bar}: length M vector (treated mean)
#'     \item \code{theta_shrunk}: J x M matrix of shrunk unit effects
#'     \item \code{nu}: length M shrinkage vector used
#'     \item \code{nu_details}: output of \code{estimate_nu_eb()} if nu was estimated
#'   }
#'
#' @export
multisynth_shrink <- function(tau_hat, K, use_k = 1:K, nu = NULL) {
  stopifnot(length(dim(tau_hat)) == 3)
  J <- dim(tau_hat)[1]
  M <- dim(tau_hat)[2]

  idx <- use_k + 1

  theta_hat <- apply(tau_hat[, , idx, drop = FALSE], c(1, 2), mean)
  theta_hat <- matrix(theta_hat, nrow = J, ncol = M)

  theta_bar <- colMeans(theta_hat, na.rm = TRUE)

  nu_details <- NULL
  if (is.null(nu)) {
    nu_details <- estimate_nu_eb(tau_hat, K = K, use_k = use_k)
    nu <- nu_details$nu
  }
  stopifnot(length(nu) == M)

  theta_shrunk <- theta_hat
  for (m in seq_len(M)) {
    theta_shrunk[, m] <- (1 - nu[m]) * theta_hat[, m] + nu[m] * theta_bar[m]
  }

  list(
    theta_hat = theta_hat,
    theta_bar = theta_bar,
    theta_shrunk = theta_shrunk,
    nu = nu,
    nu_details = nu_details
  )
}
