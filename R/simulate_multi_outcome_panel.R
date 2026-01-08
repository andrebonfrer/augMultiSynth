#' Simulate multi-outcome staggered-adoption panel data under an interactive fixed effects DGP
#'
#' @param N Integer. Number of units.
#' @param T Integer. Number of time periods.
#' @param M Integer. Number of outcomes.
#' @param F Integer. Number of latent factors.
#' @param treated_frac Numeric in (0,1). Fraction of units treated.
#' @param L Integer. Required pre-treatment periods (feasibility).
#' @param K Integer. Required post-treatment periods (feasibility).
#' @param te_scale Numeric. Scale of treatment effects.
#' @param noise_sd Numeric. SD of idiosyncratic noise.
#' @param alpha_sd Numeric. SD of unit fixed effects.
#' @param delta_sd Numeric. SD of time fixed effects.
#' @param link_tau_to_loadings Logical. If TRUE, heterogeneous treatment effects depend on factor loadings.
#' @param dynamic_te Logical. If TRUE, treatment effects evolve over event time (simple ramp).
#' @param outcome_specific_intercepts Logical. If TRUE, simulate outcome-specific unit and time intercepts
#'   \eqn{\alpha_i^{(m)}} and \eqn{\delta_t^{(m)}}. If FALSE (default), use common \eqn{\alpha_i} and \eqn{\delta_t}
#'   shared across outcomes.
#'
#' @return List with Y0, Yobs, treat_time, treated, tau_unit, tau_true_event, and DGP components.
#'
#' @export
simulate_multi_outcome_panel <- function(
    N = 4000, T = 200, M = 3, F = 3,
    treated_frac = 0.4,
    L = 30, K = 10,
    te_scale = 0.5, noise_sd = 0.3,
    alpha_sd = 0.5, delta_sd = 0.3,
    link_tau_to_loadings = FALSE,
    dynamic_te = FALSE,
    outcome_specific_intercepts = FALSE
) {
  # --- Interactive fixed effects (factors) ---
  # f_t: T x F, lambda_i: N x F
  f_t <- matrix(rnorm(T * F), nrow = T, ncol = F)
  lambda_i <- matrix(rnorm(N * F), nrow = N, ncol = F)

  # Outcome-specific mapping of factors (shared factors, different loadings-to-outcome rotations)
  # B_m: F x F (diagonal scaling by default)
  Bm <- lapply(seq_len(M), function(m) diag(runif(F, 0.5, 1.5)))

  # --- Fixed effects / intercepts ---
  # Either common across outcomes or outcome-specific
  if (isTRUE(outcome_specific_intercepts)) {
    # alpha_im: N x M, delta_tm: T x M
    alpha_im <- matrix(rnorm(N * M, sd = alpha_sd), nrow = N, ncol = M)
    delta_tm <- matrix(rnorm(T * M, sd = delta_sd), nrow = T, ncol = M)
    alpha_i <- NULL
    delta_t <- NULL
  } else {
    alpha_i <- rnorm(N, sd = alpha_sd)  # unit FE shared across outcomes
    delta_t <- rnorm(T, sd = delta_sd)  # time FE shared across outcomes
    alpha_im <- NULL
    delta_tm <- NULL
  }

  # Untreated outcomes Y0^{(m)}: N x T
  # Y0_it^(m) = alpha + delta + lambda_i' (B_m f_t) + eps
  Y0 <- lapply(seq_len(M), function(m) {
    Ft_m <- f_t %*% Bm[[m]]                # T x F
    signal <- lambda_i %*% t(Ft_m)         # N x T

    if (isTRUE(outcome_specific_intercepts)) {
      mu <- matrix(alpha_im[, m], nrow = N, ncol = T) +
        matrix(delta_tm[, m], nrow = N, ncol = T, byrow = TRUE) +
        signal
    } else {
      mu <- matrix(alpha_i, nrow = N, ncol = T) +
        matrix(delta_t, nrow = N, ncol = T, byrow = TRUE) +
        signal
    }

    mu + matrix(rnorm(N * T, sd = noise_sd), nrow = N, ncol = T)
  })

  # --- Staggered adoption times ---
  treated <- sample.int(N, size = floor(treated_frac * N))
  treat_time <- rep(Inf, N)
  possible_T <- (L + 5):(T - K - 5)
  treat_time[treated] <- sample(possible_T, length(treated), replace = TRUE)

  # --- Heterogeneous treatment effects tau_i^{(m)} ---
  tau_unit <- matrix(0, nrow = N, ncol = M)

  if (isTRUE(link_tau_to_loadings)) {
    for (m in seq_len(M)) {
      tau_unit[, m] <- te_scale * (0.7 * scale(lambda_i[, 1])[, 1] + 0.3 * rnorm(N))
    }
  } else {
    tau_unit <- matrix(rnorm(N * M, sd = te_scale), nrow = N, ncol = M)
  }

  # Optional dynamic treatment effect path over event time (k = 0..T-1)
  ramp <- function(k) if (k < 5) (k + 1) / 5 else 1
  tau_path_shape <- sapply(0:(T - 1), ramp)  # length T

  # --- Apply treatment to build observed outcomes ---
  Yobs <- Y0
  for (m in seq_len(M)) {
    Y <- Y0[[m]]
    for (j in treated) {
      Tj <- treat_time[j]
      if (is.finite(Tj) && Tj <= T) {
        if (isTRUE(dynamic_te)) {
          ks <- 0:(T - Tj)
          Y[j, Tj:T] <- Y[j, Tj:T] + tau_unit[j, m] * tau_path_shape[ks + 1]
        } else {
          Y[j, Tj:T] <- Y[j, Tj:T] + tau_unit[j, m]
        }
      }
    }
    Yobs[[m]] <- Y
  }

  # tau_true_event: N x M x T array where [i,m,k+1] is true tau at event time k (k=0..T-1)
  tau_true_event <- array(0, dim = c(N, M, T))
  for (m in seq_len(M)) {
    for (i in seq_len(N)) {
      tau_true_event[i, m, ] <- tau_unit[i, m] * (if (isTRUE(dynamic_te)) tau_path_shape else rep(1, T))
    }
  }

  list(
    Y0 = Y0,
    Yobs = Yobs,
    treat_time = treat_time,
    treated = treated,
    tau_unit = tau_unit,
    tau_true_event = tau_true_event,
    tau_path_shape = tau_path_shape,  # helpful for evaluation; avoids reconstruction
    # DGP components (useful for diagnostics)
    alpha_i = alpha_i,
    delta_t = delta_t,
    alpha_im = alpha_im,
    delta_tm = delta_tm,
    lambda_i = lambda_i,
    f_t = f_t,
    Bm = Bm,
    link_tau_to_loadings = link_tau_to_loadings,
    dynamic_te = dynamic_te,
    outcome_specific_intercepts = outcome_specific_intercepts
  )
}
