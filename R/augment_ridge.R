# Ridge augmentation (augmented synthetic control; Ben-Michael, Feller and
# Rothstein, "The Augmented Synthetic Control Method", Lemma 1).
#
# For treated unit j with synthetic-control weights w over n donors, let A be
# the P x n matrix of donor pre-treatment features and y the treated unit's
# features, exactly as used to fit w (all M outcomes stacked over the L pre
# periods, standardised and de-meaned within outcome, pre-whitened by the row
# weights). Centre the features across donors, A_c = A - rowMeans(A). Ridge
# regression of donors' outcomes on A_c gives the bias correction
#   m_hat(x_j) - sum_i w_i m_hat(x_i) = gamma_corr' Y_D,
#   gamma_corr = A_c' (A_c A_c' + lambda I_P)^{-1} (y - A w),
# so the augmented estimator is the weighting estimator with weights
#   gamma = w + gamma_corr.
# Because the rows of A_c sum to zero across donors, sum(gamma_corr) = 0 and
# gamma sums to one. gamma does not depend on which outcome is predicted, so a
# single gamma serves all M outcomes. Weights can be negative (extrapolation),
# which happens when the treated unit lies outside the donors' convex hull.
#
# lambda is given relative to tr(A_c A_c') / P (the mean eigenvalue), so
# lambda_rel = 1 penalises at the average scale of the features; Inf gives back w.

# Augmented weights for one unit, for each lambda in lambda_rel (n x length(lambda_rel)).
.ridge_augment <- function(A, y, w, lambda_rel) {
  A <- as.matrix(A)
  P <- nrow(A); n <- ncol(A)
  Ac <- A - rowMeans(A)
  r <- as.numeric(y - A %*% w)
  scale <- sum(Ac^2) / P
  if (!is.finite(scale) || scale <= 0) scale <- 1
  out <- matrix(w, n, length(lambda_rel))
  fin <- which(is.finite(lambda_rel))
  if (!length(fin)) return(out)
  if (P <= n) {
    e <- eigen(tcrossprod(Ac), symmetric = TRUE)
    Ur <- as.numeric(crossprod(e$vectors, r))
    for (k in fin) {
      v <- e$vectors %*% (Ur / (pmax(e$values, 0) + lambda_rel[k] * scale))
      out[, k] <- w + as.numeric(crossprod(Ac, v))
    }
  } else {
    e <- eigen(crossprod(Ac), symmetric = TRUE)
    Vr <- as.numeric(crossprod(e$vectors, crossprod(Ac, r)))
    for (k in fin) {
      out[, k] <- w + as.numeric(e$vectors %*% (Vr / (pmax(e$values, 0) + lambda_rel[k] * scale)))
    }
  }
  out
}

# Per-unit diagnostics for augmented weights g (vector) given A, y, w.
.augment_diag <- function(A, y, w, g) {
  c(pre_rmse_scm = sqrt(mean((y - A %*% w)^2)),
    pre_rmse_aug = sqrt(mean((y - A %*% g)^2)),
    neg_mass = sum(-g[g < 0]),
    corr_l1 = sum(abs(g - w)))
}

#' Tune the ridge augmentation penalty by in-time placebo forecasts
#'
#' For a sample of treated units, moves adoption back `h` periods, fits
#' [multiout_synth()] with ridge augmentation on the remaining `L - h` pre periods
#' for every value in `grid` (and `Inf`, i.e. no augmentation), and forecasts the
#' `h` held-out pre periods, where the true effect is zero. The criterion is the
#' mean squared pseudo effect in outcome-sd units, pooled over units, outcomes
#' and held-out periods. With `one_se = TRUE` the largest penalty within one
#' standard error of the minimum is chosen, which favours less extrapolation.
#'
#' Ben-Michael, Feller and Rothstein choose the penalty by leave-one-period-out
#' cross-validation over the pre-period; holding out the last `h` periods instead
#' scores forecasts forward in time, as the post-period counterfactual is used.
#'
#' @param Y_list,treat_time As in [multiout_synth()].
#' @param L,K Pre and post window of the main fit.
#' @param grid Relative penalties to compare (see [multiout_synth()], `augment_lambda`).
#' @param h Held-out pre periods; default `min(K, floor(L / 3))`.
#' @param treated_units Candidate treated units (default: all with a feasible window).
#' @param tune_units Number of units sampled for tuning.
#' @param one_se Use the one-standard-error rule.
#' @param seed Seed for the unit sample.
#' @param ... Further arguments passed to [multiout_synth()] (e.g. `max_donors`,
#'   `lambda`, `intercept`, `standardize_outcomes`, `parallel_backend`). Use the
#'   same settings as the main fit.
#' @return A list with `lambda` (chosen relative penalty; `Inf` means augmentation
#'   did not help), `lambda_min`, `summary` (one row per penalty: `mse`, `se`,
#'   `rmse`, `neg_mass` (mean negative weight mass)), `units`, and `h`.
#' @seealso [multiout_synth()]
#' @export
tune_augment <- function(Y_list, treat_time, L, K, grid = 10^seq(-3, 2, by = 0.5),
                         h = NULL, treated_units = NULL, tune_units = 50L,
                         one_se = TRUE, seed = 1L, ...) {
  if (is.null(h)) h <- min(K, floor(L / 3))
  h <- as.integer(h)
  if (h < 1L || L - h < 3L) stop("Need h >= 1 and L - h >= 3.")
  grid <- sort(unique(c(as.numeric(grid), Inf)))
  TT <- ncol(Y_list[[1]])
  if (is.null(treated_units)) treated_units <- which(is.finite(treat_time))
  ok <- treated_units[treat_time[treated_units] - L >= 1 & treat_time[treated_units] + K <= TT]
  if (!length(ok)) stop("No treated units with a feasible window.")
  units <- if (length(ok) <= tune_units) ok else {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed); sort(sample(ok, tune_units))
  }

  tt_b <- treat_time
  tt_b[units] <- treat_time[units] - h
  fb <- multiout_synth(Y_list, tt_b, treated_units = units, L = L - h, K = h - 1L,
                       augment = "ridge", augment_lambda = grid, ...)

  M <- length(Y_list)
  err <- array(NA_real_, dim = c(length(units), length(grid)))
  neg <- err
  for (u in seq_along(units)) {
    j <- fb$treated_units[u]; d <- fb$donors[[u]]
    pre <- (tt_b[j] - (L - h)):(tt_b[j] - 1L); held <- tt_b[j] + 0:(h - 1L)
    G <- fb$augment$path[[u]]
    for (k in seq_along(grid)) {
      g <- G[, k]
      e <- vapply(seq_len(M), function(m) {
        Y <- Y_list[[m]]
        s <- stats::sd(as.numeric(Y[d, pre])); if (!is.finite(s) || s <= 0) s <- 1
        gap <- Y[j, c(pre, held)] - as.numeric(crossprod(g, Y[d, c(pre, held), drop = FALSE]))
        mean(((gap[length(pre) + seq_len(h)] - mean(gap[seq_along(pre)])) / s)^2)
      }, numeric(1))
      err[u, k] <- mean(e)
      neg[u, k] <- sum(-g[g < 0])
    }
  }
  mse <- colMeans(err)
  se <- apply(err, 2, stats::sd) / sqrt(nrow(err))
  summ <- data.frame(lambda = grid, mse = mse, se = se, rmse = sqrt(mse), neg_mass = colMeans(neg))
  i_min <- which.min(mse)
  lam <- if (isTRUE(one_se)) max(grid[mse <= mse[i_min] + se[i_min]]) else grid[i_min]
  list(lambda = lam, lambda_min = grid[i_min], summary = summ, units = units, h = h)
}
