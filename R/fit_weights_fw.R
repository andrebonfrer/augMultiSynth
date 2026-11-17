#' Fit synthetic control weights using Frank–Wolfe on the simplex
#'
#' Solves a simplex-constrained least squares problem
#' \deqn{\min_{w \in \Delta} \|y - Xw\|^2 + \lambda \|w\|^2}
#' where \code{w} lies on the probability simplex (\code{w >= 0}, \code{sum(w)=1}).
#'
#' Designed for large donor sets with warm-start potential and low memory
#' overhead, particularly when paired with donor screening.
#'
#' @param X Numeric matrix of dimension \code{P x Nd}. Pre-treatment donor design.
#' @param y Numeric vector of length \code{P}. Pre-treatment treated target.
#' @param lambda Nonnegative scalar ridge penalty on weights. Defaults to \code{1e-3}.
#' @param max_iter Maximum number of iterations. Defaults to \code{50000}.
#' @param tol Convergence tolerance on the Frank–Wolfe duality gap, relative to
#'   \code{max(1, f(w))}. Defaults to \code{1e-8}.
#' @param verbose Logical. If \code{TRUE}, prints occasional progress messages.
#' @param w0 Optional starting weights (projected to the simplex by normalising
#'   their nonnegative part). Defaults to uniform weights.
#'
#' @details
#' Pairwise Frank–Wolfe: each iteration moves weight from the active donor with
#' the largest gradient (the "away" vertex) to the donor with the smallest
#' gradient (the Frank–Wolfe vertex), with exact line search on the quadratic.
#' The gradient is updated incrementally, so an iteration costs O(Nd) after
#' forming \eqn{Q = X'X + \lambda I} once.
#'
#' Iteration stops when the Frank–Wolfe duality gap
#' \deqn{g(w) = \nabla f(w)'(w - e_s),\quad s = \arg\min_i \nabla f(w)_i,}
#' falls below \code{tol * max(1, f(w))}. Because \eqn{f} is convex,
#' \eqn{f(w) - f^* \le g(w)}, so the gap certifies how far the objective is from
#' its minimum. (Versions before 0.3.8 used plain Frank–Wolfe steps and stopped
#' when the objective changed little between iterations, which could stop well
#' short of the optimum.)
#'
#' @return Numeric vector of length \code{Nd} with simplex weights summing to 1,
#'   with attributes \code{fw_gap} (final duality gap), \code{fw_rel_gap}
#'   (gap relative to \code{max(1, f(w))}), \code{fw_iter}, and
#'   \code{fw_converged}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200*50), 200, 50)
#' w_true <- rep(0, 50); w_true[c(3,10,25)] <- c(0.2, 0.5, 0.3)
#' y <- X %*% w_true + rnorm(200, sd = 0.1)
#' w_hat <- fit_weights_fw(X, as.numeric(y), lambda = 1e-4)
#' sum(w_hat)  # 1
#' attr(w_hat, "fw_rel_gap")
#'
#' @export
fit_weights_fw <- function(X, y, lambda = 1e-3,
                           max_iter = 50000, tol = 1e-8, verbose = FALSE, w0 = NULL) {

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(y)) stop("y must be numeric.")
  y <- as.numeric(y)
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y).")
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0) stop("lambda must be nonnegative.")

  Nd <- ncol(X)
  if (Nd < 2) stop("Need at least 2 donors (ncol(X) >= 2).")

  Q <- crossprod(X)
  diag(Q) <- diag(Q) + lambda
  b <- as.numeric(crossprod(X, y))
  yy <- sum(y^2)

  if (is.null(w0)) {
    w <- rep(1 / Nd, Nd)
  } else {
    if (length(w0) != Nd) stop("w0 must have length ncol(X).")
    w <- pmax(as.numeric(w0), 0)
    if (sum(w) <= 0) w <- rep(1 / Nd, Nd)
    w <- w / sum(w)
  }

  Qw <- as.numeric(Q %*% w)
  g <- 2 * (Qw - b)                       # gradient
  f <- sum(w * Qw) - 2 * sum(b * w) + yy  # objective
  gap <- Inf; it <- 0L; converged <- FALSE

  for (it in seq_len(max_iter)) {
    s <- which.min(g)
    gap <- sum(g * w) - g[s]
    if (gap <= tol * max(1, f)) { converged <- TRUE; break }

    act <- which(w > 0)
    v <- act[which.max(g[act])]
    if (v == s) { converged <- TRUE; break }

    # pairwise direction d = e_s - e_v, exact line search, step <= w_v
    dQd <- Q[s, s] + Q[v, v] - 2 * Q[s, v]
    slope <- g[s] - g[v]                   # = g'd < 0
    gam <- if (dQd > 0) -slope / (2 * dQd) else w[v]
    gam <- min(max(gam, 0), w[v])
    if (gam <= 0) { converged <- TRUE; break }

    w[s] <- w[s] + gam
    w[v] <- w[v] - gam
    if (w[v] < 1e-15) w[v] <- 0
    Qd <- Q[, s] - Q[, v]
    f <- f + gam * slope + gam^2 * dQd
    g <- g + 2 * gam * Qd
    if (it %% 2000L == 0L) {               # refresh to avoid accumulated rounding
      Qw <- as.numeric(Q %*% w)
      g <- 2 * (Qw - b)
      f <- sum(w * Qw) - 2 * sum(b * w) + yy
    }

    if (verbose && it %% 1000 == 0) {
      message(sprintf("iter %d: obj=%.8g, gap=%.3g", it, f, gap))
    }
  }

  w <- w / sum(w)
  attr(w, "fw_gap") <- gap
  attr(w, "fw_rel_gap") <- gap / max(1, f)
  attr(w, "fw_iter") <- it
  attr(w, "fw_converged") <- converged
  w
}
