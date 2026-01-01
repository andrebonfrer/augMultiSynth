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
#' @param max_iter Maximum number of Frank–Wolfe iterations. Defaults to \code{2000}.
#' @param tol Convergence tolerance on objective improvement. Defaults to \code{1e-6}.
#' @param verbose Logical. If \code{TRUE}, prints occasional progress messages.
#'
#' @details
#' This routine:
#' \itemize{
#'   \item Initializes weights uniformly on donors
#'   \item Computes quadratic objective and gradient in donor space
#'   \item Uses Frank–Wolfe direction to the best simplex vertex (coordinate descent-like)
#'   \item Performs exact line search for the quadratic objective
#' }
#'
#' For best results at scale, screen donors first to a manageable candidate set.
#'
#' @return Numeric vector of length \code{Nd} with simplex weights summing to 1.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200*50), 200, 50)
#' w_true <- rep(0, 50); w_true[c(3,10,25)] <- c(0.2, 0.5, 0.3)
#' y <- X %*% w_true + rnorm(200, sd = 0.1)
#' w_hat <- fit_weights_fw(X, as.numeric(y), lambda = 1e-4)
#' sum(w_hat)  # should be ~1
#'
#' @export
fit_weights_fw <- function(X, y, lambda = 1e-3,
                           max_iter = 2000, tol = 1e-6, verbose = FALSE) {

  if (!is.matrix(X)) stop("X must be a matrix.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y).")
  if (lambda < 0) stop("lambda must be nonnegative.")

  Nd <- ncol(X)
  if (Nd < 2) stop("Need at least 2 donors (ncol(X) >= 2).")

  w <- rep(1 / Nd, Nd)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)

  obj <- function(w) sum((y - X %*% w)^2) + lambda * sum(w^2)
  grad <- function(w) as.numeric(2 * (XtX %*% w - Xty) + 2 * lambda * w)

  f_prev <- obj(w)

  # Use sparse diagonal if Matrix is available; otherwise base diag
  if (requireNamespace("Matrix", quietly = TRUE)) {
    Q <- XtX + lambda * Matrix::Diagonal(Nd)
  } else {
    Q <- XtX + lambda * diag(Nd)
  }

  for (it in seq_len(max_iter)) {
    g <- grad(w)

    # best vertex: put all mass on smallest gradient coordinate
    k <- which.min(g)
    s <- rep(0, Nd); s[k] <- 1
    d <- s - w

    num <- sum(g * d)
    den <- 2 * as.numeric(crossprod(d, Q %*% d))
    gamma <- if (den > 0) min(1, max(0, -num / den)) else 1

    w_new <- w + gamma * d
    f_new <- obj(w_new)

    if (verbose && it %% 200 == 0) {
      message(sprintf("iter %d: obj=%.6f, step=%.4g", it, f_new, gamma))
    }

    if (abs(f_prev - f_new) < tol * (1 + abs(f_prev))) {
      w <- w_new
      break
    }
    w <- w_new
    f_prev <- f_new
  }

  w
}
