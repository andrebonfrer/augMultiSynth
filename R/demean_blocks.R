#' Demean within outcome blocks (per-outcome intercept) with optional time weights
#'
#' @param X Numeric matrix of dimension (M*L) x Nd (stacked by outcome blocks of size L).
#' @param y Numeric vector of length (M*L) (stacked by outcome blocks of size L).
#' @param L Integer. Block size (number of pre-periods).
#' @param time_weights Optional numeric vector length L (nonnegative). Used for weighted means within each block.
#'
#' @return A list with demeaned X and y.
#' @export
demean_within_outcome_blocks <- function(X, y, L, time_weights = NULL) {
  n <- nrow(X)
  if (length(y) != n) stop("y must have length nrow(X).")
  if (n %% L != 0) stop("nrow(X) must be a multiple of L.")
  M <- n %/% L

  if (is.null(time_weights)) {
    tw <- rep(1, L)
  } else {
    tw <- as.numeric(time_weights)
    if (length(tw) != L) stop("time_weights must have length L.")
    if (any(!is.finite(tw)) || any(tw < 0)) stop("time_weights must be finite and nonnegative.")
    if (sum(tw) <= 0) stop("time_weights must sum to > 0.")
  }

  X_out <- X
  y_out <- y

  for (m in seq_len(M)) {
    idx <- ((m - 1) * L + 1):(m * L)

    # Weighted mean across time (rows) within this outcome block
    sw <- sum(tw)

    # y block mean
    ybar <- sum(tw * y_out[idx]) / sw
    y_out[idx] <- y_out[idx] - ybar

    # X block column means (1 x Nd): weighted mean over rows
    # crossprod(tw, X_block) gives 1 x Nd
    Xbar <- as.numeric(crossprod(tw, X_out[idx, , drop = FALSE]) / sw)

    X_out[idx, ] <- X_out[idx, , drop = FALSE] -
      matrix(Xbar, nrow = L, ncol = ncol(X_out), byrow = TRUE)
  }

  list(X = X_out, y = y_out)
}
