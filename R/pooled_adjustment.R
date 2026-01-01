#' Fit pooled regression adjustment for synthetic predictions (intercept-only)
#'
#' Given unit-specific synthetic predictions for treated units in their pre-periods,
#' estimate a single pooled intercept adjustment beta0 (per outcome) via ridge:
#'   min_b sum (Y - Ysc - b)^2 + nu * b^2
#'
#' @param y_pre numeric vector of observed outcomes stacked over (i,t) in pre-periods.
#' @param ysc_pre numeric vector of synthetic predictions stacked over (i,t) in pre-periods.
#' @param nu nonnegative ridge penalty for the pooled intercept.
#'
#' @return Estimated intercept adjustment beta0.
#' @keywords internal
fit_pooled_intercept <- function(y_pre, ysc_pre, nu = 0) {
  stopifnot(length(y_pre) == length(ysc_pre))
  r <- y_pre - ysc_pre
  n <- sum(!is.na(r))
  if (n == 0) return(0)
  sum_r <- sum(r, na.rm = TRUE)
  sum_r / (n + nu)
}
