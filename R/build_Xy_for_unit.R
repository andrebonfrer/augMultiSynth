#' Build stacked pre-treatment design for one treated unit
#'
#' Constructs the stacked (multi-outcome) pre-treatment design matrix and target
#' vector for a single treated unit \code{j}. This is used to fit donor weights
#' that minimize pre-treatment mismatch across multiple outcomes simultaneously.
#'
#' @param Y_list A list of outcome matrices, one per outcome.
#'   Each element must be a numeric matrix of dimension \code{N x T} with units
#'   in rows and time in columns. All outcomes must have identical \code{N} and \code{T}.
#' @param treat_time Numeric vector of length \code{N}. \code{treat_time[i]}
#'   is the first treated period for unit \code{i}. Use \code{Inf} for never-treated units.
#' @param j Integer. Index of the treated unit for which to build the design.
#' @param donors Integer vector. Indices of donor units eligible for unit \code{j}.
#' @param L Integer. Number of pre-treatment periods to use. Uses periods
#'   \code{(treat_time[j]-L):(treat_time[j]-1)}.
#' @param standardize_outcomes Logical. If \code{TRUE}, standardizes each outcome block
#'   using donor-pool pre-treatment mean and SD (computed over donors and \code{pre_times}),
#'   then applies the same transformation to treated and donor pre-treatment values.
#'   This helps when outcomes are on different scales or have outcome-specific intercepts.
#' @param eps_sd Small positive scalar used as a floor for the SD to avoid division by
#'   near-zero values when standardizing.
#'
#' @details
#' The function stacks pre-treatment observations over outcomes and time:
#' \itemize{
#'   \item \code{y} is length \code{M*L}, containing treated unit \code{j}'s pre-treatment
#'         values for each outcome and time
#'   \item \code{X} is a \code{(M*L) x Nd} matrix where each column is a donor unit's
#'         stacked pre-treatment values over outcomes and time
#' }
#'
#' When \code{standardize_outcomes=TRUE}, each outcome \code{m} is transformed as
#' \deqn{Y^{(m)} \leftarrow (Y^{(m)} - \bar Y^{(m)}_{\text{donor,pre}})/s^{(m)}_{\text{donor,pre}}}
#' where the mean and SD are computed using donor units and pre-treatment periods only,
#' avoiding post-treatment leakage.
#'
#' @return A list with components:
#' \describe{
#'   \item{X}{Numeric matrix of dimension \code{(M*L) x Nd}, stacked donor pre outcomes.}
#'   \item{y}{Numeric vector of length \code{M*L}, stacked treated unit pre outcomes.}
#'   \item{pre_times}{Integer vector of pre-treatment time indices used.}
#'   \item{scales}{List of length \code{M} with per-outcome \code{mean} and \code{sd} used
#'     when \code{standardize_outcomes=TRUE}. If \code{FALSE}, returns mean=0 and sd=1.}
#' }
#'
#' @examples
#' N <- 50; T <- 60
#' Y1 <- matrix(rnorm(N*T), N, T)
#' Y2 <- matrix(10 + 5*rnorm(N*T), N, T)   # different scale/level
#' treat_time <- rep(Inf, N); treat_time[1:10] <- 40
#' donors <- which(treat_time > 50)
#' out <- build_Xy_for_unit(
#'   list(Y1, Y2), treat_time, j = 1, donors = donors, L = 10,
#'   standardize_outcomes = TRUE
#' )
#' str(out)
#'
#' @export
build_Xy_for_unit <- function(Y_list, treat_time, j, donors, L,
                              standardize_outcomes = FALSE,
                              eps_sd = 1e-8) {

  M <- length(Y_list)
  if (M < 1) stop("Y_list must contain at least one outcome matrix.")
  N <- nrow(Y_list[[1]])
  TT <- ncol(Y_list[[1]])

  if (length(treat_time) != N) stop("treat_time must have length N.")
  for (m in seq_len(M)) {
    if (!is.matrix(Y_list[[m]])) stop("Each element of Y_list must be a matrix.")
    if (nrow(Y_list[[m]]) != N || ncol(Y_list[[m]]) != TT) {
      stop("All outcome matrices must have identical dimensions N x T.")
    }
  }

  Tj <- treat_time[j]
  if (!is.finite(Tj)) stop("Unit j must be treated (finite treat_time[j]).")
  Tj <- as.integer(Tj)

  pre_times <- (Tj - L):(Tj - 1)
  if (min(pre_times) < 1) stop("Not enough pre periods for unit j given L.")
  if (max(pre_times) > TT) stop("Pre period indexing exceeds T.")

  y <- numeric(M * L)
  X <- matrix(0.0, nrow = M * L, ncol = length(donors))
  scales <- vector("list", M)

  row <- 1
  for (m in seq_len(M)) {
    Y <- Y_list[[m]]

    # donor and treated pre blocks
    D_pre <- Y[donors, pre_times, drop = FALSE]  # Nd x L
    y_pre <- as.numeric(Y[j, pre_times])         # length L

    if (isTRUE(standardize_outcomes)) {
      mu_m <- mean(D_pre)
      sd_m <- stats::sd(as.numeric(D_pre))
      if (!is.finite(sd_m) || sd_m < eps_sd) sd_m <- eps_sd

      D_pre <- (D_pre - mu_m) / sd_m
      y_pre <- (y_pre - mu_m) / sd_m

      scales[[m]] <- list(mean = mu_m, sd = sd_m)
    } else {
      scales[[m]] <- list(mean = 0, sd = 1)
    }

    y[row:(row + L - 1)] <- y_pre
    X[row:(row + L - 1), ] <- t(D_pre)  # L x Nd

    row <- row + L
  }

  list(X = X, y = y, pre_times = pre_times, scales = scales)
}
