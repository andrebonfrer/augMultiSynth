#' Screen donors for a treated unit using pre-treatment similarity
#'
#' Selects a subset of eligible donors based on similarity to the treated unit in
#' pre-treatment periods. This is a performance-critical step for large panels.
#'
#' @param Y_ref A numeric matrix \code{N x T} for a reference outcome used to screen donors.
#'   Typically choose the most informative or least noisy outcome.
#' @param treat_time Numeric vector length \code{N}, first treated period per unit (or \code{Inf}).
#' @param j Integer index of treated unit.
#' @param K Integer. Post window size used for donor eligibility (donors must remain untreated
#'   through \code{treat_time[j] + K}).
#' @param L Integer. Number of pre-treatment periods used for screening.
#' @param max_donors Integer. Maximum number of donors to return.
#' @param method Screening method: \code{"cor"} (correlation) or \code{"mse"} (mean squared error).
#' @param candidates Optional integer vector of unit indices. When given, only these units are
#'   considered (after the eligibility rule below); \code{NULL} (default) considers every unit.
#'   Used to impose a design-specific donor pool, e.g. donors adopting in a band of weeks around
#'   the treated unit's relative adoption week.
#'
#' @details
#' Eligibility rule implemented:
#' \itemize{
#'   \item Donors are units with \code{treat_time > treat_time[j] + K} (including never-treated).
#'   \item With \code{candidates}, the eligible set is intersected with it before screening.
#' }
#'
#' @return Integer vector of donor indices (subset of eligible donors), length \code{<= max_donors}.
#'
#' @examples
#' N <- 500; T <- 120
#' Y <- matrix(rnorm(N*T), N, T)
#' treat_time <- rep(Inf, N); treat_time[1:100] <- 80
#' donors <- donor_screen(Y, treat_time, j = 1, K = 10, L = 20, max_donors = 200)
#' length(donors)
#'
#' @export
donor_screen <- function(Y_ref, treat_time, j, K, L, max_donors = 1000,
                         method = c("cor", "mse"), candidates = NULL) {

  method <- match.arg(method)

  N <- nrow(Y_ref); T <- ncol(Y_ref)
  Tj <- treat_time[j]
  if (!is.finite(Tj)) stop("Unit j must be treated (finite treat_time[j]).")

  pre_times <- (Tj - L):(Tj - 1)
  if (min(pre_times) < 1) stop("Not enough pre periods for screening.")

  eligible <- which(treat_time > (Tj + K))
  if (!is.null(candidates)) {
    candidates <- as.integer(candidates)
    candidates <- candidates[is.finite(candidates) & candidates >= 1L & candidates <= N]
    eligible <- intersect(eligible, candidates)
  }
  eligible <- setdiff(eligible, j)
  if (length(eligible) == 0) stop("No eligible donors under the eligibility rule.")

  yj <- Y_ref[j, pre_times]
  Yd <- Y_ref[eligible, pre_times, drop = FALSE]

  score <- switch(
    method,
    cor = {
      s <- suppressWarnings(apply(Yd, 1, function(x) cor(x, yj)))
      s[is.na(s)] <- -Inf
      s
    },
    mse = {
      # negative MSE so "larger is better" for ordering
      -rowMeans((Yd - matrix(yj, nrow = nrow(Yd), ncol = length(yj), byrow = TRUE))^2)
    }
  )

  keep <- order(score, decreasing = TRUE)[seq_len(min(max_donors, length(eligible)))]
  eligible[keep]
}
