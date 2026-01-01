#' Row-demean a matrix
#'
#' Subtracts the row mean from each element of a matrix. Useful for implementing
#' intercept-shift / demeaned formulations where pre-period levels are removed
#' before fitting synthetic control weights.
#'
#' @param M A numeric matrix.
#'
#' @return A numeric matrix with each row demeaned (row mean equals zero).
#'
#' @examples
#' M <- matrix(1:9, nrow = 3)
#' demean_rows(M)
#'
#' @export
demean_rows <- function(M) {
  M - rowMeans(M)
}
