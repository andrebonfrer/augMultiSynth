#' Internal helper: run a function over X in a chosen parallel backend
#'
#' @param X Vector/list of inputs.
#' @param FUN Function of one argument.
#' @param parallel_backend One of "none","fork","psock".
#' @param n_cores Integer number of workers/cores.
#' @param export Character vector of object names to export (psock only).
#' @param envir Environment to export from (psock only).
#'
#' @return List of length length(X).
#' @keywords internal
run_parallel <- function(X, FUN,
                         parallel_backend = c("none","fork","psock"),
                         n_cores = 2L,
                         export = NULL,
                         envir = parent.frame()) {

  parallel_backend <- match.arg(parallel_backend)
  n_cores <- as.integer(n_cores)
  if (n_cores < 1L) n_cores <- 1L

  if (parallel_backend == "none" || n_cores == 1L) {
    return(lapply(X, FUN))
  }

  if (parallel_backend == "fork") {
    if (.Platform$OS.type == "windows") {
      stop("parallel_backend='fork' is not supported on Windows; use 'psock'.")
    }
    return(parallel::mclapply(
      X, FUN,
      mc.cores = n_cores,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE
    ))
  }

  # psock backend
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(export)) {
    parallel::clusterExport(cl, varlist = export, envir = envir)
  }

  # IMPORTANT: FUN passed exactly once, as the 3rd argument (positional)
  parallel::parLapply(cl, X, FUN)
}
