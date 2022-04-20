#' Evaluates a given setting of lambda and mu on a given matrix with a given PCP function. 
#' 
#' \code{eval_params} evaluates a given setting of lambda and mu on a given matrix with a given PCP function. 
#' The given matrix has a given percentage of its entries randomly corrupted as missing \code{NA} values before evaluation.
#' The parameter setting is scored by how well the given PCP function recovers those randomly corrupted values. 
#'
#' @param lambda The value of lambda to be passed to \code{pcp_func}.
#' @param mu The value of mu to be passed to \code{pcp_func}.
#' @param seed The seed to use for the random corruption of the given matrix, \code{mat}.
#' @param mat The data matrix to run PCP on.
#' @param pcp_func The PCP function to use. \emph{Note: the PCP function passed must be able to handle missing \code{NA} values.} For example: \code{root_pcp_na}.
#' @param perc_b The percentage of \code{mat} to randomly corrupt as missing. 
#' @param eval_params A character vector containing the names of the parameters under evaluation.
#' @param ... The parameters to pass on to \code{pcp_func}.
#'
#' @return A vector of length containing the values of the parameters used, along with evaluation metrics.
#' @seealso \code{\link{corrupt_mat_randomly}}, \code{\link{grid_search_cv}}, \code{\link{random_search_cv}}, \code{\link{bayes_search_cv}}
#' @examples
#' 
#' library(pcpr) # since we will be passing grid_search_cv a PCP function 
#' 
#' # simulate a data matrix:
#' 
#' n <- 50
#' p <- 10
#' data <- sim_data(sim_seed = 1, nrow = n, ncol = p, rank = 3, sigma=0, add_sparse = FALSE)
#' mat <- data$M
#'
#' # pick a parameter setting of lambda and mu to try:
#' 
#' lambda <- 1/sqrt(n)
#' mu <- sqrt(p/2)
#'
#' # evaluate that setting:
#' 
#' score <- eval_params(lambda, mu, 1, mat, root_pcp_na, .2)
#' @keywords internal
eval_params <- function(
  seed,
  mat,
  pcp_func,
  perc_b,
  eval_params,
  ...
) {

  pcp_args <- list(...)
  
  cor_mat <- corrupt_mat_randomly(mat, k = seed, perc_b = perc_b) # a list containing cor.mat and cor.mask
  mask <- cor_mat$cor.mask
  pcp_out <- do.call(pcp_func, c(pcp_args, list(D = cor_mat$cor.mat, verbose = F)))
  mat[is.na(mat)] <- 0
  score <- norm((mat - pcp_out$L - pcp_out$S)*mask, "F") / norm(mat * mask, "F")
  L.rank <- Matrix::rankMatrix(pcp_out$L, tol = 1e-04)
  S.sparsity <- sparsity(pcp_out$S, tol = 1e-04)
  its <- pcp_out$final_iter

  c(as.numeric(pcp_args[eval_params]), seed, score, L.rank, S.sparsity, its)
}