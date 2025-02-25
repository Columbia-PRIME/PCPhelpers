#' Simulate a data matrix of a given dimension and rank for PCP experiments.
#'
#' \code{rsim_data} simulates a data matrix, \code{M}, of a given dimension and rank by simulating
#' ground truth low rank (\code{L}), sparse (\code{S}) and noise (\code{Z}) components.
#' Meaning \code{M = L + S + Z}. The \code{M}, \code{L}, \code{S} and \code{Z} matrices are all accessible from this function,
#' allowing for easy evaluation of PCP's performance.
#'
#' @param sim_seed Required seed to allow for reproducible results & distinct simulated matrices.
#' @param nrow Sets number of rows in the simulated data matrix.
#' @param ncol Sets the number of columns in the simulated data matrix.
#' @param rank Sets the rank of the simulated data matrix.
#' @param sigma Sets the standard deviation of the noise matrix, \code{Z}.
#' @param add_sparse Logical that sets whether sparse noise, \code{S}, is added or not. Default=\code{FALSE}.
#' @param nonneg Logical that sets whether the final M matrix should be non-negative. Default=\code{TRUE}.
#' @param ldist Function that defines a distribution to use when simulating the L matrix. Default=\code{runif}.
#' @param sdist Function that defines a distribution to use when simulating the S matrix. Default=\code{runif}.
#' @param zdist Function that defines a distribution to use when simulating the Z matrix. Default=\code{rnorm}.
#'
#' @return list containing the simulated \code{M}, \code{L}, \code{S} and \code{Z} matrices.
#' @examples
#' data <- rsim_data(sim_seed = 1, nrow = 10, ncol = 10, rank = 3, sigma = .1, add_sparse = TRUE)
#' @export
rsim_data <- function(sim_seed, nrow, ncol, rank, sigma, add_sparse=FALSE, nonneg=TRUE, ldist=runif, sdist=runif, zdist=rnorm, ...) {
  
  set.seed(sim_seed)
  
  # gaussian noise = Z:
  Z <- matrix(zdist(prod(nrow, ncol)), nrow, ncol) * sigma

  # low rank matrix = L:
  U <- matrix(ldist(prod(nrow, rank), ...), nrow, rank)
  V <- matrix(ldist(prod(rank, ncol), ...), rank, ncol)
  L <- U%*%V # ground truth low rank matrix.

  # intermediate step to help define sparse matrix S below:
  D <- L + Z

  # sparse matrix = S:
  if (add_sparse) {
    S <- -D * ((D < 0) * 1) + (matrix(sdist(prod(nrow,ncol), ...), nrow, ncol)<0.03) * matrix(sdist(prod(nrow,ncol), ...), nrow, ncol) * 1 # Jingkai's method of simulating sparse noise (1st term ensures nonnegativity, second term adds some sparse noise to random entries)
    #S <- (-D + rand(nrow, ncol)) * ((D < 0) * 1) # Lawrence's method of simulating sparse noise (use sparse to ensure matrix is non-neg + add extra noise on top of those same entries)
  } else {
    S <- matrix(0, nrow, ncol)
  }

  # final simulated data matrix = M:
  M <- S + D
  if (nonneg) M[M < 0] <- 0 # non-negative

  list(M = M, L = L, S = S, Z = Z)
}
