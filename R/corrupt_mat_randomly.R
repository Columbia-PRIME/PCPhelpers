#' Corrupts a data matrix by knocking out a random specified percentage of entries.
#'
#' \code{corrupt_mat_randomly} corrupts a given data matrix by knocking out a random
#' specified percentage of entries as \code{NA} missing values.
#'
#' @param mat The data matrix to be corrupted.
#' @param k The random seed to use. 
#' @param perc_b The percentage of the entries in the data matrix to be corrupted. 
#'
#' @return A list containing the corrupted matrix (named \code{cor.mat}), along with a binary matrix (named \code{cor.mask}) where 1's identify entries that were corrupted, and 0
#' identify entries that were left untouched.
#' @seealso \code{\link{corrupt_mat}}
#' @examples
#' mat <- matrix(rnorm(25), 5, 5)
#' corrupted_mat <- corrupt_mat_randomly(mat, k = 1, perc_b = .2)
#' @export
corrupt_mat_randomly <- function(mat, k, perc_b) {
  nvals_to_corrupt <- floor(perc_b*prod(dim(mat)))
  mat.vec <- as.vector(mat)
  mask <- rep(0, length(mat.vec))

  pool <- which(mat.vec >= 0)
  if (nvals_to_corrupt > length(pool)) {
    corrupted <- pool
  } else {
    set.seed(k)
    corrupted <- sample(pool, nvals_to_corrupt, replace = FALSE)
  }

  mask[corrupted] <- 1
  mat.vec[corrupted] <- NA

  rows <- nrow(mat)
  cols <- ncol(mat)

  ret.mat <- matrix(mat.vec, nrow = rows, ncol = cols)
  ret.mask <- matrix(mask, nrow = rows, ncol = cols)

  list(cor.mat = ret.mat, cor.mask = ret.mask)
}