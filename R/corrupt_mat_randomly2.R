#' @export
corrupt_mat_randomly2 <- function(mat, seed, perc) {
  nvals_to_corrupt <- floor(perc*prod(dim(mat)))
  mat_vec <- as.vector(mat)
  mask <- rep(0, length(mat_vec))
  
  pool <- which(mat_vec >= 0)
  if (length(pool) == 0) stop('There is nothing in the input matrix "mat" that can be corrupted as missing.')
  if (nvals_to_corrupt > length(pool)) {
    corrupted <- pool
  } else {
    set.seed(seed)
    corrupted <- sample(pool, nvals_to_corrupt, replace = FALSE)
  }
  
  mask[corrupted] <- 1
  mat_vec[corrupted] <- NA
  
  rows <- nrow(mat)
  cols <- ncol(mat)
  
  ret_mat <- matrix(mat_vec, nrow = rows, ncol = cols)
  ret_mask <- matrix(mask, nrow = rows, ncol = cols)
  
  list(mat_tilde = ret_mat, tilde_mask = ret_mask)
}