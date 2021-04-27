#' Rough way to determine the sparsity of a given matrix.
#'
#' \code{sparsity} computes the (crude) sparsity of a given matrix.
#'
#' @param mat The matrix whose sparsity is in question.
#' @param tol (optional) a numeric that determines if an entry in \code{mat} is effectively 0. 
#' If the absolute value of an entry is below \code{tol} then it is effectively 0.
#'
#' @return The sparsity of the given matrix.
#' @examples
#' sparsity(matrix(0:8, 3, 3))
#' @export
sparsity <- function(mat, tol = .00001) {
    mat[abs(mat) < tol] <- 0
    100 - (100 * sum((mat != 0) * 1) / prod(dim(mat)))
}
