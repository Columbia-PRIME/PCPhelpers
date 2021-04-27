#' Decimates a data matrix based on the patterns the data matrix was simulated from.
#'
#' Given a data matrix, \eqn{M}, simulated from a ground-truth pattern matrix, \eqn{P},
#' \code{decimate_mat} corrupts \eqn{M} row by row to produce \eqn{M}', such that every row of \eqn{M}'
#' has every value missing besides exactly one random observation corresponding to each pattern in \eqn{P}
#' (where the patterns are the rows of \eqn{P}).
#'
#' @param mat The matrix to be decimated.
#' @param patterns The ground truth patterns the matrix was simulated from.
#' @return The decimated matrix
#' @seealso \code{\link{decimate_row}}
#' @examples
#' pattern <- t(cbind(c(1,1,1,0,0), c(0,0,0,1,1)))
#' scores <- matrix(rnorm(10), 5, 2)
#' sim_mat <- scores %*% pattern
#' decimated_mat <- decimate_mat(sim_mat,  patterns=pattern)
#' @export
decimate_mat <- function(mat, patterns) {
  col_names <- colnames(mat)
  dm <- apply(mat, MARGIN=1, FUN=decimate_row, pats=patterns)
  dm <- t(dm)
  colnames(dm) <- col_names
  dm
}
