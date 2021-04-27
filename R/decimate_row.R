#' Decimates a given row of a data matrix based on the patterns the data matrix was simulated from.
#'
#' Given a row vector, \eqn{\beta}, simulated from a ground-truth pattern matrix, \eqn{P},
#' \code{decimate_row} corrupts \eqn{\beta} to produce \eqn{\beta}', such that \eqn{\beta}'
#' has every value missing besides exactly one random observation corresponding to each pattern in \eqn{P}
#' (where the patterns are the rows of \eqn{P}).
#'
#' @param row The row to be decimated.
#' @param pats The ground truth patterns the row was simulated from.
#' @return The decimated row.
#' @seealso \code{\link{decimate_mat}}
#' @keywords internal
decimate_row <- function(row, pats) {
  ret <- rep(-1, length(row)) #or ncol(pats)
  for (r in 1:nrow(pats)) {
    block <- which(pats[r,] %in% c(1))
    pick <- sample(block, size=1)
    ret[pick] <- row[pick]
  }
  ret
}