# imputation function for the hybrid method (pre-processing step...)
# D = matrix to be hybridized
# r is the rank of D
# limit = LOD

#' Pre-processing step for PCP_hybrid
#'
#' \code{hybridize} computes the necessary information for the PCP_hybrid method.
#'
#' @param D The data matrix to hybridize.
#' @param r The rank of the data matrix.
#' @param limit The percent of each column in \code{D} to put under the LOD.
#'   Ex: if \code{limit = 0.25} then the first quartile of each column in \code{D}
#'   is put under the LOD. Avoid passing \code{limit = 0} (see \emph{Warnings} below).
#'
#' @return List containing: \code{"M_hybrid", "is_safe", "below_lod_mat"}.
#' See \strong{Methods} below for further details.
#' \describe{
#'   \item{\code{"M_hybrid"}}{The hybridized matrix.}
#'   \item{\code{"is_safe"}}{A logical vector of length \code{nrows(D)}.
#'   Entries are \code{TRUE} when the corresponding row in \code{D} is a "safe row",
#'   and \code{FALSE} when the corresponding row in \code{D} is "unsafe".}
#'   \item{\code{"below_lod_mat"}}{A binary matrix, where \code{1}'s signify the corresponding
#'   entry in \code{D} was below LOD, and \code{0}'s signify the corresponding entry was above LOD.}
#' }
#' @section Methods:
#' The main idea with PCP_hybrid is to make use of information on "safe" vs. "unsafe" rows in a data matrix.
#' A \emph{safe row} is defined as a row with at least \code{r}-many entries above the LOD.
#' An \emph{unsafe row} is a row with less than \code{r}-many entries above the LOD.
#' \strong{When hybridizing a matrix, if an entry is below the LOD in a safe row, it is imputed as \code{-1},
#' whereas if it is below the LOD in an unsafe row, it is instead imputed as \code{LOD/sqrt\{2\}}.}
#' @section Warnings:
#' Do not pass \code{limit = 0} to \code{hybridize}, as the underlying \code{\link[stats]{quantile}}
#' function will result in some values put under the lod anyway (since passing 0 to
#' \code{\link[stats]{quantile}}) results in the minimum value selected as the LOD.
#' @examples
#' data <- sim_data(sim_seed = 1, nrow = 10, ncol = 10, rank = 3, sigma=0, add_sparse = FALSE)
#' mat <- data$M
#' hybridize(mat, r = 3, limit = 0.25)
#' @export
hybridize <- function(D, r, limit) {

  delta <- as.vector(apply(D, 2, quantile, probs = limit))
  B <- t(D) <= delta
  B <- t(B) * 1 # 1 is below LOD, 0 is above
  Bv <- ncol(D) - as.vector(apply(B, 1, sum))
  is_safe <- Bv >= r

  M <- matrix(NA, nrow = nrow(D), ncol = ncol(D))
  D <- as.matrix(D)
  for (i in 1:ncol(D)) {
    q = delta[i]
    for (j in 1:nrow(D)) {
      if (D[j, i] <= q) {
        M[j, i] <- ifelse(is_safe[j], -1, q/sqrt(2))
      } else {
        M[j, i] <- D[j, i]
      }
    }
  }
  ret <- list(M_hybrid=M, is_safe=is_safe, below_lod_mat=B)
  ret
}
