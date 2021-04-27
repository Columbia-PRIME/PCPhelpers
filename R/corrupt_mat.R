#' Corrupts a data matrix by subjecting specified columns to a given limit of detection.
#'
#' \code{corrupt_mat} subjects specified columns of a given data matrix to a given
#' limit of detection (LOD). Follows a given imputation scheme when imputing values below LOD.
#'
#' @param mat The data matrix to be corrupted.
#' @param cols The columns in the data matrix to subject to an LOD.
#' @param limit The percent of each specified column to put under the LOD.
#'   Ex: if \code{limit = 0.25} then the first quartile of each specified column is put under the LOD.
#' @param fill The impute scheme to use for values below LOD. Takes one of the following:
#'   \code{"sqrt2", "-1", "mean", "NA", "status"} (Default = \code{"NA"}). See below for descriptions of each.
#'
#' @section Options for \code{fill}:
#' \describe{
#'   \item{\code{"sqrt2"}}{results in values below the LOD imputed as \code{LOD/sqrt\{2\}}.}
#'   \item{\code{"-1"}}{results in values below the LOD imputed as \code{-1}.}
#'   \item{\code{"mean"}}{imputes values below the LOD with the mean of those values below the LOD.}
#'   \item{\code{"NA"}}{replaces values below the LOD with \code{NA} missing values.}
#'   \item{\code{"status"}}{is a special option that returns a vector detailing whether or not values are below the given LOD or not.}
#' }
#' Note that the LOD is column-specific.
#' @return The corrupted matrix
#' @seealso \code{\link{impute}}
#' @examples
#' mat <- matrix(rnorm(25), 5, 5)
#' corrupted_mat <- corrupt_mat(mat, cols=1:5, limit=.4, fill="-1")
#' @export
corrupt_mat <- function(mat, cols, limit, fill="NA") {
  mat[, cols] <- apply(mat[, cols, drop=FALSE], FUN=impute, limit=limit, fill=fill, MARGIN=2)
  mat
}