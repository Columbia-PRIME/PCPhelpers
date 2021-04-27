#' Subjects a given vector to a given limit of detection.
#'
#' \code{impute} puts a given vector (usually a column vector) under a given limit of detection
#' (LOD) by imputing values below the LOD with a given impute scheme.
#'
#' @param v The vector to impute (usually a column vector in a matrix).
#' @param limit The percent of the vector to put under the LOD.
#'   Ex: if \code{limit = 0.25} then the first quartile of the vector is put under the LOD.
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
#' @return The imputed vector.
#' @seealso \code{\link{corrupt_mat}}
#' @examples
#' mat <- matrix(rnorm(25), 5, 5)
#' mat[,1] <- impute(mat[,1], limit=0.4, fill="-1")
#' @keywords internal
impute <- function(v, limit, fill="NA") {
  q = quantile(v, probs = limit)
  if (fill == "sqrt2") {
    rep_val = q/sqrt(2)
  } else if (fill == "-1") {
    rep_val = -1
  } else if (fill == "mean") {
    rep_val = mean(v[v <= q])
  } else {
    rep_val = NA
  }
  # annoying edge cases with quantile
  if (limit == 0) {
    x = ifelse(v < q, rep_val, v)
  } else {
    x = ifelse(v <= q, rep_val, v)
  }
  if (fill == "status") {
    if (limit == 0) {
      x = ifelse(v < q, "< LOD", "above")
    } else {
      x = ifelse(v <= q, "< LOD", "above")
    }
  }
  x
}
