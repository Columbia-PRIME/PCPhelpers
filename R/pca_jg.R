#' Crude version of PCA
#'
#' \code{pca_jg} performs a crude PCA given a matrix and the rank of the matrix.
#'
#' @param mat The matrix to perform PCA upon.
#' @param rank The rank of the given matrix.
#' @param true_pattern The pattern the matrix was simulated from (optional).
#' @param cov_method The method with which to calculate the covariance matrix. Only relevant when \code{true_pattern} is left \code{NULL}.
#'   Default=\code{"complete.obs"}. See \code{\link[stats]{cov}}.
#'   Will likely need to pass \code{"pairwise.complete.obs"} for matrices with many missing values.
#' @return List containing: \code{"fitted_mat", "scores", "pattern"}.
#' \describe{
#'   \item{\code{"fitted_mat"}}{The recovered matrix (after regressing rows on patterns to estimate scores).}
#'   \item{\code{"scores"}}{The estimated scores (estimated by regressing rows on patterns).}
#'   \item{\code{"pattern"}}{The estimated pattern (if none was provided).}
#' }
#' @section Methods:
#' If \code{true_pattern} is left \code{NULL} then the pattern is estimated via an
#' eigen decomposition on the covariance matrix of \code{mat}.
#' The scores are estimated by regressing rows on the pattern provided/estimated.
#' The fitted values are then estimated by multiplying the estimated scores by the provided/estimated patterns.
#' @examples
#' # simulate a matrix:
#'
#' pattern <- t(cbind(c(1,1,1,0,0), c(0,0,0,1,1)))
#' scores <- matrix(rnorm(10), 5, 2)
#' mat <- scores %*% pattern
#' 
#' # run PCA:
#' 
#' pca_jg(mat, rank = 2, true_pattern = pattern)
#' @export
#' @importFrom magrittr %>%
pca_jg = function(mat, rank = 1, true_pattern = NULL, cov_method="complete.obs") {

  I = dim(mat)[1]
  D = dim(mat)[2]

  data = replace(mat, mat == -1, NA)

  if (is.null(true_pattern)) {

    mean_vec = rep(0, length = ncol(mat))
    data.tilde = data - matrix(mean_vec, I, D, byrow = TRUE)

    # when you have enough missing values, you can't compute covariances with "complete.obs" anymore.
    eigen_decomp = cov(data.tilde, use=cov_method) %>% eigen
    #cov(data.tilde, use = "complete.obs") %>% eigen

    model_patterns = eigen_decomp$vectors[, 1:rank]

  }

  if (!is.null(true_pattern)) {
    mean_vec = rep(0, length = ncol(mat))
    data.tilde = data - matrix(mean_vec, I, D, byrow = TRUE)

    model_patterns = t(true_pattern)[, 1:rank]
  }

  est_scores = matrix(NA, I, rank)
  est_mat = matrix(NA, I, D)
  colnames(est_mat) = colnames(mat)

  for (i in 1:I) {
    est_scores[i, ] = lm(data.tilde[i, ] ~ 0 + model_patterns) %>% coef
    est_mat[i, ] = mean_vec + est_scores[i, ] %*% t(model_patterns)
  }

  ret = list(fitted_mat = est_mat, scores = est_scores, pattern = model_patterns)
  ret
}
