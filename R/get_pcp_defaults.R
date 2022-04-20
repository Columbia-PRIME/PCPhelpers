#' Retrieves the default Root PCP parameter settings for a given matrix.
#'
#' \code{get_pcp_defaults} retrieves the default Root PCP lambda and mu parameter settings for a given matrix.
#' 
#' @param mat The matrix for which the default lambda and mu are needed. 
#'
#' @return A list containing the default lambda and mu values used in most PCP functions, and the default eta value used in \code{\link{RRMC}}. labelled as "lambda", "mu", and "eta" respectively.
#' @examples
#'
#' # simulate a data matrix:
#' 
#' n <- 50
#' p <- 10
#' data <- sim_data(sim_seed = 1, nrow = n, ncol = p, rank = 3, sigma=0, add_sparse = FALSE)
#' mat <- data$M
#'
#' # get the default PCP parameters:
#' default <- get_pcp_defaults(mat)
#' default_lambda <- default$lambda
#' default_mu <- default$mu
#' @export
get_pcp_defaults <- function(mat) {
	n <- nrow(mat)
	p <- ncol(mat)
	l <- 1/sqrt(max(n, p))
	m <- sqrt(min(n, p)/2)
	e <- l / m
	list(lambda = l, mu = m, eta = e)
}