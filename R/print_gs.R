#' Visualizes the \code{formatted} output from \code{\link{bayes_search_cv}}, \code{\link{grid_search_cv}}, or \code{\link{random_search_cv}}.
#'
#' \code{print_gs} visualizes the \code{formatted} output grid searches in the form of heatmaps. 
#'
#' @param df The \code{formatted} output from \code{\link{bayes_search_cv}}, \code{\link{grid_search_cv}}, or \code{\link{random_search_cv}}.
#' @param dlambda A numeric providing the default value of lambda used in the grid search.
#' @param dmu A numeric providing the default value of mu used in the grid search.
#' @param title A character providing the title for the plot. By default, \code{title = "Lambda and Mu Grid Search"}
#' 
#' @return A \code{\link{heatmaply}} plot. 
#' @seealso \code{\link{bayes_search_cv}}, \code{\link{grid_search_cv}}, \code{\link{random_search_cv}}
#' @examples
#' 
#' library(pcpr) # since we will be passing grid_search_cv a PCP function 
#' 
#' # simulate a data matrix:
#' 
#' n <- 50
#' p <- 10
#' data <- sim_data(sim_seed = 1, nrow = n, ncol = p, rank = 3, sigma=0, add_sparse = FALSE)
#' mat <- data$M
#'
#' # pick parameter settings of lambda and mu to try:
#' 
#' lambdas <- c(1/sqrt(n), 1.25/sqrt(n), 1.5/sqrt(n))
#' mus <- c(sqrt(p/2), sqrt(p/1.5), sqrt(p/1.25))
#' param_grid <- expand.grid(lambda = lambdas, mu = mus)
#'
#' # run the grid search:
#' 
#' param_grid.out <- bayes_search_cv(mat, pcp_func = root_pcp_na, grid_df = param_grid, init_evals = 3, bayes_evals = 3, cores = 4, acquisition_function = "ei", perc_b = 0.2, runs = 20, seed = 1, verbose = TRUE, file = NULL)
#' 
#' # visualize the output:
#' 
#' print_gs(param_grid.out$formatted)
#' @export
print_gs <- function(df, dlambda = -1, dmu = -1, title="Lambda and Mu Grid Search") {
	lambdas <- sort(unique(df$lambda))
	mus <- sort(unique(df$mu))

	df.mat <- matrix(df$value, nrow = length(lambdas), ncol = length(mus))
	rownames(df.mat) <- lambdas
	colnames(df.mat) <- mus

	dl <- match(dlambda, rownames(df.mat))
  	dm <- match(dmu, colnames(df.mat))

	min_val <- which.min(df.mat)
	cn <- matrix("", nrow = nrow(df.mat), ncol = ncol(df.mat))
    cn[min_val] <- "optimal"
    cn[dl, dm] <- paste(cn[dl, dm], "default")

	heatmaply::heatmaply(df.mat, cellnote = cn, 
		Colv = F, Rowv = F, showticklabels = c(T, T), cellnote_textposition="middle center",
		label_names = c("lambda", "mu", "rel err"), xlab = "mu", ylab = "lambda", main = title)
}