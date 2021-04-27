#' Conducts a cross-validated randomized search of the parameters for Principle Component Pursuit (PCP).
#'
#' \code{random_search_cv} conducts a cross-validated randomized search of the
#' parameters for a given PCP function, given a data matrix \code{mat} and parameter settings to search through.
#' See the \strong{Methods} section below for more details.
#'
#' @param mat The data matrix to conduct the grid search on.
#' @param pcp_func The PCP function to use when grid searching. \emph{Note: the PCP function passed must be able to handle missing \code{NA} values.} For example: \code{root_pcp_na}. 
#' @param grid_df A dataframe with dimension N x P containing the N-many settings of P-many parameters to try. \emph{The columns of grid_df should be named exactly as they are in the function 
#' header of} \code{pcp_func}. For example, if \code{pcp_func = root_pcp_noncvx_na}, then the columns of \code{grid_df} should be named "lambda", "mu", and "r" 
#' (assuming you want to search all 3 parameters; if one of those parameters is constant, instead of giving it its own column in \code{grid_df}, you can simply pass it as a free argument to this method. See \code{...} below).  
#' An optional additional column named "value" can be included that contains the mean relative errors recorded by that row's parameter setting, with those rows (settings) that have not been tried left as \code{NA}. 
#' In this way, you can perform a grid search in which you already know the relative errors of some parameter settings, but would like to expand your knowledge of the unexplored parts of the grid further. 
#' Ex: conduct a a bayesian grid search, examining 10/50 settings. Then search again, looking at another 10 settings, but including the information learned from the first run.
#' @param n_evals The number of parameter settings in \code{grid_df} you would like to evaluate. 
#' @param cores The number of cores to use when parallelizing the grid search. If \code{cores = 1}, the search will be conducted sequentially. If \code{cores > 1}, then the search will be parallelized.
#' By default, \code{cores =} the maximum available cores on your machine. For optimal performance, \code{cores} should usually be set to half that.
#' @param perc_b The percent of entries of the matrix \code{mat} that will be randomly imputed as \code{NA} missing values. By default, \code{perc_b = 0.2}.
#' @param runs The number of times to test a given parameter setting. By default, \code{runs = 100}. 
#' @param seed The seed used when randomly selecting parameter settings to evaluate. By default, \code{seed = NULL} to simulate randomness. For reproducible results, set \code{seed} to some whole number.   
#' @param progress_bar An optional logical indicating if you would like a progress bar displayed or not. By default, \code{progress_bar = TRUE}.
#' @param file An optional character containing the file path used to save the output in. Should end in "\code{.Rda}". When \code{file = NULL}, the output is not saved. By default, \code{file = NULL}.
#' @param ... \emph{Any parameters required by \code{pcp_func} that were not specified in \code{grid_df}, and therefore are kept constant (not involved in the grid search). 
#' An example could be the \code{LOD} parameter for those PCP functions that require the \code{LOD} argument.} 
#' 
#' @section Methods:
#' Each hyperparameter setting is cross-validated by:
#' \enumerate{
#'   \item Randomly corrupting \code{perc_b} percent of the entries in \code{mat} as missing (i.e. \code{NA} values), yielding \code{corrupted_mat}.
#'   Done via the \code{\link{corrupt_mat_randomly}} function.
#'   \item Running a PCP function (\code{pcp_func}) on \code{corrupted_mat}, giving \code{L_hat} and \code{S_hat}.
#'   \item Recording the relative recovery errors of \code{L_hat + S_hat} compared with the raw input data matrix for only those values that were imputed as missing during the corruption step.
#'   \item Repeating steps 1-3 for a total of \code{runs} many times.
#'   \item Reporting the mean of the \code{runs}-many runs for each parameter setting.
#' }
#'
#' @return A list containing the following:
#' \describe{
#'    \item{\code{raw}}{a \code{data.frame} containing the raw statistics of each run comprising the grid search. 
#'    These statistics include the parameter settings for the run, 
#'    the random \code{seed} used for the corruption step outlined in step 1 of the \strong{Methods} section below, 
#'    the relative error for the run, the rank of the recovered L matrix, the sparsity of the recovered S matrix, 
#'    and the number of iterations PCP took to reach convergence (20,000 = Did not converge as of PCPhelpers v. 0.3.1).}
#'    \item{\code{formatted}}{A \code{data.frame} containing the summary of the grid search. 
#'    Made to easily pass on to \code{\link{print_gs}}.}
#'    \item{\code{constants}}{A list containing those arguments initially passed as constant values when calling \code{random_search_cv}.}
#' }
#' @seealso \code{\link{grid_search_cv}}, \code{\link{bayes_search_cv}}, and \code{\link{print_gs}}
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
#' param_grid.out <- random_search_cv(mat, pcp_func = root_pcp_na, grid_df = param_grid, n_evals = 4, cores = 4, perc_b = 0.2, runs = 20, seed = 1, progress_bar = TRUE, file = NULL)
#' 
#' # visualize the output:
#' 
#' print_gs(param_grid.out$formatted)
#' @export
#' @importFrom magrittr %>%
#' @importFrom foreach foreach %dopar%
random_search_cv <- function(
  mat, 
  pcp_func,
  grid_df, 
  n_evals,
  cores = NULL,
  perc_b = 0.2,
  runs = 100,
  seed = NULL,
  progress_bar = TRUE,
  file = NULL,
  ...) 
{

  # setting up the parallel programming
  if (is.null(cores)) {
    cores <- parallel::detectCores()
  } 
  cl <- snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)

  # initialization
  if (!is.null(seed)) set.seed(seed)

  metrics <- c("value", "L.rank", "S.sparsity", "iterations")
  for (metric in metrics) {
    if (!metric %in% names(grid_df)) {
      grid_df[, metric] <- as.numeric(NA)
    }
  }

  points_to_eval <- which(is.na(grid_df$value))
  if (n_evals > length(points_to_eval)) {
    n_evals <- length(points_to_eval)
  }

  metrics <- c("value", "L.rank", "S.sparsity", "iterations")
  param_names <- grid_df %>% dplyr::select(!tidyselect::all_of(metrics)) %>% colnames()
  points <- sample(points_to_eval, size = n_evals)
  params <- data.frame(grid_df[points, param_names], rep(1:runs, each = length(points)), row.names = NULL)
  colnames(params) <- c(param_names, "seed")
  constant_params <- list(...)
  
  # progress bar setup
  if (progress_bar) {
    pb <- txtProgressBar(min=0, max=nrow(params), width=50, style=3)
    progress <- function(p) setTxtProgressBar(pb, p)
    opts <- list(progress=progress)
  } else {
    opts <- list()
  }
  
  # grid search
  cv <- foreach(i = iterators::icount(nrow(params)), .options.snow=opts, .combine = cbind, .packages = c("pcpr", "PCPhelpers"), .inorder = FALSE) %dopar% {
    do.call(eval_params, c(as.list(params[i,]), constant_params, list(mat = mat, pcp_func = pcp_func, perc_b = perc_b, eval_params = param_names)))
  }

  # close the progress bar and stop cluster
  if (progress_bar) close(pb)
  snow::stopCluster(cl)

  # format the output
  rownames(cv) <- c(colnames(params), metrics)
  cv <- as.data.frame(t(cv))

  cv.formatted <- cv %>% 
    dplyr::group_by_at(param_names) %>% 
    dplyr::summarise(value = mean(value), L.rank = mean(L.rank), S.sparsity = mean(S.sparsity), iterations = mean(iterations)) %>% 
    tidyr::unite(param_setting, param_names)

  grid_df.formatted <- grid_df %>% 
    tidyr::unite(param_setting, param_names)

  grid_df.formatted[match(cv.formatted$param_setting, grid_df.formatted$param_setting), ] <- cv.formatted

  grid_df <- grid_df.formatted %>% 
    tidyr::separate(param_setting, into = param_names, sep = "_", convert = TRUE)
  
  # save
  if (!is.null(file)) save(cv, grid_df, file = file)

  # return
  list(raw = cv, formatted = grid_df, constants = constant_params)
}