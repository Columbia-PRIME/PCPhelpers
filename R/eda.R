#' A rough function to perform exploratory data analysis on a given matrix. Returns correlation matrix and PCA summary figures. 
#'
#' \code{eda} is a (very) rough function to perform exploratory data analysis on a given matrix. Returns correlation matrix and PCA summary figures. 
#' Note: must have \code{ggfortify} package installed and imported before calling this function.
#'
#' @param mat The matrix to perform the eda on.
#' @param cor_method Arguments to be passed to \code{GGally::ggcorr}'s "method" argument. By default, \code{cor_method = c("pairwise.complete.obs", "pearson")}.
#' @param cor_lbl A logical indicating if there should be labels on the correlation matrix. By default, \code{cor_lbl = FALSE}.
#' @param impute What to impute missing values with. By default, \code{impute = 0}.
#' @param scale_flag A logical indicating if the matrix should be scaled before running PCA. By default, \code{scale_flag = TRUE}. 
#' @param pcs A character vector detailing which components you would like to get the loadings plots for. By default, \code{pcs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")}. 
#' @param colgroups An optional dataframe with the column names for \code{mat} in the first column, and their grouping in the second column. By default, \code{colgroups = NULL}.
#' @param rowgroups An optional vector identifying the groups the rows of \code{mat} belong to. Its length should be = \code{nrow(mat)}. By default, \code{rowgroups = NULL}. 
#' @param rowgroups_name A character providing the overall name for the \code{rowgroups} vector.  
#' @param frame An optional logical specifying if you would like frames over the clusters in a PCA biplot. By default, \code{frame = FALSE}. Only worth passing as \code{TRUE} when \code{rowgroups} is not \code{NULL}.
#'
#' @return A list containing the following:
#' \describe{
#'    \item{\code{cor}}{A correlation matrix figure. From \code{GGally::ggcorr}.}
#'    \item{\code{var}}{A summary table/figure describing PCA's principle components and their respective proportion of variance explained. From the \code{kableExtra} package.}
#'    \item{\code{var_raw}}{The raw data underlying the \code{var} figure.}
#'    \item{\code{load}}{A figure of the loadings for the first \code{pcs}-many components PCA has to offer. With the default value of \code{pcs}, you get the loadings for the first six components.}
#'    \item{\code{biplot1}}{A PCA biplot figure of PC1 plotted against PC2.}
#'    \item{\code{biplot2}}{A PCA biplot figure of PC1 plotted against PC3. \code{NULL} if there are not at least 3 PCs returned by PCA.}
#'    \item{\code{biplot3}}{A PCA biplot figure of PC2 plotted against PC3. \code{NULL} if there are not at least 3 PCs returned by PCA.}
#' }
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
eda <- function(
  mat, 
  cor_method=c("pairwise.complete.obs", "pearson"), 
  cor_lbl = FALSE,
  impute = 0, 
  scale_flag = TRUE,
  pcs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
  colgroups = NULL,
  rowgroups = NULL,
  rowgroups_name = "cohort",
  frame = F) 
{
  
  #### CORRELATION MATRIX ####
  cor_mat <- mat %>% GGally::ggcorr(., method = cor_method, limits = F, label = cor_lbl, label_size = 3, label_alpha = T, size = 3, layout.exp = 1)
  
  #### HANDLING MISSINGNESS / SCALING ####
  mat[is.na(mat)] <- impute
  
  mat <- apply(mat, 2, scale, scale = scale_flag, center = F)
  
  #### PCA ####
  pca.ln <- prcomp(mat) 
  
  #### VARIANCE TABLE ####
  eigenvalues_ln <- matrix(pca.ln$sdev^2) # eigenvalues
  perc_variance <- round(100*matrix(pca.ln$sdev^2/sum(pca.ln$sdev^2)),1) # variance
  
  pca_summary <- data.frame("Principle Component" = 1:min(dim(mat)), "Eigenvalues" = eigenvalues_ln, "Percent Variance" = perc_variance, "Total Cumulative Variance" = purrr::accumulate(perc_variance, sum))
  
  var_tbl <- kableExtra::kbl(pca_summary, col.names = c("Principle Component", "Eigenvalues", "% Variance", "Total Cumulative Variance"), align = "c") %>% 
    kableExtra::kable_classic(full_width = F, html_font = "Cambria", position = "center") %>% 
    kableExtra::kable_styling(bootstrap_options = c("hover", "condensed"), fixed_thead = T)
  
  #### LOADINGS ####
  pca.ln.ld <- as.data.frame.matrix(pca.ln$rotation) # rotation is the loadings variable within the pca output.
  pca.ln.ld$chem <- row.names(pca.ln.ld)
  
  if (!is.null(colgroups)) {
    colgroups <- colgroups %>% dplyr::rename(chem = !!names(colgroups)[1])
  } else {
    colgroups <- data.frame(chem = colnames(mat), group = "1")
  }
  grouping <- names(colgroups)[2]
  
  plot_loadings_pca <- pca.ln.ld %>% 
    tidyr::gather(key = "PC", value = "Loading", -chem) %>% 
    tibble::as_tibble() %>% 
    dplyr::right_join(., colgroups, by = "chem")
  plot_loadings_pca$chem <- factor(as.character(plot_loadings_pca$chem), levels = unique(as.character(plot_loadings_pca$chem)))
  
  loadings <- plot_loadings_pca %>%
    dplyr::filter(PC %in% pcs) %>% 
    ggplot(aes(x = chem, y = Loading, color = !!sym(grouping))) + 
    geom_point() +
    geom_segment(aes(yend=0, xend = chem)) +
    facet_wrap(~ PC) + theme_bw() +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    geom_hline(yintercept = 0, size = 0.2) + ggtitle("Principle Component Loadings")
  
  #### BIPLOT ####
  if (is.null(rowgroups)) {
    rowgroups <- factor(rep(1, times = nrow(mat)))
  }
  
  bp.m <- mat %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(!!rowgroups_name := as.factor(rowgroups)) %>% 
    as.data.frame()
  
  biplot1 <- autoplot(pca.ln, data = bp.m, size = 0.8, colour = rowgroups_name, frame = frame, alpha = 0.5,
                     loadings = TRUE, loadings.colour = 'black',
                     loadings.label = TRUE, loadings.label.repel = T, 
                     loadings.label.size = 2.5, loadings.label.colour = 'black',
                     main = "Principal Component Analysis Biplot")
  
  if ("PC3" %in% plot_loadings_pca$PC) {
    biplot2 <- autoplot(pca.ln, x = 1, y = 3, data = bp.m, size = 0.8, colour = rowgroups_name, frame = frame, alpha = 0.5,
                        loadings = TRUE, loadings.colour = 'black',
                        loadings.label = TRUE, loadings.label.repel = T, 
                        loadings.label.size = 2.5, loadings.label.colour = 'black',
                        main = "Principal Component Analysis Biplot")
    
    biplot3 <- autoplot(pca.ln, x = 2, y = 3, data = bp.m, size = 0.8, colour = rowgroups_name, frame = frame, alpha = 0.5,
                        loadings = TRUE, loadings.colour = 'black',
                        loadings.label = TRUE, loadings.label.repel = T, 
                        loadings.label.size = 2.5, loadings.label.colour = 'black',
                        main = "Principal Component Analysis Biplot")
    
  } else {
    biplot2 <- NULL
    biplot3 <- NULL
  }
  
  list(cor = cor_mat, var = var_tbl, var_raw = pca_summary, load = loadings, biplot1 = biplot1, biplot2 = biplot2, biplot3 = biplot3)
}
