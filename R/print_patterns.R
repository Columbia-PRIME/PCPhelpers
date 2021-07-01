#' Prints a matrix of patterns in a PCA-loadings-style plot.
#'
#' \code{print_patterns} is a (very) rough function to print a matrix of patterns in a PCA-style plot.
#'
#' @param pats A matrix of patterns where each pattern is in a different column. 
#' @param colgroups An optional dataframe with the column names for \code{pats} in the first column, and their grouping in the second column. By default, \code{colgroups = NULL}.
#' @param n The number of patterns to print. By default, \code{n = 6}.
#' @param pat_type An optional character specifying the prefix for pattern names. By default, \code{pat_type = "pat"}.
#' @param title The title for the plot. By default, \code{title = ""}.
#'
#' @return A ggplot object
#' @export
#' @importFrom magrittr %>%
#' @import ggplot2
print_patterns <- function(pats, colgroups = NULL, n = 1:6, pat_type = "pat", title = "") {
  
  if (!is.null(colgroups)) {
    colgroups <- colgroups %>% dplyr::rename(chem = !!names(colgroups)[1])
  } else {
    colgroups <- data.frame(chem = rownames(pats), group = "1")
  }
  
  if (n > ncol(pats)) n <- ncol(pats)
  
  grouping <- names(colgroups)[2]
  
  colnames(pats) <- paste0(pat_type, stringr::str_pad(1:ncol(pats), width = 2, pad = "0", side = "left"))
  
  pats.df <- pats %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(chem = colgroups[[1]]) %>%
    tidyr::pivot_longer(-chem, names_to = "pattern", values_to = "loading") %>%
    dplyr::right_join(., colgroups, by = "chem")

  pats.df$chem <- factor(as.character(pats.df$chem), levels = unique(as.character(pats.df$chem)))

  loadings <- pats.df %>%
    dplyr::filter(pattern %in% paste0(pat_type, stringr::str_pad(n, width = 2, pad = "0", side = "left"))) %>%
    ggplot(aes(x = chem, y = loading, color = !!sym(grouping))) +
    geom_point() +
    geom_segment(aes(yend=0, xend = chem)) +
    facet_wrap(~ pattern) + theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    geom_hline(yintercept = 0, size = 0.2) + ggtitle(title)
  
  loadings
}
