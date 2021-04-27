#' Converts matrix to tibble dataframe
#'
#' \code{mat_to_df} returns the given matrix as a tibble dataframe.
#' As a tibble, the structure of the matrix can then be easily viewed with \code{ggplot}.
#'
#' @param mat The matrix to be converted to a tibble.
#' @param col_names vector of desired column names. By default = "chem01,
#'	 chem02, chem03, ..."
#' @return The matrix as a tibble object.
#' @examples
#' mat_as_df <- mat_to_df(matrix(rnorm(9), 3, 3))
#' @export
#' @importFrom magrittr %>%
mat_to_df <- function(mat, col_names=NA) {
  if (is.na(col_names)) {
    colnames(mat) = stringr::str_c("chem_", stringr::str_pad(1:ncol(mat), 2, pad = "0"))
  } else {
    colnames(mat) = col_names
  }
  ret =
    mat %>%
    tibble::as_tibble() %>%
    dplyr::mutate(obs_num = row_number()) %>%
    tidyr::pivot_longer(
      colnames(mat),
      names_to = "chem",
      values_to = "value"
    )

  ret

}
