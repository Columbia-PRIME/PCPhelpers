#' @keywords internal
eval_params2 <- function(
  settings, 
  test_mat,
  pcp_mat, 
  test_mask
) {
  if ("message" == names(pcp_mat)[1]) {
    stats <- cbind(settings, data.frame(rel_err = NA, L_rank = NA, S_sparsity = NA, iterations = NA, run_error = pcp_mat$message))
  } else if ("L_list" == names(pcp_mat)[3]) {
    if ("r" %in% names(settings)) settings <- settings[setdiff(names(settings), "r")]
    test_mat[is.na(test_mat)] <- 0
    stats <- purrr::imap_dfr(.x = 1:length(pcp_mat$L_list), ~ data.frame(r = .y, settings, 
      rel_err = norm((test_mat - pcp_mat$L_list[[.x]])*test_mask, "F") / norm(test_mat*test_mask, "F"), 
      L_rank = Matrix::rankMatrix(pcp_mat$L_list[[.x]], tol = 1e-04),
      S_sparsity = sparsity(pcp_mat$S_list[[.x]], tol = 1e-04),
      iterations = NA, run_error = NA
    ))
  } else {
    test_mat[is.na(test_mat)] <- 0
    stats <- settings
    stats$rel_err <- norm((test_mat - pcp_mat$L)*test_mask, "F") / norm(test_mat*test_mask, "F")
    stats$L_rank <- Matrix::rankMatrix(pcp_mat$L, tol = 1e-04)
    stats$S_sparsity <- sparsity(pcp_mat$S, tol = 1e-04)
    stats$iterations <- pcp_mat$final_iter
    stats$run_error <- NA
  }
  return(stats)
}
