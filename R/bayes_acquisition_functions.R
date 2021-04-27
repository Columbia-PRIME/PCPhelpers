#' @keywords internal
expected_improvement <- function(mu, sigma, y_best, points_to_eval) {

  acquisition <- purrr::map2_dbl(mu, sigma, function(m, s) {
  	if (s == 0) return(0)
    gamma <- (y_best - m) / s
    phi <- pnorm(gamma)  

    return(s * (gamma * phi + dnorm(gamma)))
  })
  
  points_to_eval[which.max(acquisition), ]

}

#' @keywords internal
probability_improvement <- function(mu, sigma, y_best, points_to_eval) {

  acquisition <- purrr::map2_dbl(mu, sigma, function(m, s) {
    if (s == 0) return(0)
    else return(pnorm((y_best - m) / s))
  })

  points_to_eval[which.max(acquisition), ]

}

#' @keywords internal
lower_confidence_bound <- function(mu, sigma, y_best, points_to_eval) {
  kappa <- 2 # tunable
  acquisition <- mu - kappa * sigma

  points_to_eval[which.min(acquisition), ]
}