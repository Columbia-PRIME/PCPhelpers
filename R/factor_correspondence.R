#' Factor correspondence
#'
#' WITH error handling
#' if matrices are not same size, NA
#' If nn is false, we allow *signed permutations*
#'    find a matching between two sets of factors (presented as columns of
#'    matrices A and B) that minimizes the sum of the squared errors. 
#'
#'    Assumes that A and B have columns of unit L2 norm. Looks for a
#'    permutation matrix (or signed permutation matrix) Pi such that 
#'
#'       || A - B * Pi ||_F^2 
#'
#'    is minimized.
#'
#'    Uses an (exact) Lp relaxation, solved via the cvx package.
#'
#'   Inputs
#'       A, B matrices of the same size, with unit L2 norm columns
#'       nn   boolean, optional, default true. If nn, we search for a
#'       *permutation* matrix Pi. If nn is false, we allow *signed
#'       permutations* which can multiply the factors by -1. 
#'
#'   Outputs
#'       e -- the sum of squared errors under the best calibration Pi
#'       Pi -- the permutation or signed permutation matrix
#' @import CVXR
#' @export
factor_correspondence <- function (A, B, nn = TRUE) {
  # This is all from the CVXR package, which is kind of its own language
  # for convex optimization
  G <- t(B) %*% A
  n <- nrow(G)
  
  # Step 1. Define the variable to be estimated
  # Pi -- the permutation or signed permutation matrix
  Pi <- Variable(n,n)
  
  # Step 2. Define the objective to be optimized
  objective <- Maximize(base::sum(Pi * G))
  
  if (nn) {
    # Step 2.5. Subject to these constraints
    constX = list()
    for (i in 1:nrow(G)) {
      constX <-  c(constX, list(base::sum(Pi[,i]) == 1))
      constX <-  c(constX, list(base::sum(Pi[i,]) == 1))
    }
    constX <- c(constX, Pi >= 0)
  } else {
    # % allow sign flips 
    # Step 3. vector l1 norms along rows and columns
    constX = list()
    for (i in 1:nrow(G)) {
      constX <-  c(constX, list(base::sum(abs(Pi[,i])) <= 1))
      constX <-  c(constX, list(base::sum(abs(Pi[i,])) <= 1))
    }
  }  
  # Step 3. Create a problem to solve
  problem <- Problem(objective, constraints = constX)
  
  # Step 4. Solve it!
  result <- solve(problem)
  
  # Step 5. Extract solution and objective value
  perm <- round(result$getValue(Pi), 0)
  
  e <- norm(B,'f')^2 + norm(A,'f')^2 - 2 * base::sum(perm * G)
  # e -- the sum of squared errors under the best calibration \Pi
  
  # New matrix with best order
  newB <- B %*% perm
  
  return(list(rearranged = newB, perm_matrix = perm))
}
