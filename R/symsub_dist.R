#' Symmetric Subspace Distance 
#' 
#' see: http://www.paper.edu.cn/scholar/showpdf/NUT2cN2INTz0YxeQh
#' (Wang et al., 2006) for more info
#' 
#' Because all solutions may not successfully identify the correct number of patterns, 
#' symmetric subspace distance compares matrices of different dimensions 
#' (e.g., true score matrix 1,000 x 4 vs. solution scores 1,000 x 5). 
#' Symmetric subspace distance defines a distance measure between two linear subspaces, 
#' where the number of individual vectors is constant but the number of elements they 
#' contain may differ. The symmetric distance between m-dimensional 
#' subspace U and n-dimensional subspace V is defined as 
#' \code{d(U, V) =sqrt{max (m, n)-sum_{i=1}^{m} sum_{j=1}^{n}left(u_{i}^{mathrm{T}} v_{j}right)^{2}}}
#' If two subspaces largely overlap, they will have a small distance; if they are almost orthogonal, 
#' they will have a large distance. To use this method to compare simulations with solutions, 
#' we took orthonormal bases of all included matrices and calculated the ratio between symmetric subspace
#' distance and \code{sqrt(max(m,n))}. Two subspaces are similar if and only if this ratio is \code{leq frac{1}{2}}.
#' @export
symsub_dist <- function(U, V) {
  
  # there is some hacky error-handling in here
  # to return `NA' instead of crashing
  if(any(is.na(U)) | any(is.na(V))) {return(NA)} else {
    
    U = as.matrix(U)
    V = as.matrix(V)
    
    # patterns should be columns
    if (nrow(U) < ncol(U)) {U <- t(U)}
    if (nrow(V) < ncol(V)) {V <- t(V)}
    
    # this gets an orthonormal basis
    # Both provide orthonormal bases, but svd does not give NaNs for lots of zeros
    # the distance is invariant to choice of orthonormal basis
    qrU <- qr.Q(qr(U)) # svd(U)$v # qr.Q(qr(U)) # svd(U)$u
    qrV <- qr.Q(qr(V)) # svd(V)$v # qr.Q(qr(V)) # svd(V)$u
    
    if (any(is.nan(qrU)) | any(is.nan(qrV))) {
      qrU <- svd(U)$u
      qrV <- svd(V)$u
    }
    
    m <- ncol(U)
    n <- ncol(V)
    
    # this is the formula
    # don't want the script to crash
    tryCatch(
      {dUV <- sqrt( min(m,n) - sum((t(qrU) %*% qrV)^2) ) # lawrence: changed max to min to see nesting properties.
      },
      error = function(er){
        dUV = NA
      },
      warning = function(cond){
        dUV = NA
      })
    
    # it happened somewhere that round((max(m,n) - sum((t(qrU) %*% qrV)^2)),10) was so close to 0
    # that the sqrt function gave NaN
    if(round((min(m,n) - sum((t(qrU) %*% qrV)^2)),10) == 0) {dUV = 0} # if sqrt(0), make zero # lawrence: changed max to min to see nesting properties.
    
    if (!is.na(dUV)) {ratio <- dUV/sqrt( min(m,n))} else {ratio = NA} # lawrence: changed max to min to see nesting properties.
    return(ratio)
  }
}
