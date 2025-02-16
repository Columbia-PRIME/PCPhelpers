% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsim_data.R
\name{rsim_data}
\alias{rsim_data}
\title{Simulate a data matrix of a given dimension and rank for PCP experiments.}
\usage{
rsim_data(
  sim_seed,
  nrow,
  ncol,
  rank,
  sigma,
  add_sparse = FALSE,
  nonneg = TRUE,
  ldist = runif,
  sdist = runif,
  zdist = rnorm,
  ...
)
}
\arguments{
\item{sim_seed}{Required seed to allow for reproducible results & distinct simulated matrices.}

\item{nrow}{Sets number of rows in the simulated data matrix.}

\item{ncol}{Sets the number of columns in the simulated data matrix.}

\item{rank}{Sets the rank of the simulated data matrix.}

\item{sigma}{Sets the standard deviation of the noise matrix, \code{Z}.}

\item{add_sparse}{Logical that sets whether sparse noise, \code{S}, is added or not. Default=\code{FALSE}.}

\item{nonneg}{Logical that sets whether the final M matrix should be non-negative. Default=\code{TRUE}.}

\item{ldist}{Function that defines a distribution to use when simulating the L matrix. Default=\code{runif}.}

\item{sdist}{Function that defines a distribution to use when simulating the S matrix. Default=\code{runif}.}

\item{zdist}{Function that defines a distribution to use when simulating the Z matrix. Default=\code{rnorm}.}
}
\value{
list containing the simulated \code{M}, \code{L}, \code{S} and \code{Z} matrices.
}
\description{
\code{rsim_data} simulates a data matrix, \code{M}, of a given dimension and rank by simulating
ground truth low rank (\code{L}), sparse (\code{S}) and noise (\code{Z}) components.
Meaning \code{M = L + S + Z}. The \code{M}, \code{L}, \code{S} and \code{Z} matrices are all accessible from this function,
allowing for easy evaluation of PCP's performance.
}
\examples{
data <- rsim_data(sim_seed = 1, nrow = 10, ncol = 10, rank = 3, sigma = .1, add_sparse = TRUE)
}
