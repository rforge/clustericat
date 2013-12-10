\name{partitioncat}
\alias{partitioncat}
\title{partitioncat function}
\usage{
  partitioncat(data)
}
\arguments{
  \item{data}{Input data as matrix of non-zero integers.}
}
\value{
  Return an vector of the variables partition.
}
\description{
  This function proposed a partition of the variables for
  the Gibbs chain initialization. The partition is done by
  a CAH algorithm based on the Cramer's V distance matrix
  between the couples of variables. The partition with at
  the maximum four variables in the same block which
  minimizes the number of blocks is return by this
  function.
}
\examples{
# Simple example with binary data 
data("dentist")
partitioncat(dentist)
}

