\name{strategycat}
\alias{strategycat}
\title{strategycat function}
\usage{
  strategycat(data, nb_init = 5,
    stop_criterion = 20 * ncol(data),
    partition = partitioncat(data))
}
\arguments{
  \item{data}{Input data as matrix of non-zero integers.}

  \item{nb_init}{Integer of the number of times where a
  MCMC chain is trating. By default 5 MCMC chains are
  initializated.}

  \item{stop_criterion}{Integer corresponding to the number
  of successive iterations of the MCMC chain where if no
  best model are fund then the algorithm is stopped. By
  default it takes the values of 20x d.}

  \item{partition}{Vector of integer. The element j takes
  the values of the block where the variable j is affected
  for the initialization of the MCMC chains for all the
  classes. By default the function
  \code{\link{partitioncat}} is called.}
}
\value{
  Return an object of \code{\linkS4class{strategycat}}
  class.
}
\description{
  This function allows to create a new
  \code{\linkS4class{strategycat}} object. This function
  needs a data matrix as mandatory argument.
}
\examples{
# Simple example with binary data i
data("dentist")

#creation of a strategycat object with the constructor by default
defaultstrategycat <- strategycat(dentist)

#presentation of the adjustment parameters of the strategycat object
summary(defaultstrategycat)

#creation of new object strategycat
newstrategycat <- strategycat(dentist, nb_init=100, stop_criterion=10, partition=c(1,2,3,1,2))

#presentation of the new object
summary(newstrategycat)
}

