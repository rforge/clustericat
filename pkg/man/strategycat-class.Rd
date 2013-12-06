\name{strategycat-class}
\alias{strategycat-class}
\title{strategycat input/output objects}
\description{
  This class contains the parameters estimated by the
  clustercat algorithm for a fixed number of class.
  \describe{ \item{nb_init}{Integer of the number of times
  where a Gibbs chain is trating. By default 5 Gibbs chains
  are initializated.} \item{stop_criterion}{Integer
  corresponding to the number of successive iterations of
  the Gibbs chain where if no best model are fund then the
  algorithm is stopped. By default it takes the values of
  20x d.} \item{partition}{Vector of integer. The element j
  takes the values of the block where the variable j is
  affected for the initialization of the Gibbs chains for
  all the classes. By default the function
  \code{\link{partitioncat}} is called.} }
}

