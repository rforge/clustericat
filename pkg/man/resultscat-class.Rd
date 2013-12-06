\name{resultscat-class}
\alias{resultscat-class}
\title{resultscat output objects}
\description{
  This class contains the model and the parameters
  estimated by the clustercat algorithm for a fixed number
  of class. \describe{ \item{sigma}{This matrix indicates
  the block repartition of the variables. sigma(k,j)=b if
  the variable j is affected to the block b for the
  component k.} \item{bic}{Value of the BIC criterion.}
  \item{likelihood}{Value of the log-likelihood.}
  \item{nbcluster}{Number of cluster.}
  \item{probapost}{This matrix indicates probabilities for
  each individuals to be arisen from each class.}
  \item{partition}{Partition estimated using the MAP.}
  \item{parameters}{A \code{\linkS4class{parameterscat}}
  object containing the estimated parameters.} }
}

