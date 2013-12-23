\name{clustcat-class}
\alias{clustcat-class}
\title{clustcat output objects}
\description{
  This class contains all the output of the clustercat
  function \describe{ \item{best_model}{A
  \code{\linkS4class{resultscat}} object containing the
  best model results according to the BIC criterion.}
  \item{models}{A list of \code{\linkS4class{resultscat}}
  object containing all the results sorted in ascending
  order according to the number of classes.} }
}
\examples{
  getSlots("clustcat")
}

