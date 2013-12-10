\name{summary_dependencies}
\alias{summary_dependencies}
\title{Summary Dependencies function.}
\usage{
  summary_dependencies(object)
}
\arguments{
  \item{object}{output object from
  \code{\link{clustercat}}.}
}
\description{
  This function gives the summary of the dependencies of
  the best model from a \code{clustcat} object. For each
  block, the first line presents the variables of the block
  and the correlation coefficient Rho. The other lines
  print the value of the parameters of maximum dependence
  distribution (Tau: probabilities of the modalities
  crossing, Delta: modalities crossing).
}
\examples{
# Simple example with binary data
data("dentist")


# estimation of the model for a classes number equal to 1,2.
res <- clustercat(dentist, 1:2,modal=rep(2,5))

# presentation of the best model
summary(res)

# presentation of the parameters of the conditional dependencies for the best model
summary_dependencies(res)

}

