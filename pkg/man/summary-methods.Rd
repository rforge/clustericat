\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,clustcat-method}
\alias{summary,strategycat-method}
\title{Summary function.}
\arguments{
  \item{object}{output object from \code{\link{clustercat}}
  or \code{\linkS4class{strategycat}}.}
}
\description{
  This function gives the summary of output from
  \code{clustercat} or \code{strategycat}. For a
  \code{clustercat} object, this function presents the
  repartition per block of the variables and the
  correlation coefficient Rho.
}
\examples{
# Simple example with binary data
data("dentist")




# estimation of the model for a classes number equal to 1,2.
res <- clustercat(dentist, 1:2,modal=rep(2,5))

# presentation of the best model
summary(res)
}

