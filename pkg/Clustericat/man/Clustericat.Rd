\name{Clustericat-package}
\alias{Clustericat-package}
\title{Clustericat Package}
\description{
  This package performs clustering for categorical
  data-set.
}
\details{


  This package performs clustering for categorical data-set
  by the conditional correlated mixture model. It contains
  a function clusteringcat() which perform clustering on
  categorical data-set and returns object of appropriate
  class (refer to documentation of clustcat). The package
  also provide utility function like summary(),
  summary_dependencies() and plot to summarize results, to
  present the main conditional dependencies, to visualize
  the parameters respectively.
}
\examples{
# Simple example with binary data
data("dentist")

# to define the parameters of the algorithm performing the estimation
st=strategycat(dentist,nb_init=35,stop_criterion=200)


# estimation of the model for a classes number equal to 1,2.
res <- clustercat(dentist, 1:2,modal=rep(2,5), strategy=st)

# presentation of the best model
summary(res)

# presentation of the parameters of the conditional dependencies for the best model
summary_dependencies(res)

# a plot summarizing the best best model
plot(res)
}

