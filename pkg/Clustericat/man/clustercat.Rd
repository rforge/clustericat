\name{clustercat}
\alias{clustercat}
\title{clustercat function}
\usage{
  clustercat(data, nb_cluster,modal=0,
    strategy = strategycat(data))
}
\arguments{
  \item{data}{Input data as matrix of non-zero integers.}

  \item{nb_cluster}{Integer vector specifying the number of
  classes.}
  
    \item{modal}{Vector of modalities. If modal=0, then the modalities number of a variable is equal to the number of the different observed levels.}

  \item{strategy}{An instance of the
  \code{\linkS4class{strategycat}} class which contains the
  adjustments parameters.}
}
\value{
  Return an instance of the \code{\linkS4class{clustcat}}
  class. Those two attributes will contains all outputs:
}
\description{
  This function performs clustering for categorical data
  using the block model extension of the latent class
  model.
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

