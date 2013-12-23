\docType{package}
\name{Clustericat-package}
\alias{Clustericat}
\alias{Clustericat-package}
\title{Clustericat Package}
\description{
  This package performs CLUSTERing for CATegorical data-set by using
  a mixture model taking into account the conditional dependencies between
  variables.
}
\details{


  This package performs clustering for categorical data-set
  by the conditional correlated mixture model. In this model, 
  variables are grouped into inter-independent and intra-dependent
  blocks in order to consider the main intra-class correlations. 
  The dependence between variables grouped into the same block of 
  a class is taken into account by mixing two extreme distributions, 
  which are respectively the independence and the maximum dependence ones.
  In the conditionally correlated data case, this approach is expected to 
  reduce biases involved by the latent class model and  to produce a
  meaningful dependency model with few additional parameters. 
  
  The parameters estimation by maximum likelihood is performed by an EM algorithm while a Gibbs algorithm is used for model selection to avoid combinatorial problems involved by the block structure search.
    
  Clustericat contains  a function clusteringcat() which perform clustering on
  categorical data-set and returns object of appropriate
  class (refer to documentation of \code{\linkS4class{clustcat}}). The package
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
\author{
  Matthieu Marbac and Christophe Biernacki and Vincent Vandewalle.
}
\references{
  Marbac M., Biernacki C., Vandewalle V., 2014. 
  "Model-based clustering for conditionally correlated categorical data". Rapport de recherche INRIA RR-8232.
}

