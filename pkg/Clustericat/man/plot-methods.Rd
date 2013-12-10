\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,clustcat-method}
\title{Plot function.}
\arguments{
  \item{x}{output object from \code{\linkS4class{clustcat}}.}

  \item{y}{Ignored}
}
\description{
  This function gives the graphical summary of the
  parameters of the best model estimated by the function
  \code{\linkS4class{resultscat}}.
  This plot helps the class interpretation for the best model.
   In the lines, the estimated classes are represented with respect to their proportions in decreasing order. 
   Their corresponding area depends on their proportion. 
   The cumulated proportions are indicated on the left side. 
   On the column, three indications are given. 
   The first one is the inter-variables correlations (rho_{kb}) for all the blocks of the class ordered by their strength of correlation (in decreasing order). 
   The second one is the intra-variables correlations (tau_{kb}) for each block drawn according to their strength dependencies (in decreasing order).
   The third is the variables repartition per blocks. A black cell indicates that the variable is affected to the block and a white cell indicates that, conditionally on this class, the variable is independent to the variables of this block.
}
\examples{
# Simple example with binary data
data("dentist")




# estimation of the model for a classes number equal to 1,2.
res <- clustercat(dentist, 1:2,modal=rep(2,5))


# a plot summarizing the best best model
plot(res)
}

