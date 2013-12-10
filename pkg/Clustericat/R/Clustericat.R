#' Clustericat Package  
#' 
#' This package performs clustering for categorical data-set.
#' 
#'
#'  
#' 
#' This package performs clustering for categorical data-set by the conditional correlated mixture model. It contains a function
#'    clusteringcat() which perform clustering on categorical data-set and returns
#'    object of appropriate class (refer to documentation of clustcat). The
#'    package also provide utility function like summary(),
#'    summary_dependencies() and plot to summarize results, to present the main
#'    conditional dependencies, to visualize the parameters respectively.
#' 
#' 
#' @examples 
#' # Simple example with binary data introduced by Goodman for illustrate the Latent Class model
#' #importation of the data
#' data("binaryexample")
#' 
#' #estimation of the model for a classes number equal to 1,2,3.
#' res <- clustercat(binaryexample, 1:3)
#' 
#' #presentation of the best model
#' summary(res)
#'
#' # Another example with simulated categorical data.
#' #importation of the data
#' data("simuldata")
#' 
#' #estimation of the model for a classes number equal to 1,2,3.
#' newstrategycat <- strategycat(simuldata, nb_init=3, stop_criterion=100)
#' res <- clustercat(simuldata, nb_cluster=2,strategy=newstrategycat)
#' 
#' #presentation of the best model
#' summary(res)
#' summary_dependencies(res)
#' 
#' 
#' @name Clustericat-package
#' @rdname Clustericat
#' 
NULL