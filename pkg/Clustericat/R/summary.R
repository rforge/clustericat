#'
#' Summary function.
#' 
#' This function gives the summary of output from \code{clustercat} or \code{strategycat}.
#' For a  \code{clustercat} object, this function presents the repartition per block of the variables and the correlation coefficient Rho.
#'
#'
#' @param object output object from \code{\link{clustercat}} or \code{\linkS4class{strategycat}}.
#' 
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' 
#' 
NULL

#' @rdname summary-methods
#' @aliases summary summary,clustcat-method
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
#' 
#' #creation of a strategycat object with the constructor by default
#' defaultstrategycat <- strategycat(binaryexample)
#' 
#' #presentation of the adjustment parameters of the strategycat object
#' summary(defaultstrategycat)


setMethod(
		f="summary",
		signature = c("clustcat"),
		definition = function(object,...) {
				cat("**************************************************************************************\n")
				cat("Number of classes:",object@best_model@nbcluster,"      BIC value:",object@best_model@bic,"      log-Likelihood value:",object@best_model@likelihood," \n")
        
				cat("*************************************\n")
				cat("Proportions:",object@best_model@parameters@proportions,"\n");
				cat("\n")
				cat("*************************************\n")
				
				for (k in 1:object@best_model@nbcluster){
					cat("Blocks repartition of the variables for the class ", k,":\n",sep="")
					valeur=data.frame(Variables=rep(NA,length(object@best_model@parameters@alpha[[k]])),Rho=rep(0,length(object@best_model@parameters@alpha[[k]])))
					rownames(valeur)=paste("Block ",1:length(object@best_model@parameters@alpha[[k]])," ",sep="")
					for (b in 1:length(object@best_model@parameters@alpha[[k]])){
					  valeur[b,1]=as.character(toString(row.names(object@best_model@parameters@alpha[[k]][[b]])))
					  valeur[b,2]=object@best_model@parameters@rho[[k]][[b]]
					}
					print(valeur)
					cat("\n")
				}
				cat("\n**************************************************************************************\n")
				
		}
)


#' @rdname summary-methods
#' @aliases summary summary,strategycat-method

setMethod(
  f="summary",
  signature = c("strategycat"),
  definition = function(object,...) {
    cat("**************************************************************************************\n")
    cat("Number of the Gibbs chains:",object@nb_init,"\n")
    cat("Stopping criterion (q_max):",object@stop_criterion,"\n")
    cat("Initial partition of variables:",object@sigma,"\n")
    cat("\n**************************************************************************************\n")
    
  }
)
