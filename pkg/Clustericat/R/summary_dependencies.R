#'
#' Summary Dependencies function.
#' 
#' This function gives the summary of the dependencies of the best model from a \code{clustcat} object.
#' For each block, the first line presents the variables of the block and the correlation coefficient Rho. The other lines
#' print the value of the parameters of maximum dependence distribution (Tau: probabilities of the modalities crossing, Delta: modalities crossing).
#' 
#' @param object output object from \code{\link{clustercat}}.
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
#' summary_dependencies(res)
#'
#' # Another example with simulated categorical data.
#' #importation of the data
#' data("simuldata")
#' 
#' #estimation of the model for a classes number equal to 1,2,3.
#' res <- clustercat(simuldata, nb_cluster=2)
#' 
#' #presentation of the best model
#' summary_dependencies(res)
#' 
#' 


#' @export
#' 
#' @exportPattern "^[[:alpha:]]+"
#' @useDynLib Clustericat
#' 
#'
#'
summary_dependencies<-function(object){
	cat("**************************************************************************************\n")
	cat("**************************************************************************************\n")
	for (k in 1:object@best_model@nbcluster){
		cat("Blocks repartition of the variables for the class ", k,":\n",sep="")
		cat("\n")
		for (b in 1:length(object@best_model@parameters@alpha[[k]])){
			cat("Block ",b," contains the variables: ",row.names(object@best_model@parameters@alpha[[k]][[b]])," with Rho=",object@best_model@parameters@rho[[k]][[b]],"\n");		
			cat("\n")
			if (object@best_model@parameters@rho[[k]][[b]]>0){
				result=matrix(NA,length(object@best_model@parameters@tau[[k]][[b]]),length(which(object@best_model@sigma[k,]==b))+1)
				colnames(result)=c("Tau",row.names(object@best_model@parameters@alpha[[k]][[b]]))
				result[,1]=ceiling(object@best_model@parameters@tau[[k]][[b]]*1000000)/1000000
				for (h in 1:length(object@best_model@parameters@tau[[k]][[b]])){
				  result[h,-1]=object@best_model@parameters@delta[[k]][[b]][[h]]
				}
				result=result[order(result[,1],decreasing=TRUE),]
				print(result)				
				cat("\n")
				cat("\n")
			}
		}
		cat("\n**************************************************************************************\n")
	}
	cat("**************************************************************************************\n")

}

