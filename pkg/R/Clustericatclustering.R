#'
#' @include optionclasses.R
#' @include tools.R
#' 
NULL 



#' clustercat function
#' 
#' This function performs clustering for categorical data using the block model extension of the latent class model.
#'  
#' 
#' 
#' @param data Input data as matrix of non-zero integers.
#' @param nb_cluster Integer vector specifying the number of classes.
#' @param strategy An instance of the \code{\linkS4class{strategycat}} class which contains the adjustments parameters.
#' @return Return an instance of the \code{\linkS4class{clustcat}} class. Those two attributes will contains all outputs:
#' 
#' 
#' @examples
#'  
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
#' nvstrategycat <- strategycat(simuldata, nb_init=3, stop_criterion=100)
#' res <- clustercat(simuldata, nb_cluster=2, strategy=nvstrategycat)
#' 
#' #presentation of the best model
#' summary(res)
#' summary_dependencies(res)
#'
#' @export
#' 
#' @name clustercat
#' @useDynLib Clustericat
#'

clustercat <- function(data,nb_cluster,modal=0,strategy=strategycat(data)){
  if (is.matrix(data)==FALSE){data<-as.matrix(data);}
 
  codage<-list();
  nom<-colnames(data)
  if (is.character(data)){
  	passe<-matrix(0,nrow(data),ncol(data))
	for (j in 1:ncol(data)){
	  codage[[j]]<- sort(unique(as.factor(data[,j])))
	  passe[,j]<- as.numeric(as.factor(data[,j]))
	}
	rm(data)
    data<-as.matrix(passe)
    colnames(data)<-nom 
  }
  if (is.null(colnames(data))){
	  colnames(data)<-paste("X",1:ncol(data)," ",sep="")
	}
  if (is.numeric(data)){
      erreur=test_input(data,nb_cluster,strategy)
	  if (erreur==0){
	  		if (any(modal==0)){
	  			modal=rep(0,ncol(data))
	  		 	for (j in 1:ncol(data)){modal[j] <- max(data[,j])}
	  		 }
	      data<-ordonne_variables(data);
	      valeur_bic_models<-length(nb_cluster)
	      reponse=list();
	      differents=list();
	      for (nbcl in 1:length(nb_cluster)){
	    	  ret = .Call( "BMEclustering", data , modal,nb_cluster[nbcl], strategy@sigma, strategy@nb_init, strategy@stop_criterion,PACKAGE = "Clustericat" )
	          solution<-exploite_resultats(ret,data);
	    	  param_temp=new("parameterscat",rho=solution$parameters$rho,tau=solution$parameters$tau,delta=solution$parameters$delta,alpha=solution$parameters$alpha,proportions=solution$parameters$proportions)
	    	  if (length(codage)>0){
	    	  	for (k in 1:nb_cluster[nbcl]){
	    	  		for (b in 1:length(param_temp@delta[[k]])){
	    	  			qui<-which( solution$sigma[k,]==b);
	    	  			if ((length(qui)>1)&&(param_temp@rho[[k]][b]>0)){	    	  			
	    	  				for (h in 1:length(param_temp@delta[[k]][[b]])){
	    	  						passe<-param_temp@delta[[k]][[b]][[h]]
	    	  						param_temp@delta[[k]][[b]][[h]]<-as.character(rep(NA,length(passe)))
	    	  					for (loc in 1:length(passe)){
	    	  						param_temp@delta[[k]][[b]][[h]][loc]<- paste(codage[[qui[loc]]][passe[loc]])
	    	  					}
	    	  				}
	    	  			}
	    	  		}
	    	  	}
	    	  }
	    	  differents[[nbcl]]=new("resultscat",sigma=solution$sigma,bic=solution$bic,likelihood=solution$likelihood,probapost=as.matrix(solution$probapost),nbcluster=solution$nbcluster,partition=solution$partition,parameters=param_temp);
	    	  valeur_bic_models[nbcl]=solution$bic
	      }
	      resultats=new("clustcat",best_model=differents[[which(valeur_bic_models==max(valeur_bic_models))[1]]],models=differents)
	      return(resultats);
  	  }
  }else{
  	print("The format of the data set is not allowed");
  }
}
