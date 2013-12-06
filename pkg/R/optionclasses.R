#' parameterscat output objects 
#' 
#' This class contains the parameters estimated by the clustercat algorithm for a fixed number of class.
#' \describe{
#' \item{rho}{List where rho[[k]][[b]] indicated the proportion of the maximum dependence distribution in the block b of the class k.}
#' \item{tau}{List where tau[[k]][[b]] indicated parameters of the first block variable multinomial of the maximum dependence distribution in the block b of the class k.}
#' \item{delta}{List where delta[[k]][[b]] indicated modalities relations involving by the maximum dependence distribution in the block b of the class k.}
#' \item{alpha}{List where rho[[k]][[b]] indicated the parameters of the independence distribution in the block b of the class k.}
#' \item{proportions}{Vector of classes proportions.}
#' }
#' 
#' @name parameterscat-class
#' @exportClass parameterscat
#' 
#'
#'
#'

setClass(
	Class="parameterscat",
	representation=representation(
		rho="list",
		tau="list",
		delta="list",
		alpha="list",
		proportions="numeric"
	),
	prototype = prototype(
		rho=list(),
		tau=list(),
		delta=list(),
		alpha=list(),
		proportions=numeric(0)
	)
)


#'
#' resultscat output objects
#' 
#' This class contains the model and the parameters estimated by the clustercat algorithm for a fixed number of class.
#' \describe{
#' \item{sigma}{This matrix indicates the block repartition of the variables. sigma(k,j)=b if the variable j is affected to the block b for the component k.}
#' \item{bic}{Value of the BIC criterion.}
#' \item{likelihood}{Value of the log-likelihood.}
#' \item{nbcluster}{Number of cluster.}
#' \item{probapost}{This matrix indicates probabilities for each individuals to be arisen from each class.}
#' \item{partition}{Partition estimated using the MAP.}
#' \item{parameters}{A \code{\linkS4class{parameterscat}} object containing the estimated parameters.} 
#' }
#' 
#' @name resultscat-class
#' @exportClass resultscat
#' 
#'
#'

setClass(
	Class="resultscat",
	representation=representation(
		sigma="matrix",
		bic="numeric",
		likelihood="numeric",
		probapost="matrix",
		nbcluster="numeric",
		partition="numeric",
		parameters="parameterscat"
	),
	prototype = prototype(
		sigma=matrix(nrow=0,ncol=0),
		bic=numeric(0),
		likelihood=numeric(0),
		probapost=matrix(nrow=0,ncol=0),
		partition=numeric(0),
		nbcluster=numeric(0)
	)
)

#'
#' clustcat output objects
#' 
#' This class contains all the output of the clustercat function
#' \describe{
#' \item{best_model}{A \code{\linkS4class{resultscat}} object containing the best model results according to the BIC criterion.}
#' \item{models}{A list of \code{\linkS4class{resultscat}} object containing all the results sorted in ascending order according to the number of classes.}
#' }
#' @name clustcat-class
#' @exportClass clustcat
#' 
#'
#'

setClass(
		Class = "clustcat",
		representation = representation(
				best_model = "resultscat",
				models="list"
		),
		prototype = prototype(
				models=list()
		)
)

#' strategycat input/output objects 
#' 
#' This class contains the parameters estimated by the clustercat algorithm for a fixed number of class.
#' \describe{
#' \item{nb_init}{Integer of the number of times where a Gibbs chain is trating. By default 5 Gibbs chains are initializated.}
#' \item{stop_criterion}{Integer corresponding to the number of successive iterations of the Gibbs chain where if no best model are fund then the algorithm is stopped. By default it takes the values of 20x d.}
#' \item{partition}{Vector of integer. The element j takes the values of the block where the variable j is affected for the initialization of the Gibbs chains for all the classes. By default the function \code{\link{partitioncat}} is called.}
#' }
#' 
#' @name strategycat-class
#' @exportClass strategycat
#' 
#'
#'
#'

setClass(
  Class = "strategycat",
  representation = representation(
    nb_init = "numeric",
    stop_criterion="numeric",
    sigma="numeric"
  ),
  prototype = prototype(
    nb_init = numeric(0),
    stop_criterion = numeric(0),
    sigma=numeric(0)
  )
)




#' partitioncat function
#' 
#' This function proposed a partition of the variables for the Gibbs chain initialization. The partition is done by a CAH algorithm based on the Cramer's V distance matrix between the couples of variables.
#' The partition with at the maximum four variables in the same block which minimizes the number of blocks is return by this function.
#'  
#' 
#' 
#' @param data Input data as matrix of non-zero integers.
#' @return Return an vector of the variables partition.
#' 
#' 
#' @examples
#'  
#' # Simple example with binary data introduced by Goodman for illustrate the Latent Class model
#' #importation of the data
#' data("simuldata")
#' partitioncat(simuldata)
#' @export

partitioncat<-function(data){
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
  if (is.numeric(data)){
    dist<-matrix(NA,ncol(data),ncol(data));
    for (j in 1:(ncol(data)-1)){
      for (j2 in (j+1):ncol(data)){
      	distance <- table(data[,j],data[,j2])
      	distance <- distance/sum(distance)
     	dist[j,j2]<- sum((distance - (as.matrix(rowSums(distance)))%*%t(as.matrix(colSums(distance))))^2 / ( (as.matrix(rowSums(distance)))%*%t(as.matrix(colSums(distance))))) / (min(dim(distance))-1)
      
      #  dist[j,j2]<- chisq.test(table(data[,j],data[,j2]),correct=FALSE)$stat/(nrow(data))
        dist[j2,j]<-dist[j,j2];
      }
      dist[j,j]<-1;
    }
    dist[ncol(data),ncol(data)]<-1;
    dist<-(1-dist)*100
    distance<-as.dist(dist)
    test<-hclust(distance,method="ward");
    for (model in 1:ncol(data)){partition<-cutree(test,k=model);if(all(table(partition)<5)){break;}}
  }else{
    print("The format of the data set is not allowed");
    partition=NULL;
  }
  return(partition)
}




#' strategycat function
#' 
#' This function allows to create a new \code{\linkS4class{strategycat}} object. This function needs a data matrix as mandatory argument.
#'  
#' 
#' 
#' @param data Input data as matrix of non-zero integers.
#' @param nb_init Integer of the number of times where a MCMC chain is trating. By default 5 MCMC chains are initializated.
#' @param stop_criterion Integer corresponding to the number of successive iterations of the MCMC chain where if no best model are fund then the algorithm is stopped. By default it takes the values of 20x d.
#' @param partition Vector of integer. The element j takes the values of the block where the variable j is affected for the initialization of the MCMC chains for all the classes. By default the function \code{\link{partitioncat}} is called.

#' @return Return an object of \code{\linkS4class{strategycat}} class.
#' 
#' 
#' @examples
#'  
#' # Simple example with binary data introduced by Goodman for illustrate the Latent Class model
#' #importation of the data
#' data("simuldata")
#' 
#' #creation of a strategycat object with the constructor by default
#' defaultstrategycat <- strategycat(simuldata)
#' 
#' #presentation of the adjustment parameters of the strategycat object
#' summary(defaultstrategycat)
#' 
#' #creation of new object strategycat
#' newstrategycat <- strategycat(simuldata, nb_init=100, stop_criterion=10, partition=c(1,2,3,1))
#' 
#' #presentation of the new object
#' summary(newstrategycat)
#'@export

strategycat<-function(data,nb_init=5,stop_criterion=20*ncol(data),partition=partitioncat(data)){
 if (is.matrix(data)==FALSE){data<-as.matrix(data);}
  if (is.numeric(data) || is.character(data)){
    param_temp=new("strategycat",nb_init=nb_init,stop_criterion=stop_criterion,sigma=partition);
  }else{
    print("The format of the data set is not allowed");
    param_temp=NULL;
  }
  return(param_temp);
}

#'
#'
NULL