
# test_input function
# 
# This function test the validity of the input parameters for the function clustericat
#  
# @param data Input data as matrix of non-zero integers.
# @param nb_cluster Integer vector specifying the number of classes.
# @param strategycat An instance of the \code{\linkS4class{strategycat}} class which contains the adjustments parameters.
# @return a boolean if all is allright erreur=0 and erreur=1 otherwise
# 
# @name test_input
# @rdname test_input-methods

setGeneric (
    name= "test_input" ,
    def = function (data,nb_cluster,strategycat){ standardGeneric("test_input")}
) 

setMethod(
    f = "test_input" ,
    signature(data="matrix",nb_cluster="numeric",strategycat="strategycat") ,
    definition = function (data,nb_cluster,strategycat){
      erreur=0;
      if (any (data!=floor(data))){erreur=1;}
      if (any(data<=0)){erreur=1;}
      if (is.integer(nb_cluster)==FALSE){nb_cluster=as.integer(nb_cluster);}
      if (is.integer(strategycat@nb_init)==FALSE){strategycat@nb_init=as.integer(strategycat@nb_init);}
      if (is.integer(strategycat@stop_criterion)==FALSE){strategycat@stop_criterion=as.integer(strategycat@stop_criterion);}
      return(erreur)
    }
)



ordonne_variables <- function (data){
      maxi=rep(0,ncol(data));
      for (j in 1:ncol(data)){maxi[j]=max(data[,j])}
      data=data[,order(maxi,decreasing=TRUE)]
      return(data);
}



# exploite_resultats function
# 
# This function ordered the parameters
#  
# @param ret Input parameters of the MCMC algorithm
# @param data  Input data as matrix of non-zero integers.
# @return solution is the best parameters for a fix number of classes
#
#
#

filtre_resultats_si_1classe<-function(ret,data){
  if (ret[[6]]>1){
    ret[[4]]=as.matrix(ret[[4]])
    ret[[4]]=ret[[4]][ret[[5]]+1,];
    ret[[5]]=rep(1,nrow(data));
    for (k in 1:ret[[6]]){
      ret[[5]][which(rowSums(ret[[4]]<=ret[[4]][,k])==ncol(ret[[4]]))]=k;
    }
  }else{
    ret[[4]]=ret[[4]][ret[[5]]+1];
    ret[[5]]=rep(1,nrow(data))
  }
  return(ret)
}

##on renome les éléments de la liste
renomme<-function(ret){
  names(ret$parameters$tau)=paste("component",1:length(ret$parameters$tau),sep="")
  for (k in 1:length(ret$parameters$tau)){
    names(ret$parameters$tau[[k]])=paste("block",1:length(ret$parameters$tau[[k]]),sep="")   
  }
  names(ret$parameters$rho)=paste("component",1:length(ret$parameters$rho),sep="")
  names(ret$parameters$delta)=paste("component",1:length(ret$parameters$delta),sep="")
  for (k in 1:length(ret$parameters$delta)){
    names(ret$parameters$delta[[k]])=paste("block",1:length(ret$parameters$delta[[k]]),sep="")   
  }
  names(ret$parameters$alpha)=paste("component",1:length(ret$parameters$alpha),sep="")
  for (k in 1:length(ret$parameters$alpha)){
    names(ret$parameters$alpha[[k]])=paste("block",1:length(ret$parameters$alpha[[k]]),sep="") 	
  }
  return(ret);
}

# a modifier
filtre_identifiabilite<-function(ret,data){
  for (k in 1:ret$nbcluster){
    for (h in 1:length(ret$parameters$rho[[k]])){
      if (sum(ret$sigma[k,]==(h-1))<2){
        ret$parameters$rho[[k]][h]=0;
      }
      if ((sum(ret$sigma[k,]==(h-1))==2)&&(max(data[,which(ret$sigma[k,]==(h-1))[2]])==2)){
        #ret<-identifiabilise(ret,k,h)
      }
    }
  }
  return(ret)
}


setGeneric (
  name= "exploite_resultats" ,
  def = function (ret,data){ standardGeneric("exploite_resultats")}
) 

setMethod(
  f = "exploite_resultats" ,
  signature(ret="list",data="matrix") ,
  definition = function (ret,data){
    ret<-filtre_resultats_si_1classe(ret,data);
    ret<-renomme(ret)
    
    #on réordonne les classes en fonction des proportions
    passe<-ret;
    cl<-order(ret[["parameters"]][["proportions"]],decreasing=TRUE);
  	ret[["parameters"]][["proportions"]]=ret[["parameters"]][["proportions"]][cl];
  	
    for (loc in 1:length(cl)){
      	ret[["parameters"]][["rho"]][[loc]]<-passe[["parameters"]][["rho"]][[cl[loc]]];
      	ret[["parameters"]][["tau"]][[loc]]<-passe[["parameters"]][["tau"]][[cl[loc]]];
      	ret[["parameters"]][["delta"]][[loc]]<-passe[["parameters"]][["delta"]][[cl[loc]]];
      	ret[["parameters"]][["alpha"]][[loc]]<-passe[["parameters"]][["alpha"]][[cl[loc]]];
      	ret[["partition"]][which(passe[["partition"]]==cl[loc])]<-loc;
    }
    if (length(cl)>1){ret[["probapost"]]=ret[["probapost"]][,cl];}
    if (length(cl)>1){ret$sigma=ret$sigma[cl,];}
    rm(passe);
    
    ##on réodronne les blocks en fonction de leurs corrélations intra    
    noms<-colnames(data)
    ret<-filtre_identifiabilite(ret,data);
    solution=ret
    for (k in 1:ret$nbcluster){
      qui=1
      organise=order(ret[["parameters"]][["rho"]][[k]],decreasing=TRUE)
      for (loc in 1:length(organise)){
        solution$sigma[k,which(ret$sigma[k,]==(organise[loc]-1))]=loc
        if (sum(ret$sigma[k,]==(organise[loc]-1))>1){
          solution$parameters$alpha[[k]][[loc]]=matrix(0,length(which(ret$sigma[k,]==(organise[loc]-1))),length(ret$parameters$alpha[[k]][[organise[loc]]][[1]]));
          for (h in 1:length(which(ret$sigma[k,]==(organise[loc]-1)))){
            solution$parameters$alpha[[k]][[loc]][h,1:length(ret$parameters$alpha[[k]][[organise[loc]]][[h]])]=ret$parameters$alpha[[k]][[organise[loc]]][[h]]
          }
          rownames(solution$parameters$alpha[[k]][[loc]])=noms[which(ret$sigma[k,]==(organise[loc]-1))]
          solution$parameters$tau[[k]][[loc]]=ret$parameters$tau[[k]][[organise[loc]]]
         
          comment<-order(solution$parameters$tau[[k]][[loc]],decreasing=TRUE);
          solution$parameters$tau[[k]][[loc]]= solution$parameters$tau[[k]][[loc]][comment]
        
          solution$parameters$delta[[k]][[loc]]=ret$parameters$delta[[k]][[organise[loc]]]
          for (h in 1:length(solution$parameters$delta[[k]][[loc]])){
            solution$parameters$delta[[k]][[loc]][[h]]=ret$parameters$delta[[k]][[organise[loc]]][[comment[h]]]+1
          }
          
          #solution$parameters$delta[[k]][[loc]]=solution$parameters$delta[[k]][[loc]][comment,]
          
          solution$parameters$rho[[k]][loc]=ret$parameters$rho[[k]][organise[loc]]
          qui=qui+1
        }else{
          break;
        }
      }
      if (qui<=length(unique(ret$sigma[k,]))){
        ensemble=organise[qui:length(organise)]
        for (h in 1:length(ensemble)){if (any(ret$sigma[k,]==(ensemble[h]-1))){solution$sigma[k,which(ret$sigma[k,]==(ensemble[h]-1))]=qui}}
        solution$parameters$tau[[k]][[qui]]=0
        solution$parameters$delta[[k]][[qui]]=0
        solution$parameters$rho[[k]][qui]=0
        colone=qui
        while((all(ret$sigma[k,]!=(organise[colone]-1)))&&(colone<=length(organise))){colone=colone+1;}
        if (colone<=length(organise)){
          solution$parameters$alpha[[k]][[qui]]=matrix(0,length(which(solution$sigma[k,]==qui)),max(data))#length(ret$parameters$alpha[[k]][[organise[colone]]][[1]]))
          for (j in 1:nrow(solution$parameters$alpha[[k]][[qui]])){
            solution$parameters$alpha[[k]][[qui]][j,1:length(ret$parameters$alpha[[k]][[(ret$sigma[k,which(solution$sigma[k,]==qui)])[j]+1]][[1]])]=ret$parameters$alpha[[k]][[(ret$sigma[k,which(solution$sigma[k,]==qui)])[j]+1]][[1]]
          }
          length(solution$parameters$rho[[k]])=qui
          length(solution$parameters$delta[[k]])=qui
          length(solution$parameters$tau[[k]])=qui
          length(solution$parameters$alpha[[k]])=qui
          rownames(solution$parameters$alpha[[k]][[qui]])=noms[which(solution$sigma[k,]==qui)]
        }
      }else{
        length(solution$parameters$rho[[k]])=qui-1
        length(solution$parameters$delta[[k]])=qui-1 
        length(solution$parameters$tau[[k]])=qui-1
        length(solution$parameters$alpha[[k]])=qui-1
      }
    }
    colnames(solution$sigma)=colnames(data);
    return(solution)
  }
)






extrait_loi_jointe<-function(alpha,rho,delta,tau){
  loi_jointe<-((alpha[[1]])%*%t(alpha[[2]])) * (1-rho)
  for (loc in 1:length(delta)){
    loi_jointe[delta[[loc]][1]+1, delta[[loc]][2]+1]<-loi_jointe[delta[[loc]][1]+1, delta[[loc]][2]+1] + rho*tau[loc]
  }
  return(loi_jointe);
}

loi_binaire<-function(loi_jointe,delta){
  binaire=matrix(0,2,2);
  for (loc in 1:length(delta)){
    binaire[delta[[loc]][2]+1,delta[[loc]][2]+1]=binaire[delta[[loc]][2]+1,delta[[loc]][2]+1]+loi_jointe[delta[[loc]][1]+1,delta[[loc]][2]+1]
    binaire[2-delta[[loc]][2],delta[[loc]][2]+1]=binaire[2-delta[[loc]][2],delta[[loc]][2]+1]+loi_jointe[delta[[loc]][1]+1,2-delta[[loc]][2]]
  }
  return(binaire)
}

second.degre<-function(a,b,c){
  delta=b^2-4*a*c
  if (delta>=0){
    sol=c((-b-sqrt(delta))/(2*a),(-b+sqrt(delta))/(2*a))
  }else{
    sol= -10
  }
  return(sol)
}

check_new_rho<-function(loi_jointe,binaire,param){
  a2=second.degre(1,binaire[2,1]/(1-param$rho)-1-binaire[1,2]/(1-param$rho),binaire[1,2]/(1-param$rho))
  param$ok=0;
  param$alpha[[2]]=c(a2[1],1-a2[1])
  for (loc in 1:length(param$delta)){
    if (param$delta[[loc]][2]==1){
      param$alpha[[1]][loc]=loi_jointe[loc,1]/((1-param$rho)*param$alpha[[2]][1])
    }else{
      param$alpha[[1]][loc]=loi_jointe[loc,2]/((1-param$rho)*param$alpha[[2]][2])
    }
    if (param$delta[[loc]][2]==1){
      param$tau[loc]=(loi_jointe[loc,2]-((1-param$rho)*param$alpha[[2]][2]*param$alpha[[1]][loc]))/param$rho
    }else{
      param$tau[loc]=(loi_jointe[loc,1]-((1-param$rho)*param$alpha[[2]][1]*param$alpha[[1]][loc]))/param$rho
    }
  }
  test<-extrait_loi_jointe(param$alpha,param$rho,param$delta,param$tau)
  #if (sum(abs(test - loi_jointe))>0.001){print("pb")}
  if ( (all(is.na(c(param$alpha[[1]],param$alpha[[2]],param$tau,param$rho))==FALSE)) && (all(c(param$alpha[[1]],param$alpha[[2]],param$tau,param$rho)>=0)) && (all(c(param$alpha[[1]],param$alpha[[2]],param$tau,param$rho)<=1))){param$ok=1}
  if ((param$ok==0)&&(length(a2)==2)){
    param$alpha[[2]]=c(a2[2],1-a2[2])
    for (loc in 1:length(param$delta)){
      if (param$delta[[loc]][2]==1){
        param$alpha[[1]][loc]=loi_jointe[loc,1]/((1-param$rho)*param$alpha[[2]][1])
      }else{
        param$alpha[[1]][loc]=loi_jointe[loc,2]/((1-param$rho)*param$alpha[[2]][2])
      }
      if (param$delta[[loc]][2]==1){
        param$tau[loc]=(loi_jointe[loc,2]-((1-param$rho)*param$alpha[[2]][2]*param$alpha[[1]][loc]))/param$rho
      }else{
        param$tau[loc]=(loi_jointe[loc,1]-((1-param$rho)*param$alpha[[2]][1]*param$alpha[[1]][loc]))/param$rho
      }
    }
    test<-extrait_loi_jointe(param$alpha,param$rho,param$delta,param$tau)
    #if (sum(abs(test - loi_jointe))>0.001){print("pb")}
 	 if ( (all(is.na(c(param$alpha[[1]],param$alpha[[2]],param$tau,param$rho))==FALSE)) && (all(c(param$alpha[[1]],param$alpha[[2]],param$tau,param$rho)>=0)) && (all(c(param$alpha[[1]],param$alpha[[2]],param$tau,param$rho)<=1))){param$ok=1}
  }
  return(param)
}

estim_param_min<-function(loi_jointe,binaire,param){
  propose=param
  param$ok=1;
  propose$rho=param$rho/2;
  propose$ok=0;
  
  while((abs(param$rho- propose$rho)>10^(-6))||(propose$ok==0)){
   # print(c(propose$rho,param$rho))
    propose<-check_new_rho(loi_jointe,binaire,propose)
    if (propose$ok==1){
      sauv=propose
      propose$rho=propose$rho-abs(propose$rho-param$rho)/2
      param=sauv;
    }else{
      propose$rho=(propose$rho+param$rho)/2
    }
  }
  return(propose)
}

estim_param_max<-function(loi_jointe,binaire,param){
  propose=param
  param$ok=1;
  propose$rho=(1+param$rho)/2;
  propose$ok=0;
  
  while((abs(param$rho- propose$rho)>10^(-6))||(propose$ok==0)){
   # print(c(propose$rho,param$rho))
    propose<-check_new_rho(loi_jointe,binaire,propose)
    if (propose$ok==1){
      sauv=propose
      propose$rho=propose$rho+abs(propose$rho-param$rho)/2
      param=sauv;
    }else{
      propose$rho=propose$rho-abs(propose$rho-param$rho)/2
    }
  }
  return(propose)
}

copie<-function(ret,k,h,param){
  ret$parameters$alpha[[k]][[h]]=param$alpha
  ret$parameters$rho[[k]][h]=param$rho
  ret$parameters$tau[[k]][[h]]=param$tau
  return(ret)
}


identifiabilise<-function(ret,k,h){
  loi_jointe<-extrait_loi_jointe(ret$parameters$alpha[[k]][[h]],ret$parameters$rho[[k]][h],ret$parameters$delta[[k]][[h]],ret$parameters$tau[[k]][[h]])
  binaire<-loi_binaire(loi_jointe,ret$parameters$delta[[k]][[h]])
  param<-list(alpha=ret$parameters$alpha[[k]][[h]],rho=ret$parameters$rho[[k]][h],delta=ret$parameters$delta[[k]][[h]],tau=ret$parameters$tau[[k]][[h]])
  #param_min<-estim_param_min(loi_jointe,binaire,param);
  param_max<-estim_param_max(loi_jointe,binaire,param);
  #if ((1-param_max$rho) < param_min$rho){
    ret=copie(ret,k,h,param_max)
  #}else{
  #  ret=copie(ret,k,h,param_min)
  #}
  return(ret)
}
