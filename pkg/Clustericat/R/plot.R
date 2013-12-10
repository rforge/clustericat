#'
#' Plot function.
#' 
#' This function gives the graphical summary of the parameters of the best model of an object of code{clustcat} S4 class.
#' 
#' @param x output object from \code{\link{clustcat}}.
#' @param y Ignored
#' 
#' @importFrom graphics plot
#' @name plot
#' @rdname plot-methods
#' @docType methods
#' @exportMethod plot
#' 
#' 
NULL


#' @rdname plot-methods
#' @aliases plot plot,plot clustcat-method
#' 

#' @examples
#' # Simple example with binary data introduced by Goodman for illustrate the Latent Class model
#' #importation of the data
#' data("simuldata")
#' 
#' #estimation of the model for a classes number equal to 1,2,3.
#' res <- clustercat(simuldata, 2)
#' 
#' #presentation of the best model
#' plot(res)
#'
#' 



setMethod(
		f="plot",
		signature = c("clustcat"),
		definition = function(x,y,...) {
			cumpi=c(0,cumsum(x@best_model@parameters@proportions))
		    cumpi=ceiling(cumpi*100)/100
		    par(mar=c(0,0,0,0))
		    plot(NA,xlim=c(-.8,4.2+log(ncol(x@best_model@sigma))),ylim=c(-0.02,1.05),axes=FALSE,xlab="",ylab="")
		    lines(c(1,1)+0.09,c(-0,1.08))
		    lines(c(0,4.2+log(ncol(x@best_model@sigma))),c(0,0))
		    lines(c(0,2+log(ncol(x@best_model@sigma)))+2.2,c(1.01,1.01))
		    
		    
		    lines(c(0,1)+0.08,c(1.01,1.01))
		    lines(c(0.25,0.25)+0.08,c(1.01,1.015))
		    lines(c(0.5,0.5)+0.08,c(1.01,1.015))
		    lines(c(0.75,0.75)+0.08,c(1.01,1.015))
		    lines(c(0,0)+0.08,c(1.01,1.02))
		    lines(c(1,1)+0.08,c(1.01,1.02))
		    text(0+0.08,1.03,1,cex=0.6)
		    text(0.25+0.08,1.03,0.75,cex=0.5)
		    text(0.50+0.08,1.03,0.50,cex=0.6)
		    text(0.75+0.08,1.03,0.25,cex=0.5)
		    text(1+0.065,1.03,0,cex=0.6)
		    
	    
		    for (k in 1:1){
		      lines(c(0,1)+k*1.1,c(1.01,1.01))
		      lines(c(0.25,0.25)+k*1.1,c(1.01,1.015))
		      lines(c(0.5,0.5)+k*1.1,c(1.01,1.015))
		      lines(c(0.75,0.75)+k*1.1,c(1.01,1.015))
		      lines(c(0,0)+k*1.1,c(1.01,1.02))
		      lines(c(1,1)+k*1.1,c(1.01,1.02))
		      text(0.015+k*1.1,1.03,0,cex=0.6)
		      text(0.25+k*1.1,1.03,0.25,cex=0.5)
		      text(0.5+k*1.1,1.03,0.50,cex=0.6)
		      text(0.75+k*1.1,1.03,0.75,cex=0.5)
		      text(1+k*1.1,1.03,1,cex=0.6)
		    }
		    
		    text(0.5,1.06,expression(rho[kb]),cex=0.9)
		    text(0.5+1.1,1.06,expression(tau[kb]),cex=0.9)
		    lines(c(-.05,-.05),c(0,1))
		    for (h in 1:length(cumpi)){
		      lines(c(-.05,-.09),1-c(cumpi[h],cumpi[h]))
		      if (h!=1){text(-0.08,1-cumpi[h],cumpi[h],cex=0.7,pos=2)}else{text(-0.08,1,0,cex=0.7,pos=2)}
		    }
		    
		    for (h in 1:(length(cumpi)-1)){
		      text(-0.6,1 - 0.5*(cumpi[h]+cumpi[h+1]) ,paste("Class",h))
		      for (loc in 1:max(x@best_model@sigma[h,])){
		        polygon(c(1,1,1-x@best_model@parameters@rho[[h]][loc],1-x@best_model@parameters@rho[[h]][loc])+0.08,(1-cumpi[h])-x@best_model@parameters@proportions[h]*c((loc-1)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc-1)/max(x@best_model@sigma[h,])),col="black",border="white")
		        if (x@best_model@parameters@rho[[h]][loc]!=0){
		          for (mod in 1:length(x@best_model@parameters@tau[[h]][[loc]])){
		            polygon(1.1+c(0,0,x@best_model@parameters@tau[[h]][[loc]][mod],x@best_model@parameters@tau[[h]][[loc]][mod]),(1-cumpi[h])-x@best_model@parameters@proportions[h]*c( (loc-1+mod/length(x@best_model@parameters@tau[[h]][[loc]]))/max(x@best_model@sigma[h,]),(loc-1+(mod-1)/length(x@best_model@parameters@tau[[h]][[loc]]))/max(x@best_model@sigma[h,]),(loc-1+(mod-1)/length(x@best_model@parameters@tau[[h]][[loc]]))/max(x@best_model@sigma[h,]),(loc-1+(mod)/length(x@best_model@parameters@tau[[h]][[loc]]))/max(x@best_model@sigma[h,])),col="black",border="white")
		          }
		        }
		        lines(c(2.2,2.2),c(1.01,1.02))
		        for (j in 1:ncol(x@best_model@sigma)){
		            if (x@best_model@sigma[h,j]==loc){
		              polygon((c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,(1-cumpi[h])-x@best_model@parameters@proportions[h]*c((loc-1)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc-1)/max(x@best_model@sigma[h,])),col="black",border="white")
		            }else{
		              polygon((c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,(1-cumpi[h])-x@best_model@parameters@proportions[h]*c((loc-1)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc-1)/max(x@best_model@sigma[h,])),col="white",border="black")
		            }
		            if ((length(colnames(x@best_model@sigma)[j]!=2))||(any(colnames(x@best_model@sigma)[j]!=paste("X",j,sep="")))){
			            text(mean(c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,1.078,substr(colnames(x@best_model@sigma)[j],1,1),cex=0.5)
			            text(mean(c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,1.06,substr(colnames(x@best_model@sigma)[j],2,2),cex=0.5)
			            text(mean(c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,1.042,substr(colnames(x@best_model@sigma)[j],3,3),cex=0.5)
			        }else{
					    text(mean(c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,1.06,colnames(x@best_model@sigma)[j],cex=0.5)            
		            }
		            lines(c(j,j)*(2+log(ncol(x@best_model@sigma)))/ncol(x@best_model@sigma)+2.2,c(1.01,1.02))
		        }
		      }
		      lines(c(0,3.2),1-c(cumpi[h],cumpi[h]));
		    }
		    
		    
		    for (h in 1:(length(cumpi)-1)){
		      for (loc in 1:max(x@best_model@sigma[h,])){
		        for (j in 1:ncol(x@best_model@sigma)){
		          if (x@best_model@sigma[h,j]==loc){
		            polygon((c(j,j,j-1,j-1)*(2+log(ncol(x@best_model@sigma))))/ncol(x@best_model@sigma)+2.2,(1-cumpi[h])-x@best_model@parameters@proportions[h]*c((loc-1)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc)/max(x@best_model@sigma[h,]),(loc-1)/max(x@best_model@sigma[h,])),col="black",border="white",lwd=1)
		          }
		        }
		      }
		    }
				
				
			

		}
)