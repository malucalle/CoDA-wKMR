###########################################################################
#
# CoDA weighted Kernel Machine Regression 
# file: coda_wKMR_functions.R
#
# Developer: A. Susin and Malu Calle (2019)
########################################################################

coda_wKMR <- function(X,y, covar=NULL, plotting_var=NULL, indcol=1:ncol(X), seqweights, fdr=0.05){
    
  # zero imputation
  
  X_min<-min(as.vector(X[X>0]))
  X = X + 0.5*X_min-min(X, 0.5*X_min)
  
  # default values

  if (missing(seqweights)){ # if not changed by user
    seqweights<-c(seq(0.1,0.7,length.out = 3),seq(1,5,length.out = 5))
  }
  lseqw<-length(seqweights)
  
  if (missing(indcol)){  # if not changed by user
    indcol<-1:ncol(X)  # analisy of all variables
  }
  lindcol<-length(indcol)
  
  slope<-matrix(0,nrow=lindcol,ncol=3)
  outseq<-matrix(0,nrow=lindcol,ncol=lseqw)
  for(i in (1:lindcol)){
    index<-indcol[i]
    print(index)
    
    for (j in (1:lseqw)){
      weights<-rep(1,ncol(X))
      weights[index]<-seqweights[j]
      pval<-wKMR(X, y, cov=covar, w=weights)
      outseq[i,j]<-pval
    }
    p_ref<-wKMR(X, y, cov=covar)
    
    sl<-(summary(lm(-log10(outseq[i,])~seqweights))$coefficients)[2,1]
    psl<-(summary(lm(-log10(outseq[i,])~seqweights))$coefficients)[2,4]
    slope[i,1]<-index
    slope[i,2]<-as.numeric(sl)
    slope[i,3]<-as.numeric(psl)
  }
  Xslope<-as.data.frame(slope)
  # write.csv(X.slope, "xslope.csv")
  Xslope[,4]<-p.adjust(Xslope[,3], method="fdr")
  colnames(Xslope)<-c("taxa","slope","p-value","adj.p-value")
  rownames(Xslope)<-colnames(X[,indcol])
  Xsig<-Xslope[(Xslope[,2]>0)&(Xslope[,4]<fdr),]
  osig<-order(Xsig[,2],decreasing = T )
  Xsig[osig,]
  oslope<-order(Xslope[,2],decreasing = T )

  result <- list("Importance_table"=Xslope,"Significative_table"=Xsig[osig,], 
                 "seqweights"=seqweights, "p_values"=outseq, "p_ref"=p_ref)
  return(result)
  
}


w_geoM <- function(w, x){
  y <- x/w       # Redefine y
  sw <- sum(w)   # Sum of weights
  w_g <- exp(sum(w*log(y))/sw)   # Weighted geometric mean
  return(w_g)    # Return the value
}



log_w_geoM <- function(w, x){
  y <- x/w   # Redefine y
  sw <- sum(w)   # Sum of weights
  log_w_g <- sum(w*log(y))/sw   # Weighted geometric mean

  return(log_w_g)   # Return the value
}


w_dist <- function(w, x1, x2){
  w_d <- sqrt(sum(w*(log(x1)-log_w_geoM(w,x1) - log(x2)+log_w_geoM(w,x2))^2))
  return(w_d)
}



wKMR <- function(x, y, w = rep(1, ncol(x)), cov = NULL){

  #------------------------------------------------------------------------------#
  # Build distance matrix
  #------------------------------------------------------------------------------#

  # Definte dist as dist{proxy}
    dist <- proxy::dist

  # wDM: weigthed Distance Matrix
  # Build the matrix
  wDM <- as.matrix(dist(x,method = w_dist, w=w ))
  # 
  # m = nrow(x);
  # n = m;
  # wDMT=matrix(0.0, nrow=m, ncol=n);
  # 
  # for (i in (1:m)){
  #   row_i = x[i,];
  #     for (j in (i:n)){
  #       row_j = x[j,];
  #       wDMT[i,j]=w_dist(w,row_i,row_j);
  #       wDMT[j,i]=wDMT[i,j];
  #   }
  # }
  # 
  # err=max(wDM-wDMT);
  # cat(sprintf('err = %d \n', err));
  #----------------------------------------#
  # P-value for weighted data
  #----------------------------------------#

  # Distance to Kernel
  K1<-D2K(wDM)
  # Extract the weighted p-value
  wp <- MiRKAT(y = as.numeric(y),
               X = cov,
               Ks = K1)
  return(wp)

}


#------------------------------------------------------------------------------#


plot_slope <- function( idvar, seqweights, outseq, Xslope, p_ref) {

# Load library ggplot2
#  suppressMessages(library(ggplot2))
  plotData<-data.frame(seqweights,-log10(outseq[idvar,]))
  titol<-rownames(Xslope[idvar,])
  slope_plot <- ggplot(plotData, aes(x=plotData[,1],y=plotData[,2]))+geom_point()+ 
    geom_hline(yintercept=-log10(p_ref), size=0.5, color="black")+
    ggtitle(titol) + xlab("weights") + ylab("-log10(pvalue)") #subtitle="From midwest dataset"


  return(slope_plot)
}


plotVarImportance_coda_wKMR <- function(Xslope) {
  
# Load library tidyverse
#  suppressMessages(library(tidyverse))
  
  plotData<-data.frame(x=rownames(Xslope),y=Xslope[,2])
  titol<-"variable importance"
  plotData= plotData %>% arrange(y) %>% mutate(x=factor(x,x))
  varImportance <- ggplot(plotData, aes(x,y))+geom_point()+ 
    geom_hline(yintercept=0, size=0.5, color="black")+
    geom_segment(aes(x=x,xend=x,y=min(y), yend=y), color='blue')+
    ggtitle(titol)+ ylab("slope")+ xlab(" ") + coord_flip() 
  
  return(varImportance)
}


plotVarSlope_coda_wKRM <- function(plotting_var, seqweights, outseq, Xslope, p_ref) {
  
if (!is.null(plotting_var)){
  lplotting_var<-length(plotting_var)
  for (i in 1:lplotting_var){
    idvar= i  #variable_id ordered in Xslope
    slopePlot <- plot_slope( idvar, seqweights, outseq, Xslope, p_ref)
    fileout=paste("slopePlot_var",toString(idvar),sep = "")
    filename=paste(fileout,".pdf",sep = "")
    ggsave(filename)
  }
}

}





