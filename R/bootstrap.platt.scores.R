#--------------------------------------------
#' @title bootstrap.platt.scores
#' @description Estimate svm Platt-score null and non-null distributions by group-wise bootstrap. Serial code.
#' 
#' @details Estimate svm Platt-score null (known non-matching) and non-null (known matching) distributions via a 
#' group-wise bootstrap. Adapted from Storey and Tibshirani permutation method in PNAS. The SVM code used is from 
#' e1071 package. In an effort to decreases depencence between scores and also makes the null distribution of scores 
#' a bit more conservative, the code randomly selects which KNM score we keep for a given observation. Also
#' it only keeps scores from observations NOT in the BS sample. CAUTION: Because the algorithm bootstraps each group 
#' independently, each group should have a fairly large number (>10) of samples.
#'
#' @param dat.mat     A n by p data frame or matrix of variables. One row per pattern to classify.
#' @param lbls        Data ID labels
#' @param nbs         Number of bootstrap iterations.
#' @param svmtyp      Support vector machine type. See e1071 package documentation.
#' @param kern        Kernel type. See e1071 package documentation.
#' @param pparams     Kernel parameters. See e1071 package documentation.
#' 
#' @return list of null and non-null scores
#' 
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
bootstrap.platt.scores<-function(dat.mat,lbls,nbs,svmtyp="C-classification",kern="linear",pparams=0.1) {
  
  #parse lbls into groups:
  lbls.idxs<-lapply(1:nlevels(lbls),function(x){which(lbls==levels(lbls)[x])})
  
  #Initialize score arrays this time instead of iteratively building them up. Make make things run faster
  #   null.vec<-rep(-1.0,nbs*nrow(dat.mat))    #To keep track of KM (non-null scores)
  #   nonnull.vec<-rep(-1.0,nbs*nrow(dat.mat)) #Save only 1 random KNM score per obs-iteration
  null.vec<-NULL
  nonnull.vec<-NULL
  for(i in 1:nbs) {
    
    #Grab a bootstrap sample, bootstrapping out of each GROUP invididually:
    bsidx<-unlist(as.vector(sapply(1:length(lbls.idxs),function(x){sample(lbls.idxs[[x]],length(lbls.idxs[[x]]), replace=TRUE )})))
    bsdat<-dat.mat[bsidx,]
    
    rownames(bsdat)<-NULL #Shuts off an annoying warning
    bslbl<-lbls[bsidx]      
    
    #Fit SVM to the bootstrap sample
    svm.model<-svm(bsdat, bslbl, scale=FALSE, type=svmtyp, kernel=kern, cost=pparams, fitted=TRUE, probability=TRUE)
    
    #Get Platt scores from the whole data set
    pred<-predict(svm.model,dat.mat,probability=TRUE)
    mod.bsidx<-unique(sort(bsidx))
    iter.prob.mat<-attr(pred, "probabilities")[-mod.bsidx,] #Keep only scores not contained in the BS sample
    #rownames(iter.prob.mat)<-NULL
    #print(paste("Iter:",i))
    #print(iter.prob.mat)
    
    #Build up vector of known match (non-null) scores:
    iter.score.km<-iter.prob.mat[ cbind(1:nrow(iter.prob.mat),as.numeric(lbls[-mod.bsidx])) ]
    #print(length(iter.score.km))
    nonnull.vec<-c(nonnull.vec,as.numeric(iter.score.km)) #here, there is only one km score per observation
    
    #Build up vector of known non-match (null) scores:
    iter.score.knm<-NULL
    for(j in 1:nrow(iter.prob.mat)) {
      
      tmp.prob.vec<-iter.prob.mat[j,]
      tmp.prob.vec<-sample(tmp.prob.vec[-as.numeric(lbls[-mod.bsidx][j])],1)  #Just randomly grab one of the KNM scores
      iter.score.knm<-c(iter.score.knm,tmp.prob.vec)
      
    }
    iter.score.knm<-as.numeric(iter.score.knm)
    null.vec<-c(null.vec,iter.score.knm)
    
    print(paste("B.S. Iter: ",i))
  }
    
  print(paste("Maximum null score    : ",max(null.vec)))
  print(paste("Minimum non-null score: ",min(nonnull.vec)))
  
  return(list(null.vec,nonnull.vec))
  
}