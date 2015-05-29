bootstrap.platt.scores.iteration <- function(dat.mat, lbls, lbls.idxs, svmtyp, kern, pparams) {
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
  iter.prob.mat<-attr(pred, "probabilities")[-mod.bsidx,] #Keep only scores NOT contained in the BS sample
  
  #Get the chunk of known match (non-null) scores from this process.
  #This will be returned to be concatenated with km chunks from other processes:
  iter.score.km<-iter.prob.mat[ cbind(1:nrow(iter.prob.mat),as.numeric(lbls[-mod.bsidx])) ]
  nonnull.vec.chunk <- as.numeric(iter.score.km) #here, there is only one km score per observation
  
#   #Build up vector of known non-match (null) scores.
#   #This will be returned to be concatenated with chunks from other processes:
#   #C++ loop this.
#   #YUK. GET RID OF THIS: 
  iter.score.knm<-NULL
#   #*****On C++ side allocate a iter.score.knm NumericVector of length nrow(iter.prob.mat)
#   #send in iter prob mat. return iter.score.knm or null.vec.chunk
  for(j in 1:nrow(iter.prob.mat)) {
    
    tmp.prob.vec<-iter.prob.mat[j,]
    tmp.prob.vec<-sample(tmp.prob.vec[-as.numeric(lbls[-mod.bsidx][j])],1)  #Just randomly grab one of the KNM scores
    iter.score.knm<-c(iter.score.knm,tmp.prob.vec)
    
  }
  null.vec.chunk <- as.numeric(iter.score.knm)

  return(cbind(null.vec.chunk, nonnull.vec.chunk))

}