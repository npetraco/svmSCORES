#--------------------------------------------
#' @title bootstrap.platt.scores.parallel
#' @description Estimate svm Platt-score null and non-null distributions by group-wise bootstrap. Parallel code.
#' 
#' @details Estimate svm Platt-score null (known non-matching) and non-null (known matching) distributions via a 
#' group-wise bootstrap. Adapted from Storey and Tibshirani permutation method in PNAS. The SVM code used is from 
#' e1071 package. In an effort to decreases depencence between scores and also makes the null distribution of scores 
#' a bit more conservative, the code randomly selects which KNM score we keep for a given observation. Also
#' it only keeps scores from observations NOT in the BS sample. CAUTION: Because the algorithm bootstraps each group 
#' independently, each group should have a fairly large number (>10) of samples.
#' 
#' The "parallelization" strategy of this routine is simple, spread the bootstrap iterations across a specified 
#' number of processes. Usually the number of processes is the number of cores or logical cores (some cores can 
#' handle multiple processes) avalible on the computer you are using.
#' 
#' It's recommended to do a few "dry runs" varying the number of processes while using only a few bootstrap 
#' iterations (10-100) to see what kind of speed up profile cqn be expected. The time of the full run can be
#' approximated as: time-for-dry-run/nbs-for-dry-run * nbs-for-full-run.
#'
#' @param dat.mat     A n by p data frame or matrix of variables. One row per pattern to classify.
#' @param lbls        Data ID labels
#' @param nbs         Number of bootstrap iterations.
#' @param svmtyp      Support vector machine type. See e1071 package documentation.
#' @param kern        Kernel type. See e1071 package documentation.
#' @param pparams     Kernel parameters. See e1071 package documentation.
#' @param timerQ      Spit out timing information?
#' 
#' @return column matrix of null and non-null scores
#' 
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
bootstrap.platt.scores.parallel<-function(num.processes=1,dat.mat,lbls,nbs=10,svmtyp="C-classification",kern="linear",pparams=0.1, timerQ=FALSE) {
  
  #parse lbls into groups:
  lbls.idxs<-lapply(1:nlevels(lbls),function(x){which(lbls==levels(lbls)[x])})
  
  if(timerQ==TRUE){
    t1<-Sys.time()
  }
  
  #The parallelization. Just spread the bootstrap iterations across the processes:
  registerDoMC(num.processes)
  print(paste("Using",getDoParWorkers(), "processes."))

  null.nonnull.scores <- foreach(i=1:nbs, .combine="rbind") %dopar% {  
    bootstrap.platt.scores.iteration2(dat.mat, lbls, lbls.idxs, svmtyp, kern, pparams)
  }
  
  if(timerQ==TRUE){
    t2<-Sys.time()
    print(paste("Time: ",difftime(t2,t1)))
  }
  
  colnames(null.nonnull.scores) <- c("null","nonnull")
  
  return(null.nonnull.scores)
  
}