bootstrap.platt.scores.parallel<-function(num.processes,dat.mat,lbls,nbs,svmtyp="C-classification",kern="linear",pparams=0.1) {
  
  #parse lbls into groups:
  lbls.idxs<-lapply(1:nlevels(lbls),function(x){which(lbls==levels(lbls)[x])})
  
  #Initialize score arrays this time instead of iteratively building them up. Make make things run faster
  #   null.vec<-rep(-1.0,nbs*nrow(dat.mat))    #To keep track of KM (non-null scores)
  #   nonnull.vec<-rep(-1.0,nbs*nrow(dat.mat)) #Save only 1 random KNM score per obs-iteration
#   null.vec<-NULL
#   nonnull.vec<-NULL
#   for(i in 1:nbs) {
#     
#   }
  
  
  registerDoMC(num.processes)
  print(paste("Using",getDoParWorkers(), "processes."))
  
  a <- seq(1:10)
  t1<-Sys.time()
  junk <- foreach(i=1:nbs, .combine="rbind") %dopar% {  
    bootstrap.platt.scores.iteration2(dat.mat, lbls, lbls.idxs, svmtyp, kern, pparams)
  }
  t2<-Sys.time()
  print(paste("Time: ",difftime(t2,t1)))

 return(junk)
  
}