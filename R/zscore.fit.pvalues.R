#--------------------------------------------------------------------------
#Internal function. Compute requisite p-values in the z-score fit process.
#--------------------------------------------------------------------------
zscore.fit.pvalues <- function(
                       score.null.vec.standin, 
                       score.null.vec.val, 
                       score.nonnull.vec.val,
                       pvalue.method="empirical",
                       distribution = NULL,
                       log.null.fit = NULL,
                       printQ = FALSE,
                       plotQ = FALSE) {
  
  #Use the above null score stand in right away to get p and z from the validation sets and see if there are 
  #any problems.
  
  cat(sep="\n\n")
  print(paste(" Computing", pvalue.method, "p-values on VALIDATION SET log null and non-null Platt Scores..."))
  
  if(pvalue.method=="integral" & is.null(log.null.fit)) {
    stop("Indicated integral p-values but you forgot to send in the parameters for the distribution. Specify log.null.fit!")
  }
  
  #VALIDATION SET Null p-values and z-values wrt substituting the fit sampled (log) null as emiprical null:
  if(pvalue.method=="empirical") {
    
    #VALIDATION SET NULL empirical p-values wrt substituting the fit sampled (log) null as emiprical null:
    p.null <- empiricalPvalues(score.null.vec.standin, score.null.vec.val)
    #VALIDATION SET NONNULL empirical p-values wrt substituting the fit sampled (log) null as emiprical null:
    p.nonnull <- empiricalPvalues(score.null.vec.standin, score.nonnull.vec.val)
    
  } else if(pvalue.method=="integral") {
    
    if(distribution=="gev") {
      
      #VALIDATION SET NULL/NONNULL integral p-values wrt the gev fit bootstrapped (training) log null:
      p.null <- pgev(score.null.vec.val, loc=log.null.fit[[1]][1], scale=log.null.fit[[1]][2], shape = log.null.fit[[1]][3], lower.tail = F)
      p.nonnull <- pgev(score.nonnull.vec.val, loc=log.null.fit[[1]][1], scale=log.null.fit[[1]][2], shape = log.null.fit[[1]][3], lower.tail = F)
      
    } else if (distribution=="nig") {
      
      #VALIDATION SET NULL/NONNULL integral p-values wrt the nig fit bootstrapped (training) log null:
      p.null <- 1 - pnig(score.null.vec.val, alpha=log.null.fit[[1]][1], beta=log.null.fit[[1]][2], delta=log.null.fit[[1]][3], mu=log.null.fit[[1]][4])
      p.nonnull <- 1 - pnig(score.nonnull.vec.val, alpha=log.null.fit[[1]][1], beta=log.null.fit[[1]][2], delta=log.null.fit[[1]][3], mu=log.null.fit[[1]][4])
      
    } else if (distribution=="sn") {
      
      #VALIDATION SET NULL/NONNULL integral p-values wrt the sn fit bootstrapped (training) log null:
      p.null <- 1 - psn(score.null.vec.val, xi=log.null.fit[[1]][1], omega=log.null.fit[[1]][2], alpha=log.null.fit[[1]][3])
      p.nonnull <- 1 - psn(score.nonnull.vec.val, xi=log.null.fit[[1]][1], omega=log.null.fit[[1]][2], alpha=log.null.fit[[1]][3])
      
    } else if (distribution=="lg") {
      
      #VALIDATION SET NULL/NONNULL integral p-values wrt the lg fit bootstrapped (training) log null:
      p.null <- pnorm(score.null.vec.val, mean=log.null.fit[[1]][1], sd=log.null.fit[[1]][2], lower.tail = F)
      p.nonnull <- pnorm(score.nonnull.vec.val, mean=log.null.fit[[1]][1], sd=log.null.fit[[1]][2], lower.tail = F)
      
    } else {
      stop("Specified integral p-values, but forgot to specify a distribution!")
    }
    
  } else {
    stop("Pick a method to compute p-values. Your choices are empirical or integral!")
  }
  
  p.values <- list(pvalue.method, p.null, p.nonnull)
  names(p.values) <- c("calculation.method","null.p.values", "non.null.p.values")
  
  return(p.values)
  
}