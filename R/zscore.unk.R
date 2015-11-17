#--------------------------------------------
#' @title zscore.unk
#' 
#' @description Compute z-scores on a set of test feature vectors (can be only one). 
#' 
#' @details At this pont Platt scores on a test set should have been obtained from an SVM fit to a training set. That should have
#' occured in the first phase of \code{posterior.probs.for.unks}. The Platt scores on the test set are transformed into z-scores
#' (one for each unknown) via the distribution fit to the bootstrapped null. If the distribution name and fit info is supplied
#' integral p-values are computed. If the validation set of null Platt scores is supplied, empirical p-values will be computed.
#' The p-values on the test set are then checked for invalid values (< 0 or > 1) and then transformed to z-values. 
#'
#' @param test.score.vec    A set of Platt scores on a test set. One score for each "unknown"
#' @param pvalue.method     "empirical" or "integral"
#' @param null.vec.training Validation set of null scores needed if pvalue.method = "empirical"
#' @param distribution      Name of distribution info pvalue.method = "integral"
#' @param dist.fit.info     Distribution fit info if pvalue.method = "integral"
#' 
#' @return a vector of z-values on the test set.
#'
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
zscore.unk <- function(test.score.vec, pvalue.method, null.vec.training=NULL, distribution=NULL, dist.fit.info=NULL) {
    
  if(pvalue.method=="integral" & is.null(dist.fit.info)) {
    stop("Indicated integral p-values but you forgot to send in the parameters for the distribution. Specify dist.fit.info!")
  }

  if(pvalue.method=="integral" & is.null(distribution)) {
    stop("Indicated integral p-values but you forgot to specify the distribution name! Your choices are: gev, nig, sn or lg")
  }
  
  if(pvalue.method=="empirical" & is.null(null.vec.training)) {
    stop("Indicated empirical p-values but you forgot to send in the a null score sample. Specify null.vec.training!")
  }
  
  #VALIDATION SET Null p-values and z-values wrt substituting the fit sampled (log) null as emiprical null:
  if(pvalue.method=="empirical") {
    
    pvalues <- empiricalPvalues(null.vec.training, test.score.vec)
    
  } else if(pvalue.method=="integral") {
    
    if(distribution=="gev") {
      
      pvalues <- pgev(test.score.vec, loc=dist.fit.info[[1]][1], scale=dist.fit.info[[1]][2], shape = dist.fit.info[[1]][3], lower.tail = F)
      
    } else if (distribution=="nig") {
      
      pvalues <- 1 - pnig(test.score.vec, alpha=dist.fit.info[[1]][1], beta=dist.fit.info[[1]][2], delta=dist.fit.info[[1]][3], mu=dist.fit.info[[1]][4])
      
    } else if (distribution=="sn") {
      
      pvalues <- 1 - psn(test.score.vec, xi=dist.fit.info[[1]][1], omega=dist.fit.info[[1]][2], alpha=dist.fit.info[[1]][3])
      
    } else if (distribution=="lg") {
      
      pvalues <- pnorm(test.score.vec, mean=dist.fit.info[[1]][1], sd=dist.fit.info[[1]][2], lower.tail = F)
      
    } else {
      stop("Specified integral p-values, but forgot to specify a distribution!")
    }
    
  } else {
    stop("Pick a method to compute p-values. Your choices are empirical or integral!")
  }
  
  #Check for p-values < 0 or > 1. This happens especially for nig fits: 
  pvalues <- check.for.invalid.pvalues.unk(pvalues, printQ=TRUE)
  
  return(qnorm(pvalues))
  
}