#--------------------------------------------
#' @title XXXX
#' 
#' @description XXXXX
#' 
#' @details XXXX
#'
#' @param XXXX
#' 
#' @return XXXX
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
  
  return(qnorm(pvalues))
  
}