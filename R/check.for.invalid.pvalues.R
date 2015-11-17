#--------------------------------------------
#' @title check.for.invalid.pvalues
#' 
#' @description Check to see if p-values < 0 or > 1 are encountered in a list of null and non-null p-values.
#' 
#' @details We noticed that some of the distributions fit to the log null Platt scores can yield invalid p-values 
#' (mostly slightly negative) deep into the tails. This happens especially for the NIG distributions. This is a handy 
#' utility function to check a for bad p-values. It's general to check for bad p-values amongst null and non-null
#' sets, assuming they have been returned from the \code{zscore.fit.pvalues} function. If bad p-values are encountered, 
#' they are replaced with 0 (for negative p-values) or 1 (for positive p-values).
#' 
#' @param pvalue.info    A two element list. Each element is a vector of p-values. Cf. \code{zscore.fit.pvalues}
#' @param printQ         Diagnostic info?  
#' 
#' @return A potentially modified pvalue.info type list with bad p-valued replaced.
#'
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
check.for.invalid.pvalues <-function(pvalue.info, printQ=FALSE) {
  
  pnulls <- pvalue.info$null.p.values
  pnonnulls <- pvalue.info$non.null.p.values
    
  neg.null.p.idxs <- which(pnulls < 0)
  neg.nonnull.p.idxs <- which(pnonnulls < 0)
  
  gt1.null.p.idxs <- which(pnulls > 1)
  gt1.nonnull.p.idxs <- which(pnonnulls > 1)

  if(length(neg.null.p.idxs) > 0 ) {
    
    if(printQ == TRUE) {
      print("   *!*!*!*!*!*!*!*!* P-VALUE PROBLEM!! These null p-values are less than zero and will be replaced as 0:")
      print(neg.null.p.idxs)
      print("The problem p-values are:")
      print(pnulls[neg.null.p.idxs])
    }
    
    pnulls[neg.null.p.idxs] <- 0
    
  }
  
  if(length(neg.nonnull.p.idxs) > 0 ) {
    
    if(printQ == TRUE) {
      print("   *!*!*!*!*!*!*!*!* P-VALUE PROBLEM!! These non-null p-values are less than zero and will be replaced as 0:")
      print(neg.nonnull.p.idxs)
      print("The problem p-values are:")
      print(pnonnulls[neg.nonnull.p.idxs])
    }
    
    pnonnulls[neg.nonnull.p.idxs] <- 0
    
  }
  
  if(length(gt1.null.p.idxs) > 0 ) {
    
    if(printQ == TRUE) {
      print("   *!*!*!*!*!*!*!*!* P-VALUE PROBLEM!! These null p-values are greater than 1 and will be replaced as 1:")
      print(gt1.null.p.idxs)
      print("The problem p-values are:")
      print(pnulls[gt1.null.p.idxs])
    }
    
    pnulls[gt1.null.p.idxs] <- 1
    
  }
  
  if(length(gt1.nonnull.p.idxs) > 0 ) {
    
    if(printQ == TRUE) {
      print("   *!*!*!*!*!*!*!*!* P-VALUE PROBLEM!! These non-null p-values are greater than 1 and will be replaced as 1:")
      print(gt1.nonnull.p.idxs)
      print("The problem p-values are:")
      print(pnonnulls[gt1.nonnull.p.idxs])
    }
    
    pnonnulls[gt1.nonnull.p.idxs] <- 1
    
  }
  
  p.values <- list(pvalue.info$pvalue.method, pnulls, pnonnulls)
  names(p.values) <- c("calculation.method","null.p.values", "non.null.p.values")
  
  return(p.values)
  
}

#P-values for unknowns can be treated by default as non-null p-values (they are the p-values corresponding to the svm rendered ID)
check.for.invalid.pvalues.unk <-function(unk.pvalues, printQ=FALSE) {
  
  pnonnulls <- unk.pvalues
  
  neg.nonnull.p.idxs <- which(pnonnulls < 0)
  
  gt1.nonnull.p.idxs <- which(pnonnulls > 1)
    
  if(length(neg.nonnull.p.idxs) > 0 ) {
    
    if(printQ == TRUE) {
      print("   *!*!*!*!*!*!*!*!* P-VALUE PROBLEM!! These p-values for unknowns are less than zero and will be replaced as 0:")
      print(neg.nonnull.p.idxs)
      print("The problem p-values are:")
      print(pnonnulls[neg.nonnull.p.idxs])
    }
    
    pnonnulls[neg.nonnull.p.idxs] <- 0
    
  }
    
  if(length(gt1.nonnull.p.idxs) > 0 ) {
    
    if(printQ == TRUE) {
      print("   *!*!*!*!*!*!*!*!* P-VALUE PROBLEM!! These p-values for unknowns are greater than 1 and will be replaced as 1:")
      print(gt1.nonnull.p.idxs)
      print("The problem p-values are:")
      print(pnonnulls[gt1.nonnull.p.idxs])
    }
    
    pnonnulls[gt1.nonnull.p.idxs] <- 1
    
  }
    
  return(pnonnulls)
  
}