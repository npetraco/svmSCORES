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