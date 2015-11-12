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
posterior.probs.for.unks2 <- function(training.dmat, training.labels, C.param, test.dmat, 
                                      standardizeQ=FALSE, mean.log.null.score=NULL, sd.log.null.score=NULL,
                                      pvalue.method="empirical",
                                      null.vec.training=NULL,
                                      distribution.name=NULL, distribution.fit.info=NULL,
                                      pp.point.est.func, pp.uci.est.func, pp.lci.est.func) {
  
  if(standardizeQ == TRUE & is.null(mean.log.null.score)) {
    stop("Standardization of Test Platt scores requested, but no reference mean null log score was provided! Specify std.mean.log.null.score.")
  }
  if(standardizeQ == TRUE & is.null(sd.log.null.score)) {
    stop("Standardization of Test Platt scores requested, but no reference sd null log score was provided! Specify std.sd.log.null.score.")
  }
  
  #Train the SVM:
  training.svm.model <- svm(training.dmat, training.labels, scale=FALSE, type="C-classification", kernel="linear", cost=C.param, fitted=TRUE, probability=TRUE)
  pred.info <- predict(training.svm.model, test.dmat, probability=TRUE)
  
  #Extract predicted labels and Platt Scores for test set:
  test.preds <- predict(training.svm.model, test.dmat, probability=TRUE)
  
  #Predicted labels:
  test.pred.lbls <- test.preds[1:length(test.preds)]
  
  #Platt scores for predicted labels:
  test.log.platt.scores <- log(attr(test.preds, "probabilities")[cbind(1:length(test.pred.lbls), as.numeric(test.pred.lbls))])
  
  #Standardize scores wrt input mean and sd if requested:
  if(standardizeQ == TRUE) {
    test.log.platt.scores <- (test.log.platt.scores - mean.log.null.score)/sd.log.null.score
  }
    
  #Get z-scores on the test set:
  z.unk <- zscore.unk(test.score.vec = test.log.platt.scores, 
                      pvalue.method = pvalue.method, 
                      null.vec.training = null.vec.training, 
                      distribution = distribution.name, 
                      dist.fit.info = distribution.fit.info)
  
  if(length(which(z.unk==-Inf))>0) {
    print("*********Smearing small p-values of test set!")
    print("This situation is not uncommon, especially for empirical, nig, sn and gev p-values.")
    print("HOWEVER, we are doing this for p-values on unknowns, so take the resulting interpolated posterior probabilities with a grain of salt.")
    print("For these estimated probabilities, all we really can say is that they are smaller than the smallest posterior probability we could actually get a p-value for.")
    z.unk <- smear.extreme.nonnull.zvalues(pnorm(z.unk),
                                           upper.set.zvalue = (-12), mu.factor = (-0.5), 
                                           p.factor=0.99, plotQ=F)[[3]]
  }
  
  pps <- pp.point.est.func(z.unk)
  ucis <- pp.uci.est.func(z.unk)
  lcis <- pp.lci.est.func(z.unk)
  
  posterior.prob.mat <- data.frame(test.pred.lbls, z.unk, ucis, pps, lcis)
  colnames(posterior.prob.mat) <- c("Pred.Label", "z-score", "Upper.CI", "Post.Prob.Est.", "Lower.CI")
  
  return(posterior.prob.mat)
  
}