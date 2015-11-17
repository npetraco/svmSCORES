#--------------------------------------------
#' @title posterior.probs.for.unks
#' 
#' @description Determines labels and computes posterior probabilities for unknown feature vectors.
#' 
#' @details XXXX
#'
#' @param training.dmat         Training data matrix
#' @param training.labels       Labels for training data
#' @param C.param               Penalty parameter for linear SVM
#' @param test.dmat             Test data matrix
#' @param standardizeQ          Standardize all feature vectors?
#' @param mean.log.null.score   The mean to standardize by.
#' @param sd.log.null.score     The sd to standardize by.
#' @param pvalue.method         "empirical" or "integral"
#' @param null.vec.training     Training log null values from boostrap calculation. Cf. \code{bootstrap.platt.scores.parallel} or Cf \code{bootstrap.platt.scores}
#' @param distribution.name     Name of a distribution log null boostrapped Platt scores was fit too. Current choices are 
#' "gev" (generalized extreme value), "nig" (normal inverse gaussian), "sn" (skew normal) and "lg" (gaussian) 
#' @param distribution.fit.info Fit information
#' @param pp.point.est.func     The posterior probability point estimate interpolation function.
#' @param pp.uci.est.func       The posterior probability upper credibility/confidence bound interpolation function from training/validation set fit
#' @param pp.lci.est.func       The posterior probability lower credibility/confidence bound interpolation function from training/validation set fit
#' 
#' @return An info matrix containing for each input test feature vector: the predicted label, the z-score,  posterior prob. upper CI, posterior prob. mean,  posterior prob. lower CI
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