#--------------------------------------------
#' @title Generate non-null (matching) Platt scores from a validation set.
#' @description Estimate non-null Platt-scores using a trained SVM and a validation set not used in training.
#' 
#' @details User can opt to get non-null Platt scores off a single run on a held out validation set. These, combined with a set of (perferabley IID) 
#' null scores can go on to be used to fit an fdr model. Null Platt scores (non-match) on the input validation set are NOT returned with this function 
#' as they are known to be dependent. To obtain a set of null scores for an fdr fit, its suggested to either (1) fit a bootstraped set of null Platt
#' scrores to a known parametric form (e.g. a Gaussian or NIG) and sample a set if IID null scores from that fit distribution, (2) use the bootstrap 
#' procedure on the held out validation set to obtain a null AND a non-null set of Platt scores (i.e. don't use the non-null scores obtained from 
#' this function).
#'
#' @param dat.tr    Training data set.
#' @param lbls.tr   Training set labels.
#' @param dat.val   Validation set data set.
#' @param lbls.val  Validation set labels.
#' @param svmtyp    SVM type. See e1071 documentation.
#' @param kern      Kernel type. See e1071 documentation.
#' @param pparams   Penalty parameters. See e1071 documentation.
#' 
#' @return An array of non-null scires on the validation set.
#' 
#' @examples
#' XXXX
#--------------------------------------------
nonnull.logplatt.on.validation.set <- function(dat.tr, lbls.tr, dat.val, lbls.val, svmtyp, kern, pparams) {
  
  #Fit SVM to the training sample
  svm.model<-e1071::svm(dat.tr, lbls.tr, scale=FALSE, type=svmtyp, kernel=kern, cost=pparams, fitted=TRUE, probability=TRUE)
  
  #Get Platt scores on the validation set, not used in training.
  pred<-predict(svm.model,dat.val,probability=TRUE)
  val.prob.mat<-attr(pred, "probabilities")
  
  #Get the known match (non-null) scores only. Discard the known non-match (null) scores with the philisophy that they are dependent.
  val.score.nonnull<-as.numeric(val.prob.mat[ cbind(1:nrow(val.prob.mat),as.numeric(lbls.val)) ])
  
  return(val.score.nonnull)
}