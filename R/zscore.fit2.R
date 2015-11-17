library(fdrID)

#--------------------------------------------
#' @title Construct z-scores using a parametric fit distribution to the null Platt scores
#' 
#' @description Construct z-scores using a parametric fit distribution to the null Platt scores
#' 
#' @details A wrapper for the procedure to get null and non-null z-scores from a parametric fit around the bootstrapped
#' null (log) Platt scores. The log of the boostrapped null scores will be taken automatically.
#'
#' @param precomputed.null.scores   Optional vector of precomputed Null Platt scores. If left NULL, bootstrap calculation for null Platt scores will be computed.
#' @param training.dmat             Training set feature vector matrix 
#' @param validation.dmat           Validation set feature vector matrix
#' @param training.labels           Training set labels vector
#' @param validation.labels         Validation set labels vector
#' @param distribution="lg"         Parametric distribution to fit training set log Null Platt scores to. Choices are: "gev" (generalized.extreme.value), "nig" (normal.inverse.gaussian), "sn" (skew.normal) or "lg" (gaussian)
#' @param num.processes             Number of processes to use if computing Null Platt scores to use in the z-score generation process.
#' @param standardizeQ              Standardize the Log Platt scores wrt training set log Null Platt scores?
#' @param num.bs.iter               Number of iterations to use if computing Null Platt scores to use in the z-score generation process.
#' @param C.param                   C-Penalty parameter for the linear kernel SVM steps.
#' @param pvalue.method             Method to compute requisite p-values. Choices are "empirical" or "integral". See details below for warnings.
#' @param printQ                    Print out intermediate values/diagnostic info?
#' @param plotQ                     Display diagnostic plots after each major step in the z-score generation process?
#' 
#' @return all.score.calc.info      A list containing: standardization.flag, fit.distribution.name, fit.info.and.diagnostics,
#' boostrapped.null.platt.scores, score.null.training, score.null.validation, score.nonnull.validation, null.pvalues, 
#' nonnull.pvalues (UN-SMEARED!!!!!!), null.zvalues, nonnull.zvalues.smeared (SMEARED!!!!!!)
#'
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
zscore.fit2 <- function(precomputed.null.scores = NULL, training.dmat, validation.dmat, training.labels, validation.labels,
                       distribution="gaussian", num.processes,
                       standardizeQ = TRUE, 
                       num.bs.iter = 2000,
                       C.param = 0.1,
                       pvalue.method="empirical",
                       printQ = FALSE,
                       plotQ = FALSE) {
  
  
  #Non-repeated protions of z-score construction.
  #Computes:
  # BS Null Platt scores if requested
  # Parametric fit to (possibly standardized) log of BS Null Platt scores
  # Log Non-Null Platt scores on the vlaidation set.
  # Diagnostic plots if requested.
  
  #Perhaps in Ver 3 do this over all dists in parallel????
  
  prep.info <- zscore.fit.prep(precomputed.null.scores, 
                               training.dmat, validation.dmat, 
                               training.labels, validation.labels,
                               distribution,
                               num.processes, 
                               standardizeQ, 
                               num.bs.iter, 
                               C.param, 
                               printQ, plotQ)
  
  #*****NOTE: From here on the Bootstrapped (Pre-Training) Platt score null and Validation set nonnull have been log-ed and may be 
  #(probably are) standardized. They can be found inside prep.info. 
  #
  #The Next step simulates an acceptable Training set Platt score null and Validation set Platt score null utilizing
  #the fit to the Log Bootstrap (Pre-Training) Platt score null.
  
  
  #-----------------------------------------SCORE-NULL on the VALIDATION set------------------------------
  # First obtain a random sample from the fit distributions, the same size as the boostrapped null, to stand 
  # in for the boostrapped training null when computing empirical p-values of the NULL and NONNULL. This new set
  # of simulated log null Platt scores (drawn from the fit distribution) will serve as the training log null
  # Platt scores.
  #------------------------------------------------------------------------------------------------------
    
  #This process may take a few tries if bad (Infinite) z values are encountered. Put the process in a while loop
  #so it repeats a few times if need be. After it repeats a few times (currently 20), if bad p/z values are still
  #being encountered, stop. There is probably another problem. Throw an error.
  #
  #NOTE: score.null.vec.standin will serve as the "training set" log null Platt socres.
  continueQ <- FALSE
  max.val.calc.iters <- 20 #See above note.
  val.calc.iter <- 1
  while(continueQ == FALSE) {
    
    sim.info <- zscore.fit.sim(
      score.null.vec = prep.info$boostrapped.log.null.platt.scores, 
      distribution = prep.info$fit.distribution.name,
      log.null.fit = prep.info$fit.info.and.diagnostics, 
      nrow.validation.dmat = nrow(validation.dmat), 
      nlevels.validation.labels = nlevels(validation.labels),
      printQ = printQ, 
      plotQ = plotQ)
    
    #SPLIT EMP AND INT:?? Decided no for now unless we were computing for all dists at once. We can try that in Ver 3
    
    pvals <- zscore.fit.pvalues(
      score.null.vec.standin = sim.info$score.null.training, 
      score.null.vec.val = sim.info$score.null.validation, 
      score.nonnull.vec.val = prep.info$log.nonnull.platt.scores,
      pvalue.method = pvalue.method,
      distribution = prep.info$fit.distribution.name,
      log.null.fit = prep.info$fit.info.and.diagnostics,
      printQ = printQ,
      plotQ = plotQ)
    
    #CHECK P-VALUES!!!!!
    pvals <- check.for.invalid.pvalues(pvals, printQ=TRUE)
    print(" Done.")
    
    cat(sep="\n\n")
    print("#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#")
    print("                                    Computing Validation Set Z-Values                                  ")
    print("#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#")
    
    cat(sep="\n\n")
    print(paste(" Transforming", pvalue.method, "p-values to z-values..."))
    
    p.null <- pvals$null.p.values
    p.nonnull <- pvals$non.null.p.values
    
    #Possible Null z-values:
    z.null <- qnorm(p.null)    
    
    #Possible Non-null z-values:
    #"Smear" nonnull zero p-values wrt null used if necessary. Return POTENTIAL z values. These still may need some tweeking after
    #this routine finshes.  
    if(length(which(qnorm(p.nonnull)==-Inf)) == 0) {
      print(" No smearing performed for non-null z-values...")
      z.nonnull.smeared <- qnorm(p.nonnull)
      z.smearedQ <- "No"
      
    } else {
      print(" **Smearing negative infinity (potential) non-null z-values...")
      z.nonnull.smeared <- smear.extreme.nonnull.zvalues(p.nonnull, 
                                                         upper.set.zvalue = (-12), 
                                                         mu.factor = (-0.5), 
                                                         p.factor=0.99, 
                                                         printQ=printQ, plotQ=F)[[3]]
      z.smearedQ <- "Yes"
    }
    
    #If problem p and z values are encountered at this NULL construction stage, stop and go get new 
    #randomsamples. Below we perfom further checks on the p-values and z-values. 
    
    #Before moving on, check null and nonnull pvalues for problems.
    if(plotQ == TRUE) {
      pz.diagnostic.info <- check.ps.and.zs(p.null, pnorm(z.nonnull.smeared), printQ=printQ, plotQ=T)
      invisible(readline(prompt="Press [enter] to continue")) 
    } else {
      pz.diagnostic.info <- check.ps.and.zs(p.null, pnorm(z.nonnull.smeared), printQ=printQ, plotQ=F)
    }
    
    problem.idxs <- c(pz.diagnostic.info[[1]], pz.diagnostic.info[[2]], pz.diagnostic.info[[3]], pz.diagnostic.info[[4]])
    if(length(problem.idxs) == 0) {
      print(" ***************************")
      print(paste(" Iteration:", val.calc.iter))
      print(" ***************************")
      print(" No problems were found with either Null or Non-Null z-values obtained with respect to the Training (Stand-in) Log Null Platt Scores. Continuing.")
      continueQ <- TRUE
    } else {
      
      cat(sep="\n\n")
      print(" ***************************")
      print(paste(" Iteration:", val.calc.iter))
      print(" ***************************")
      print(" Problem z-values were found:")
      print(" Null: Positive INF z-value for the indices:")
      print(pz.diagnostic.info[[1]])
      print(" Null: Negative INF z-value for the indices:")
      print(pz.diagnostic.info[[2]])
      print(" Non-Null Positive INF z-value for the indices:")
      print(pz.diagnostic.info[[3]])
      print(" Non-Null: Negative INF z-value for the indices:")
      print(pz.diagnostic.info[[4]])
      
      print(" *************Getting new samples for Training (Stand-in) Null and Validation Null.**********")
      
      val.calc.iter <- val.calc.iter + 1
      
      if(val.calc.iter == max.val.calc.iters) {
        print(" ***************************")
        print(paste(" ******************Iteration:", val.calc.iter))
        print(" ***************************")
        
        stop("Maximum number of iterations reached in trying to get a problem free set of z-values. There is probably something wrong!")
      }
    }     
    #continueQ <- T
  }# end while
  cat(sep="\n\n")
  print(paste("                 ****Acceptable Null Log Platt Score Validation Set Obtained After:", val.calc.iter, "Iterations****                "))
  cat(sep="\n\n")
  
  #Before moving on, eyball null and nonnull pvalues. 
  if(plotQ == TRUE) {
    split.screen(figs=c(1,2))
    screen(1)
    hist(c(z.null, z.nonnull.smeared), main="All z-scores",xlab="z")
    screen(2)
    hist(c(z.nonnull.smeared), main="Smeared Non-Null z-scores", xlab="z")
    close.screen(all=T)
    invisible(readline(prompt="Press [enter] to continue")) 
  }
  
  #Shut off plot device screen:
  if(plotQ == TRUE) {
    gj <- dev.off()
  }
  
  cat(sep="\n\n")
  print("Computation/simulation of necessary z-values done.")
  
  if(pvalue.method == "empirical") {
    cat(sep="\n\n")
    print("NOTE: All p-values computed are empirical wrt the Training Log Null Simulated Platt Scores.")
  }
  
  if(pvalue.method == "integral") {
    cat(sep="\n\n")
    print("NOTE: P-values were computed wrt the fit parametric distribution to the bootstapped Log Null Platt Scores.")
    print("      We have noticed numerical instabilities for very small/large p-values using this method!")
  }
  
  cat(sep="\n\n")
  print("**CAUTION: Double check z-values graphically and with check.ps.and.zs function,")
  print("           especially if you used a GEV fit to the Bootstrap Null!!!!!!!!!!!!!!")  
  print("=============================================================================")
  
  all.score.calc.info <- list(
    standardization.used=standardizeQ,           # Flag if standardization was used.
    fit.distribution.name=distribution,          # Record of parametric distribution name use for fit and Null simulations
    prep.info$fit.info.and.diagnostics,          # Fit info/diagnostics for chosen parametric fit to Boostrapped Log Null Platt Scores
    prep.info$boostrapped.log.null.platt.scores, # Boostrapped Null Platt Scores
    sim.info$score.null.training,                # Simulated stand-in for boostrapped Log Null Platt Scores  
    sim.info$score.null.validation,              # Simulated Validation Log Null Platt Scores  
    prep.info$log.nonnull.platt.scores,          # Validation Log NonNull Platt Scores  
    p.null,                                      # Null p-values from score.null.vec.val, and wrt simulated stand-in for boostrapped Log Null Platt Scores  
    p.nonnull,                                   # Un-smeared non-null p-values from score.nonnull.vec.val, and wrt simulated stand-in for boostrapped Log Null Platt Scores
    z.null,                                      # Null z-values from Null p-values.
    z.smearedQ,                                  # Flag if non-null z-values were smeared.  
    z.nonnull.smeared)                           # Potential (smeared) Non Null z-values. Check acceptability especially if gev fit to bs-null used.
  
  # ******No non-null (officially) un-smeard z-values are returned because they need to be modified. See p.nonnull if you want them!
  
  
  names(all.score.calc.info) <-c(
    "standardization.flag",
    "fit.distribution.name",
    "fit.info.and.diagnostics",
    "boostrapped.null.platt.scores",
    "score.null.training",
    "score.null.validation",
    "score.nonnull.validation",
    "null.pvalues",
    "nonnull.pvalues",         # UN-SMEARED!!!!!!
    "null.zvalues",
    "Were.nonnull.zvalues.smeared?",
    "nonnull.zvalues.smeared"  # Possiblly SMEARED!!!!!!
  )
  
  return(all.score.calc.info)

}