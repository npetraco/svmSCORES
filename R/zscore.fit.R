#--------------------------------------------
#' @title Construct z-scores using a parametric fit distribution to the null Platt scores
#' 
#' @description Construct z-scores using a parametric fit distribution to the null Platt scores
#' 
#' @details A wrapper for the procedure to get null and non-null z-scores from a parametric fit around the bootstrapped
#' null (log) Platt scores. The log of the boostrapped null scores will be taken automatically.
#'
#' @param XXXX XXXXXXX
#' 
#' @return XXXX
#'
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
zscore.fit <- function(PSNTMP, training.dmat, validation.dmat, training.labels, validation.labels, distribution="gaussian", num.processes, 
                       standardizeQ=TRUE, 
                       num.bs.iter=2000,
                       C.param = 0.1,
                       plotQ=FALSE) {
  
  if(!(distribution %in% c("gev", "nig", "sn", "lg")) ) {
    stop("Your choices in distribution are: generalized.extreme.value (gev), normal.inverse.gaussian (nig), skew.normal (sn) or gaussian (lg).")
  }
  
  if(distribution=="gev") {
    dist.nme <- "generalized extreme value"
  } else if (distribution=="nig") {
    dist.nme <- "normal.inverse.gaussian"
  } else if (distribution=="sn") {
    dist.nme <- "skew.normal"
  } else if (distribution=="lg") {
    dist.nme <- "gaussian"
  }
  
  
  print("Beginning Bootstrap Computation of the Null Platt Scores")
  #platt.scores <- bootstrap.platt.scores.parallel(
  #  num.processes=num.processes, 
  #  dat.mat=training.dmat, 
  #  lbls=training.labels, 
  #  nbs=num.bs.iter, 
  #  svmtyp="C-classification", 
  #  kern="linear", 
  #  pparams=C.param,
  #  timerQ=T)
  platt.scores <- PSNTMP
  print("========================================================")
  
  score.null.vec<-platt.scores[,1]
  score.nonnull.vec<-platt.scores[,2]
  
  print(paste("Number of Bootstrapped Null Platt Scores:", length(score.null.vec)))
  print(paste("Extent of Null right tail:               ", max(score.null.vec)))
  print(paste("Extent of Non-Null left tail:            ", min(score.nonnull.vec)))

  # We don't need score.nonnull anymore:
  remove(score.nonnull.vec)

  #-----------------------------------------------------------------------------------------------------
  # Parametric fits for the standardized log null (KNM) bootstrapped Platt scores
  #-----------------------------------------------------------------------------------------------------
  cat(sep="\n\n")
  print(paste("Fitting log null Platt score bootstrap sample to a", dist.nme, "distribution:"))
  if(distribution=="gev") {
    
    # Examine null fit for the extreme value distribution:
    log.null.fit <- null.logplatt.gevd.fit(score.null.vec, standardizeQ=standardizeQ, plotQ=F)
    
  } else if(distribution=="nig") {
    
    # Examine null fit for the extreme value distribution:
    log.null.fit <- null.logplatt.nig.fit(score.null.vec, alpha.init=1, beta.init=0.5, delta.init=1, mu.init=0, standardizeQ=standardizeQ, plotQ=F)
    
  } else if(distribution=="sn") {
    
    # Examine null fit for the skewed normal distribution:
    log.null.fit <- null.logplatt.sn.fit(score.null.vec, standardizeQ=standardizeQ, plotQ=F)
    
  } else if(distribution=="lg") {
    
    # Examine null fit for the normal distribution:
    log.null.fit <- null.logplatt.gaussian.fit(score.null.vec, standardizeQ = standardizeQ, plotQ = F)
    
  } else {
    stop("Bad specification of a parametric distribution for the Platt Null Scores!")
  }
  print("  **Done. See fit.diagnostics for fit information.")
  
  
  
  if(standardizeQ == TRUE) {
    #-----------------------------------------------------------------------------------------------------
    # Log and standardize the score.null
    #-----------------------------------------------------------------------------------------------------
    cat(sep="\n\n")
    print("Standardizing the Log of the Null Platt Scores...")
    log.s.null.mean <- mean(log(score.null.vec))                            # Mean of the log-null
    log.s.null.sd <- sd(log(score.null.vec))                                # SD of the log-null
    score.null.vec <- (log(score.null.vec) - log.s.null.mean)/log.s.null.sd # Standarize the boostrapped empirical log-null
    #hist(sln)
    print("  Done.")
    
  } else {
    #-----------------------------------------------------------------------------------------------------
    # Just Log the score.null
    #-----------------------------------------------------------------------------------------------------
    score.null.vec <- log(score.null.vec)
  }
  
  #*****NOTE: From here on the Platt score null, score.null.vec, has been log-ed and maybe (probably) standardized.

  
  #-------------------------------------------SCORE-NONNULL on the VALIDATION set-----------------------
  # Obtain the non null platt scores on the validation set.
  #-----------------------------------------------------------------------------------------------------
  cat(sep="\n\n")
  print("Computing Validation Set Log Non-Null Platt Scores...")
  score.nonnull.vec.val <- log(nonnull.logplatt.on.validation.set(training.dmat, training.labels, validation.dmat, validation.labels, svmtyp="C-classification", kern="linear", pparams=C.param))
  
  if(standardizeQ == TRUE) {
    #S tandardize the log non-null scores obtained from the validation set wrt the bootstrapped log null distribution:
    print("  Standardizing the Log of the Non-Null Platt Scores wrt Log Null Platt Scores...")
    score.nonnull.vec.val <- (score.nonnull.vec.val - log.s.null.mean)/log.s.null.sd    
  }
  print("  Done.")


  #-----------------------------------------SCORE-NULL on the VALIDATION set------------------------------
  #First obtain a random sample from the fit distributions, the same size as sln, to stand in for sln 
  #when computing p-values of the NULL and NONNULL
  #------------------------------------------------------------------------------------------------------

  cat(sep="\n\n")
  print("Beginning Simulation for Null Platt Scores Validation Set...")

  num.null.sims <- length(score.null.vec)
  
  if(distribution=="gev") {
    
  } else if (distribution=="nig") {
    
  } else if (distribution=="sn") {
    
  } else if (distribution=="lg") {
    
  }
  

  #**validation platt score log null sample from the GEV fit.
  #gev.log.s.null.vec <- rgev(num.null.sims, loc=ev.fit[[1]][1], scale=ev.fit[[1]][2], shape = ev.fit[[1]][3])

}