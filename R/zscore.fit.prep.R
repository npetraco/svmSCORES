zscore.fit.prep <- function(precomputed.null.scores = NULL, 
                       training.dmat, validation.dmat, 
                       training.labels, validation.labels,
                       distribution="gaussian", 
                       num.processes,
                       standardizeQ = TRUE, 
                       num.bs.iter = 2000,
                       C.param = 0.1,
                       printQ = FALSE,
                       plotQ = FALSE) {
  
  print("#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#")
  print("                                      Prep For Z-Score Construction                                    ")
  print("#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#")
  cat(sep="\n\n")
  
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
  
  #---------------------Compute a bootstrap set of Null (and throw-away NonNull) Platt Scores if requested:----------
  print("Beginning Bootstrap Computation of the Null Platt Scores")
  if(is.null(precomputed.null.scores)) {
    
    platt.scores <- bootstrap.platt.scores.parallel(
      num.processes=num.processes, 
      dat.mat=training.dmat, 
      lbls=training.labels, 
      nbs=num.bs.iter, 
      svmtyp="C-classification", 
      kern="linear", 
      pparams=C.param,
      timerQ=T)
    
  } else {
    platt.scores <- precomputed.null.scores 
  }
  print(" Done.")
  
  print("========================================================")
  
  score.null.vec<-platt.scores[,1]
  score.nonnull.vec<-platt.scores[,2]
  
  print(paste("Number of Bootstrapped Null Platt Scores:", length(score.null.vec)))
  print(paste("Extent of Null right tail:               ", max(score.null.vec)))
  print(paste("Extent of Non-Null left tail:            ", min(score.nonnull.vec)))
  
  #PAUSE FOR PLOT HERE
  if(plotQ == TRUE) {
    fdrID::scores.histograms(score.null.vec, score.nonnull.vec)
    invisible(readline(prompt="Press [enter] to continue")) 
  }
  #
  
  # We don't need the boostrapped score.nonnull.vec anymore:
  remove(score.nonnull.vec)
  
  #-----------------------------------------------------------------------------------------------------
  # Parametric fits for the log null (KNM) bootstrapped (training) Platt scores
  #-----------------------------------------------------------------------------------------------------
  cat(sep="\n\n")
  print(paste("Fitting log null Platt score bootstrap sample to a", dist.nme, "distribution..."))
  if(distribution=="gev") {
    
    # Examine null fit for the extreme value distribution:    
    if(plotQ == TRUE) {
      #PAUSE FOR PLOT
      log.null.fit <- null.logplatt.gevd.fit(score.null.vec, standardizeQ=standardizeQ, plotQ=T)
      invisible(readline(prompt="Press [enter] to continue")) 
    } else {
      #No plot
      log.null.fit <- null.logplatt.gevd.fit(score.null.vec, standardizeQ=standardizeQ, plotQ=F)
    }
    
  } else if(distribution=="nig") {
    
    # Examine null fit for the normal inverse gaussian distribution:
    if(plotQ == TRUE) {
      #PAUSE FOR PLOT
      log.null.fit <- null.logplatt.nig.fit(score.null.vec, alpha.init=1, beta.init=0.5, delta.init=1, mu.init=0, standardizeQ=standardizeQ, plotQ=T)
      invisible(readline(prompt="Press [enter] to continue")) 
    } else {
      #No plot
      log.null.fit <- null.logplatt.nig.fit(score.null.vec, alpha.init=1, beta.init=0.5, delta.init=1, mu.init=0, standardizeQ=standardizeQ, plotQ=F)
    }
    
  } else if(distribution=="sn") {
    
    # Examine null fit for the skewed normal distribution:
    if(plotQ == TRUE) {
      #PAUSE FOR PLOT
      log.null.fit <- null.logplatt.sn.fit(score.null.vec, standardizeQ=standardizeQ, plotQ=T)
      invisible(readline(prompt="Press [enter] to continue")) 
    } else {
      #No plot
      log.null.fit <- null.logplatt.sn.fit(score.null.vec, standardizeQ=standardizeQ, plotQ=F)
    }
    
  } else if(distribution=="lg") {
    
    # Examine null fit for the normal distribution:
    if(plotQ == TRUE) {
      #PAUSE FOR PLOT
      log.null.fit <- null.logplatt.gaussian.fit(score.null.vec, standardizeQ = standardizeQ, plotQ = T)
      invisible(readline(prompt="Press [enter] to continue")) 
    } else {
      #No plot
      log.null.fit <- null.logplatt.gaussian.fit(score.null.vec, standardizeQ = standardizeQ, plotQ = F)
    }
    
  } else {
    stop("Bad specification of a parametric distribution for the Platt Null Scores!")
  }
  print("  **Done. See fit.diagnostics for fit information.")
  
  
  
  if(standardizeQ == TRUE) {
    #-----------------------------------------------------------------------------------------------------
    # Log and standardize the training set (bootstrapped) log score.null
    #-----------------------------------------------------------------------------------------------------
    cat(sep="\n\n")
    print("Standardizing the Log of the Bootstrapped Null Platt Scores...")
    log.s.null.mean <- mean(log(score.null.vec))                            # Mean of the log-null
    log.s.null.sd <- sd(log(score.null.vec))                                # SD of the log-null
    score.null.vec <- (log(score.null.vec) - log.s.null.mean)/log.s.null.sd # Standarize the boostrapped empirical log-null
    
    print("  Done.")
    
  } else {
    #-----------------------------------------------------------------------------------------------------
    # Just Log the training set (bootstrapped) log score.null
    #-----------------------------------------------------------------------------------------------------
    score.null.vec <- log(score.null.vec)
  }
  
  #*****NOTE: From here on the Platt score null, score.null.vec, has been log-ed and maybe (probably) standardized.
  
  
  #-------------------------------------------SCORE-NONNULL on the VALIDATION set-----------------------
  # Obtain the non null platt scores on the validation set of feature vectors.
  #-----------------------------------------------------------------------------------------------------
  cat(sep="\n\n")
  print("**Computing VALIDATION SET LOG NON-NULL Platt Scores...")
  score.nonnull.vec.val <- log(nonnull.logplatt.on.validation.set(training.dmat, training.labels, validation.dmat, validation.labels, svmtyp="C-classification", kern="linear", pparams=C.param))
  
  if(standardizeQ == TRUE) {
    #Standardize the log non-null scores obtained from the validation set wrt the bootstrapped log null distribution:
    print("  Standardizing the Log of the Non-Null Platt Scores wrt Log of the boostrapped Null Platt Scores...")
    score.nonnull.vec.val <- (score.nonnull.vec.val - log.s.null.mean)/log.s.null.sd    
  }
  print("  Done.")
  
  score.calc.prep.info <- list(
    standardization.used <- standardizeQ, # Flag if standardization was used.
    fit.distribution.name <- distribution,    # Record of parametric distribution name use for fit and Null simulations
    log.null.fit,                          # Fit info/diagnostics for chosen parametric fit to Boostrapped Log Null Platt Scores
    score.null.vec,                        # Boostrapped Null Platt Scores
    score.nonnull.vec.val)                 # Validation Log Null Platt Scores 
  
  
  names(score.calc.prep.info) <-c(
    "standardization.flag",
    "fit.distribution.name",
    "fit.info.and.diagnostics",
    "boostrapped.log.null.platt.scores",
    "log.nonnull.platt.scores")
  
  return(score.calc.prep.info)
  
  
}