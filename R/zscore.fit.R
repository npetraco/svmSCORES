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
zscore.fit <- function(precomputed.null.scores = NULL, training.dmat, validation.dmat, training.labels, validation.labels,
                       distribution="gaussian", num.processes,
                       standardizeQ = TRUE, 
                       num.bs.iter = 2000,
                       C.param = 0.1,
                       printQ = FALSE,
                       plotQ = FALSE) {
  
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

  # We don't need score.nonnull anymore:
  remove(score.nonnull.vec)

  #-----------------------------------------------------------------------------------------------------
  # Parametric fits for the standardized log null (KNM) bootstrapped Platt scores
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
    # Log and standardize the score.null
    #-----------------------------------------------------------------------------------------------------
    cat(sep="\n\n")
    print("Standardizing the Log of the Bootstrapped Null Platt Scores...")
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
  print("**Computing VALIDATION SET LOG NON-NULL Platt Scores...")
  score.nonnull.vec.val <- log(nonnull.logplatt.on.validation.set(training.dmat, training.labels, validation.dmat, validation.labels, svmtyp="C-classification", kern="linear", pparams=C.param))
  
  if(standardizeQ == TRUE) {
    #S tandardize the log non-null scores obtained from the validation set wrt the bootstrapped log null distribution:
    print("  Standardizing the Log of the Non-Null Platt Scores wrt Log Null Platt Scores...")
    score.nonnull.vec.val <- (score.nonnull.vec.val - log.s.null.mean)/log.s.null.sd    
  }
  print("  Done.")


  #-----------------------------------------SCORE-NULL on the VALIDATION set------------------------------
  #First obtain a random sample from the fit distributions, the same size as the boostrapped null, to stand 
  #in for it when computing empirical p-values of the NULL and NONNULL
  #------------------------------------------------------------------------------------------------------

  cat(sep="\n\n")
  print("**Beginning Simulation for Null Platt Scores Validation Set...")
  
  print(" Sampling the Distribution Fit to Bootstrapped Log Null Platt Scores to Serve as a Model Empirical Score Null...")

  num.null.sims <- length(score.null.vec)
  
  #This process may take a few tries if bad (Infinite) z values are encountered. Put the process in a while loop
  #so it repeats a few times if need be. After it repeats a few times (currently 20), if bad p/z values are still
  #being encountered, stop. There is probably another problem. Through an error.
  continueQ <- FALSE
  max.val.calc.iters <- 20
  val.calc.iter <- 1
  while(continueQ == FALSE) {

    if(distribution=="gev") {
    
      score.null.vec.standin <- rgev(num.null.sims, loc=log.null.fit[[1]][1], scale=log.null.fit[[1]][2], shape = log.null.fit[[1]][3])
      
    } else if (distribution=="nig") {
      
      score.null.vec.standin <- rnig(num.null.sims, alpha=log.null.fit[[1]][1], beta=log.null.fit[[1]][2], delta=log.null.fit[[1]][3], mu=log.null.fit[[1]][4])
      
    } else if (distribution=="sn") {
      
      score.null.vec.standin <- rsn(num.null.sims, xi=log.null.fit[[1]][1], omega=log.null.fit[[1]][2], alpha=log.null.fit[[1]][3])
      
    } else if (distribution=="lg") {
      
      score.null.vec.standin <- rnorm(num.null.sims, mean=log.null.fit[[1]][1], sd=log.null.fit[[1]][2])
      
    }
    #???? ADD AN OPTION FOR A NEW EMPIRICAL SAMPLE ?????
    print("  Done.")
    
    print(paste("----",dist.nme,"Stand in Sample for the Bootstrapped Empirical Log Platt Score Null-----"))
    print(paste(" Number of Training (Fit Sampled) Null Platt Scores:", length(score.null.vec.standin)))
    print(paste(" Number of Bootstrapped Null Platt Scores:          ", length(score.null.vec)))
    print(paste(" Extent of Fit Sampled Null right tail:             ", max(score.null.vec.standin)))  
    print(paste(" Extent of Bootstrapped Empirical Null right tail:  ", max(score.null.vec)))
    print(paste(" Extent of Fit Sampled Null left tail:              ", min(score.null.vec.standin)))  
    print(paste(" Extent of Bootstrapped Empirical Null left tail:   ", min(score.null.vec)))
    
    #PAUSE FOR PLOT
    if(plotQ == TRUE) {
      fdrID::scores.histograms(score.null.vec, score.null.vec.standin, main="BS Null(Red), Sim-From-Fit Null (Green)")
      invisible(readline(prompt="Press [enter] to continue")) 
    }
    
    
    # Next generate a random sample from the fit distribution to serve as the validation set score-null. 
    cat(sep="\n\n")
    print(" Sampling Distribution Fit to Bootstrapped Log Null Platt Scores to Serve As the VALIDATION SET LOG NULL Platt Scores...")
    
    num.null.sims <- nrow(validation.dmat) * (nlevels(validation.labels)-1) # There are this many null values for a validation set
    
    if(distribution=="gev") {
      
      score.null.vec.val <- rgev(num.null.sims, loc=log.null.fit[[1]][1], scale=log.null.fit[[1]][2], shape = log.null.fit[[1]][3])
      
    } else if (distribution=="nig") {
      
      score.null.vec.val <- rnig(num.null.sims, alpha=log.null.fit[[1]][1], beta=log.null.fit[[1]][2], delta=log.null.fit[[1]][3], mu=log.null.fit[[1]][4])
      
    } else if (distribution=="sn") {
      
      score.null.vec.val <- rsn(num.null.sims, xi=log.null.fit[[1]][1], omega=log.null.fit[[1]][2], alpha=log.null.fit[[1]][3])
      
    } else if (distribution=="lg") {
      
      score.null.vec.val <- rnorm(num.null.sims, mean=log.null.fit[[1]][1], sd=log.null.fit[[1]][2])
      
    }
    print("  Done.")
    
    print(paste(" ----Characteristics of Sample From",dist.nme,"Fit to Bootstrapped Log Platt Score Null.-----"))
    print(paste(" Number of Fit Validation Set Null Platt Scores:      ", length(score.null.vec.val)))
    print(paste(" Extent of Fit Validation Set Null Sample right tail: ", max(score.null.vec.val)))  
    print(paste(" Extent of Bootstrapped Empirical Null right tail:    ", max(score.null.vec)))
    print(paste(" Extent of Fit Validation Set Null Sampled left tail: ", min(score.null.vec.val)))  
    print(paste(" Extent of Bootstrapped Empirical Null left tail:     ", min(score.null.vec)))
    
    
    #PAUSE FOR PLOT
    #fdrID::scores.histograms(sln, gev.log.s.null.vec.val, ylim.max=0.45, xlim.min = -4.5, xlim.max = 4)
    if(plotQ == TRUE) {
      fdrID::scores.histograms(score.null.vec, score.null.vec.val, main="BS Null(Red), Validation Set Null (Green)")
      invisible(readline(prompt="Press [enter] to continue")) 
    }
    
    
    #Use the above null score stand in right away to get p and z from the validation sets and see if there are 
    #any problems.
    
    cat(sep="\n\n")
    print(" Computing z-values on VALIDATION SET LOG NULL Platt Scores...")
    
    #VALIDATION SET Null p-values and z-values wrt substituting the fit sampled (log) null as emiprical null:
    p.null <- empiricalPvalues(score.null.vec.standin, score.null.vec.val)
    z.null <- qnorm(p.null)
    
    #VALIDATION SET NONNULL p-values wrt substituting the fit sampled (log) null as emiprical null:
    p.nonnull <- empiricalPvalues(score.null.vec.standin, score.nonnull.vec.val)
    #
    #"Smear" nonnull zero p-values wrt null used. Return POTENTIAL z values. These still may need some tweeking after
    #this routine finshes.  
    print(" Test smearing negative infinity non-null z-values...")
    z.nonnull.smeared <- smear.extreme.nonnull.zvalues(p.nonnull, 
                                                       upper.set.zvalue = (-12), 
                                                       mu.factor = (-0.5), 
                                                       p.factor=0.99, 
                                                       printQ=printQ, plotQ=F)[[3]]
    #print(z.nonnull.smeared)
    print(" Done.")
    
    
    #If problem p and z values are encountered at this NULL construction stage, stop and go get new 
    #randomsamples at ll.79 and/or ll.99
    
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
  }
  print(" **Acceptable Null Log Platt Score Validation Set Obtained.")
  
  #Before moving on, check null and nonnull pvalues for problems.
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
  cat(sep="\n\n")
  print("NOTE: All p-values computed are empirical wrt the Training Log Null Simulated Platt Scores.")
  print("      Integration of the Fit PDF was not used because we have noticed numerical instabilities for very small/large p-values.")
  cat(sep="\n\n")
  print("**CAUTION double check z-values graphically and with check.ps.and.zs function,")
  print("          especially if you used a GEV fit to the Bootstrap Null!!!!!!!!!!!!!!")  
  print("=============================================================================")
    
  all.score.calc.info <- list(
    standardization.used=standardizeQ, # Flag if standardization was used.
    fit.distribution.name=distribution,    # Record of parametric distribution name use for fit and Null simulations
    log.null.fit,                          # Fit info/diagnostics for chosen parametric fit to Boostrapped Log Null Platt Scores
    score.null.vec,                        # Boostrapped Null Platt Scores
    score.null.vec.standin,                # Simulated stand-in for boostrapped Log Null Platt Scores  
    score.null.vec.val,                    # Validation Log NonNull Platt Scores  
    score.nonnull.vec.val,                 # Validation Log Null Platt Scores  
    p.null,                                # Null p-values from score.null.vec.val, and wrt simulated stand-in for boostrapped Log Null Platt Scores  
    p.nonnull,                             # Un-smeared non-null p-values from score.nonnull.vec.val, and wrt simulated stand-in for boostrapped Log Null Platt Scores
    z.null,                                # Null z-values from Null p-values.
    z.nonnull.smeared)                     # Potential (smeared) Non Null z-values. Check acceptability especially if gev fit to bs-null used.
                                           # **No non-null un-smeard z-values are returned because they need to be modified. See p.nonnull if you want them
  
  
  names(all.score.calc.info) <-c(
    "standardization.flag",
    "fit.distribution.name",
    "fit.info.and.diagnostics",
    "boostrapped.null.platt.scores",
    "score.null.training",
    "score.null.validation",
    "score.nonnull.validation",
    "null.pvalues",
    "nonnull.pvalues", #UN-SMEARED!!!!!!
    "null.zvalues",
    "nonnull.zvalues.smeared" #SMEARED!!!!!!
    )
  
  return(all.score.calc.info)
  
}