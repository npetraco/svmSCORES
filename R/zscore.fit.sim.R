#--------------------------------------------------------------------------
#Internal function. The simulation portions of the z-score fitting process
#--------------------------------------------------------------------------
zscore.fit.sim <- function(score.null.vec, distribution, log.null.fit, 
                           nrow.validation.dmat, nlevels.validation.labels, 
                           printQ=FALSE, plotQ=FALSE) {
  
  #-----------------------------------------SCORE-NULL on the VALIDATION set------------------------------
  # First obtain a random sample from the fit distributions, the same size as the boostrapped null, to stand 
  # in for the boostrapped training null when computing empirical p-values of the NULL and NONNULL. This new set
  # of simulated log null Platt scores (drawn from the fit distribution) will serve as the training log null
  # Platt scores.
  #------------------------------------------------------------------------------------------------------

  print("#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#")
  print("                     Simulations for Training and Validation Set Log Null Platt Scores                 ")
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
  
  print("**Beginning Simulation for Null Platt Scores Validation Set...")
  
  print(" Sampling the Distribution Fit to Bootstrapped Log Null Platt Scores to Serve as a Model Empirical Score Null...")
  
  num.null.sims <- length(score.null.vec)
  
  #This process may take a few tries if bad (Infinite) z values are encountered. Put the process in a while loop
  #so it repeats a few times if need be. After it repeats a few times (currently 20), if bad p/z values are still
  #being encountered, stop. There is probably another problem. Throw an error.
  #
  #NOTE: score.null.vec.standin will serve as the "training set" log null Platt socres.
  
  if(distribution=="gev") {
    
    score.null.vec.standin <- rgev(num.null.sims, loc=log.null.fit[[1]][1], scale=log.null.fit[[1]][2], shape = log.null.fit[[1]][3])
    
  } else if (distribution=="nig") {
    
    score.null.vec.standin <- rnig(num.null.sims, alpha=log.null.fit[[1]][1], beta=log.null.fit[[1]][2], delta=log.null.fit[[1]][3], mu=log.null.fit[[1]][4])
    
  } else if (distribution=="sn") {
    
    score.null.vec.standin <- rsn(num.null.sims, xi=log.null.fit[[1]][1], omega=log.null.fit[[1]][2], alpha=log.null.fit[[1]][3])
    
  } else if (distribution=="lg") {
    
    score.null.vec.standin <- rnorm(num.null.sims, mean=log.null.fit[[1]][1], sd=log.null.fit[[1]][2])
    
  }
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
  #
  #NOTE: score.null.vec.val will serve as the "validation set" log null Platt socres.
  cat(sep="\n\n")
  print(" Sampling Distribution Fit to Bootstrapped Log Null Platt Scores to Serve As the VALIDATION SET LOG NULL Platt Scores...")
  
  # There are this many null values for a validation set:
  num.null.sims <- nrow.validation.dmat * nlevels.validation.labels-1
  
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
  if(plotQ == TRUE) {
    fdrID::scores.histograms(score.null.vec, score.null.vec.val, main="BS Null(Red), Validation Set Null (Green)")
    invisible(readline(prompt="Press [enter] to continue")) 
  }
  
  #Use the above null score stand in right away to get p and z from the validation sets and see if there are 
  #any problems.

  sim.score.calc.info <- list(
    fit.distribution.name=distribution,    # Record of parametric distribution name use for fit and Null simulations
    score.null.vec.standin,                # Simulated stand-in for boostrapped Log Null Platt Scores  
    score.null.vec.val)                    # Validation Log NonNull Platt Scores


  names(sim.score.calc.info) <-c(
    "fit.distribution.name",
    "score.null.training",
    "score.null.validation"
  )
  
  return(sim.score.calc.info)

}