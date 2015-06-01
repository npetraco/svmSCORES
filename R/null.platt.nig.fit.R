#--------------------------------------------
#' @title null.platt.nig.fit
#' 
#' @description Normal inverse Gaussian fit for log null Platt scores.
#' 
#' @details We noticed that the log of the null (non-match) Platt scores typically looks close to Gaussian.
#' However we've also noticed that this distribution can show some skew. The normal inverse Gaussian (nig) 
#' distribution is a flexible four parameter distribution which usually gives a good fit to the null log 
#' Platt score distribution estimated by the bootstrapping procedure. This routine calls the \code{nigFit} function
#' in "mle" mode from the fBasics package.
#' 
#' It is not necessary to fit the log of the null (non-match) Platt scores as p-values with respect to them can 
#' be obtained empirically (cf. Storey and Tibshirani). Also other fits are possible. It is up to the user if 
#' they choose to use this function to help obtain p-values. 
#' 
#' In practice, we use this function because (1) it avoids p-values = 0 (always encountered when obtaining empirical p-values 
#' for non-null (known-matching) scores), (2) it is faster than computing empirical p-values from hundres of 
#' thousands or millions of bootstrapped null scores and (3) it cuts down on the spagetti code for the 
#' distribution's fit when we do it step-by-step.
#'
#' @param platt.null.vec     Vector of null (non-match) Platt scores from an SVM
#' @param alpha.init         Initial guess of alpha parameter for NIG
#' @param beta.init          Initial guess of beta parameter for NIG
#' @param delta.init         Initial guess of delta parameter for NIG
#' @param mu.init            Initial guess of mu parameter for NIG
#' @param standardizeQ       Whether or not to standardize the log null scores
#' @param plotQ              Diagnostic plots?  
#' 
#' @return A list with the fitted parameters, fit info and chi-square goodness of fit test results
#'
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
null.platt.nig.fit <- function(platt.null.vec, alpha.init=NULL, beta.init=NULL, delta.init=NULL, mu.init=NULL, standardizeQ=FALSE, plotQ=FALSE) {
  
  #Take the log:
  lgs <- log(platt.null.vec)
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  nigFit.lgs <- fBasics::nigFit(lgs, alpha = alpha.init, beta = beta.init, delta = delta.init, mu = mu.init, doplot = FALSE, method="mle")
  
  alp.est <- nigFit.lgs@fit$estimate[1]
  bet.est <- nigFit.lgs@fit$estimate[2]
  del.est <- nigFit.lgs@fit$estimate[3]
  mu.est <- nigFit.lgs@fit$estimate[4]
  
  #This is needed for both plots and fit diagnostics:
  lgs.hist.info <- hist(lgs,plot=F)
  
  if(plotQ==TRUE){
    
    print("Rendering diagnostic plots...")
    
    split.screen( figs = c( 1, 2 ) )
    
    screen(1)
    dens <- fBasics::dnig(lgs, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)
    
    ylim.max <- max(dens,lgs.hist.info$density)
    xlim.max <- max(lgs,lgs.hist.info$breaks)
    xlim.min <- min(lgs,lgs.hist.info$breaks)
    
    if(standardizeQ==TRUE){
      fittitle <- "NIG fit to standardized log(null)"
    } else {
      fittitle <- "NIG fit to log(null)"
    }
    
    plot(lgs, dens, xlim=c(xlim.min,xlim.max), ylim=c(0,ylim.max), typ="l", col="blue", lwd=3, xlab="", ylab="")
    par(new=T)
    hist(lgs, probability=T, xlim=c(xlim.min,xlim.max), ylim=c(0,ylim.max), xlab="STD log(KNM)", main=fittitle)
    
    #Q-Q Plot:
    #t-axis:
    tmax<-1000
    tax<-seq(1,tmax,1)/(tmax+1)
    cemp<-ecdf(lgs)(lgs)
    #Empirical quantile function (inverse CDF):
    qemp<-splinefun(cemp,lgs)
    #Empirical quantiles:
    Zt<-qemp(tax)
    #Quantiles from fit:
    Zt.hat <- fBasics::qnig(tax,alpha=alp.est,beta=bet.est,delta=del.est,mu=mu.est)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for NIG fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for NIG fit to log(null) hist"
    }
    
    #Q-Q plot:
    screen(2)
    plot(Zt,Zt.hat, xlab="empirical quantiles", ylab="fit quantiles",main=QQtitle)
    abline(0,1)
    
    close.screen( all = TRUE )
    
  }
  
  freq.obs <- lgs.hist.info$counts
  freq.expec <- rep(-1,length(lgs.hist.info$mids))
  fit.probs <- rep(-1,length(lgs.hist.info$mids))
  
  print("Computing fit interquantile probabilities...")
  for(i in 1:(length(lgs.hist.info$breaks)-1)){
    
    upi <- lgs.hist.info$breaks[i+1]
    loi <- lgs.hist.info$breaks[i]
    
    fit.probs[i] <- fBasics::pnig(upi, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est) - fBasics::pnig(loi, alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)
    freq.expec[i] <- fit.probs[i] * length(platt.null.vec)
  }
  
  plt <- fBasics::pnig(lgs.hist.info$breaks[1], alpha=alp.est, beta=bet.est, delta=del.est, mu=mu.est)
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  fit.info <- cbind(fit.probs, c(plt*length(platt.null.vec), freq.expec, prt*length(platt.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- c(alp.est, bet.est, del.est, mu.est)
  names(fit.params) <- c("alpha.hat", "beta.hat", "delta.hat", "mu.hat")
  
  return(list(fit.params, fit.info, chisq.results))
  
}