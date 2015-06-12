#--------------------------------------------
#' @title Fit log of null Platt scores to a Gaussian
#' 
#' @description Gaussian fit for log null Platt scores.
#' 
#' @details Simple wrapper to get the Gaussian fit around the log of the null Platt scores.
#'
#' @param platt.null.vec     Vector of null (non-match) Platt scores from an SVM
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
null.platt.gaussian.fit <- function(platt.null.vec, standardizeQ=FALSE, plotQ=FALSE) {
    
  #Take the log:
  lgs <- log(platt.null.vec)
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  normfit <- MASS::fitdistr(lgs,"normal")
  #If data was standardized, these should be pretty close to 0 and 1:
  mu <- normfit$estimate[1]
  sig <- normfit$estimate[2]

  
  #This is needed for both plots and fit diagnostics:
  lgs.hist.info <- hist(lgs,plot=F)
  
  if(plotQ==TRUE){
    
    print("Rendering diagnostic plots...")
    
    split.screen( figs = c( 1, 2 ) )
    
    screen(1)
    dens <- dnorm(lgs, mean=mu, sd=sig)
    
    ylim.max <- max(dens,lgs.hist.info$density)
    xlim.max <- max(lgs,lgs.hist.info$breaks)
    xlim.min <- min(lgs,lgs.hist.info$breaks)
    
    if(standardizeQ==TRUE){
      fittitle <- "Gaussian fit to standardized log(null)"
    } else {
      fittitle <- "Gaussian fit to log(null)"
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
    Zt.hat<-qnorm(tax, mean=mu, sd=sig)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for Gaussian fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for Gaussian fit to log(null) hist"
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
    
#    print(paste("Bin:",i,
#                "Upper:",upi,
#                "Lower:",loi,
#                "p-upper:",pssd(upi, ssdfit),
#                "p-lower:",pssd(loi, ssdfit),
#                "prob:", pssd(upi, ssdfit) - pssd(loi, ssdfit)
#    )
#    )
    
    fit.probs[i] <- pnorm(upi, mean=mu, sd=sig) - pnorm(loi, mean=mu, sd=sig)
    freq.expec[i] <- fit.probs[i] * length(platt.null.vec)
  }
  
  plt <- pnorm(lgs.hist.info$breaks[1], mean=mu, sd=sig)
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  fit.info <- cbind(fit.probs, c(plt*length(platt.null.vec), freq.expec, prt*length(platt.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- c(mu,sig)
  names(fit.params) <- c("mu.est","sig.est")
  
  return(list(fit.params, fit.info, chisq.results))
  #return(list(fit.params, fit.info))
  
}