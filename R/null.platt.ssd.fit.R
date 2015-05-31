#---------------------------------------------------------------------
#Experiment to use a smoothing spline distribution to fit log null Platt scores.
#This fit is (1) SLOW!, (2) appears to be too tight aroun the histogram. IE, fit's support is only really
#around the histogram and dies off quick. Leading to 0 p-values at left tail and >1 p-values on the right tail...
#---------------------------------------------------------------------
null.platt.ssd.fit <- function(platt.null.vec, standardizeQ=FALSE, plotQ=FALSE) {
  
  stop("DON'T USE ME. Density fit too tight. Leads to probelm p-values at tails. 0 p-values on the left and >1 p-values on the right.")
  
  #Take the log:
  lgs <- log(platt.null.vec)
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  ssdfit <- ssdFit(lgs)
  
  #This is needed for both plots and fit diagnostics:
  lgs.hist.info <- hist(lgs,plot=F)
  
  if(plotQ==TRUE){
    
    print("Rendering diagnostic plots...")
    
    split.screen( figs = c( 1, 2 ) )
    
    screen(1)
    dens <- dssd(lgs, ssdfit)
    
    ylim.max <- max(dens,lgs.hist.info$density)
    xlim.max <- max(lgs,lgs.hist.info$breaks)
    xlim.min <- min(lgs,lgs.hist.info$breaks)
    
    if(standardizeQ==TRUE){
      fittitle <- "SSD fit to standardized log(null)"
    } else {
      fittitle <- "SSD fit to log(null)"
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
    Zt.hat<-qssd(tax, ssdfit)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for SSD fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for SSD fit to log(null) hist"
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
    
    print(paste("Bin:",i,
                "Upper:",upi,
                "Lower:",loi,
                "p-upper:",pssd(upi, ssdfit),
                "p-lower:",pssd(loi, ssdfit),
                "prob:", pssd(upi, ssdfit) - pssd(loi, ssdfit)
                )
          )
    
    fit.probs[i] <- pssd(upi, ssdfit) - pssd(loi, ssdfit)
    freq.expec[i] <- fit.probs[i] * length(platt.null.vec)
  }
  
  plt <- pssd(lgs.hist.info$breaks[1], ssdfit)
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  fit.info <- cbind(fit.probs, c(plt*length(platt.null.vec), freq.expec, prt*length(platt.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  #chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- ssdfit
  
  #return(list(fit.params, fit.info, chisq.results))
  return(list(fit.params, fit.info))
  
}