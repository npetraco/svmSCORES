#--------------------------------------------
#' @title null.logplatt.sn.fit
#' 
#' @description Skewed normal fit for log null Platt scores.
#' 
#' @details XXXX 
#' 
#' @param XXXX
#' 
#' @return XXXX
#'
#' @references Storey JD, Tibshirani R. Statistical significance for genomewide studies. PNAS 100(16):9440-9445 (2003)
#'
#' @examples
#' XXXX
#--------------------------------------------
null.logplatt.sn.fit <- function(platt.null.vec, standardizeQ=FALSE, plotQ=FALSE) {
  
  #Take the log:
  lgs <- log(platt.null.vec)
  
  if(standardizeQ==TRUE){
    lgs <- (lgs - mean(lgs))/sd(lgs)
  }
  
  lgs <- sort(lgs)
  
  sn.fit <- selm(lgs~1)
  dp.vec <- coef(sn.fit, param.type = "DP")
  xi.est <- dp.vec[1]
  om.est <- dp.vec[2]
  alp.est <- dp.vec[3]
  
  
  #This is needed for both plots and fit diagnostics:
  lgs.hist.info <- hist(lgs,plot=F)
  
  if(plotQ==TRUE){
    
    print("Rendering diagnostic plots...")
    
    split.screen( figs = c( 1, 2 ) )
    
    screen(1)
    dens <- dsn(lgs, xi=xi.est, omega=om.est, alpha=alp.est, tau=0)
    
    ylim.max <- max(dens,lgs.hist.info$density)
    xlim.max <- max(lgs,lgs.hist.info$breaks)
    xlim.min <- min(lgs,lgs.hist.info$breaks)
    
    if(standardizeQ==TRUE){
      fittitle <- "SN fit to standardized log(null)"
    } else {
      fittitle <- "SN fit to log(null)"
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
    Zt.hat <- qsn(tax, xi=xi.est, omega=om.est, alpha=alp.est, tau=0)
    
    if(standardizeQ==TRUE){
      QQtitle <- "Q-Q plot for SN fit to standardized log(null) hist"
    } else {
      QQtitle <- "Q-Q plot for SN fit to log(null) hist"
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
    
    fit.probs[i] <- psn(upi, xi=xi.est, omega=om.est, alpha=alp.est, tau=0) - psn(loi, xi=xi.est, omega=om.est, alpha=alp.est, tau=0) 
    freq.expec[i] <- fit.probs[i] * length(platt.null.vec)
  }
  
  plt <- psn(lgs.hist.info$breaks[1], xi=xi.est, omega=om.est, alpha=alp.est, tau=0)
  prt <- 1-sum(c(plt,fit.probs))
  fit.probs <- c(plt,fit.probs,prt)
  
  fit.info <- cbind(fit.probs, c(plt*length(platt.null.vec), freq.expec, prt*length(platt.null.vec)), c(0,freq.obs,0))
  colnames(fit.info) <- c("interquant.probs", "interquant.exp.cts", "interquant.obs.cts")
  chisq.results <- chisq.test(c(0,freq.obs,0), p = fit.probs)
  
  fit.params <- c(xi.est, om.est, alp.est)
  names(fit.params) <- c("xi.hat", "omega.hat", "alpha.hat")
    
  return(list(fit.params, fit.info, chisq.results))
  
}