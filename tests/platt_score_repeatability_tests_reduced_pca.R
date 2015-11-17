
#***********Think about using avg 1-vs-1 dec score of pred and to replave platts. 
#The platt scores seem to not be gauranteed repatable past two decimal places...........
#See below

library(fdrID)
library(svmSCORES)
library(doMC)
source("/Users/npetraco/codes/R/profiles/sourceme.R")
source("/Users/npetraco/codes/R/error_rate_utilities/sourceme.R")

fv.mat <- read.binary.profiles("/Users/npetraco/latex/papers/posterior_error/comparison_examples/NIST_Martin_screwdriver/NFIToolmark_Biasotti_Murdock_features_15-30_44-SUBSET_3DX3P.bin")

pca.model<-prcomp(fv.mat,center=TRUE,scale=FALSE)
plot(pca.model)
summary(pca.model)

M<-45                                                         #Pick dimension
Z<-predict(pca.model)[,1:M]                                   #Grab PCA scores
lbl <- gl(44,30)

#Training set:
trainobs<-createDataPartition(lbl,p=(1/3),list=F)
#count.group.replicates(lbl[trainobs])
Ztr<-Z[trainobs,]
lbltr<-lbl[trainobs]
#Tmp to split val and test:
Zvalte<-Z[-trainobs,]
lblvalte<-lbl[-trainobs]
#Validation set:                                    
valobs<-createDataPartition(lblvalte,p=0.5,list=F)
Zval<-Zvalte[valobs,]
lblval<-lblvalte[valobs]
#Test Set:
Zte<-Zvalte[-valobs,]
lblte<-lblvalte[-valobs]

Z <- Ztr
lbl <- lbltr
rownames(Z)<-NULL
colnames(Z)<-NULL
pen.param<-0.1
ind.vec<-NULL
for(i in 1:nrow(Z))
{
  #i <- 1
  Z.heldout<-t(as.matrix(Z[i,]))
  lbl.heldout<-lbl[i]
  
  Z.kept<-Z[-i,]
  lbl.kept<-lbl[-i]
  svm.model<-svm(Z.kept,lbl.kept,scale=FALSE,type="C-classification",kernel="linear",cost=pen.param,fitted=F,probability=F)
  pred<-predict(svm.model,Z.heldout)
  #print(pred==lbl.heldout)
  #  prob.vec<-attr(pred, "probabilities")[,]
  ind.vec<-c(ind.vec,pred==lbl.heldout)
  curr.err.perc <- 100 - (sum(ind.vec)/length(ind.vec) * 100)
  print(paste("Done iter:",i, "Correct??:","Pred:",pred,"Actual:",lbl.heldout,pred==lbl.heldout, ". Current err %:", curr.err.perc))
  
}
(1-sum(ind.vec)/nrow(Z))*100
#~2.3% error on Ztr at 30D

#Get a set of bootstrapped null Platt scores:
pscs <- bootstrap.platt.scores.parallel(
  num.processes=8, 
  dat.mat=Ztr, 
  lbls=lbltr, 
  nbs=2000, 
  svmtyp="C-classification", 
  kern="linear", 
  pparams=0.1,
  timerQ=T)

#Construct a set of z-scores:
zscs.info <- zscore.fit2(pscs, Ztr, Zval, lbltr, lblval, 
                         pvalue.method = "integral",
                         distribution="lg",
                         num.processes=8, 
                         standardizeQ=T, 
                         num.bs.iter=2000, 
                         C.param = 0.1, 
                         printQ=F, plotQ=F)

#Checks:
names(zscs.info)
zscs.info$fit.info.and.diagnostics
z.null <- zscs.info$null.zvalues
z.nonnull <- zscs.info$nonnull.zvalues.smeared
p.nonnull <- zscs.info$nonnull.pvalues

# z.nonnull.pot <- smear.extreme.nonnull.zvalues(
#   zscs.info$nonnull.pvalues,
#   upper.set.zvalue = (-12), 
#   mu.factor = (-3.8), 
#   p.factor=0.9999, plotQ=T)[[3]]
# z.nonnull <- z.nonnull.pot

cbind(p.nonnull,z.nonnull)
check.ps.and.zs(null.p.values = pnorm(z.null), nonnull.p.values = pnorm(z.nonnull), printQ = T, plotQ = T)
hist(z.nonnull, main="Validation z Non-Null")
hist(z.null, main="Validation z Null")
hist(c(z.null,z.nonnull), main="All z")

#Fit lfdr models:
p.vals <- c(pnorm(z.null), pnorm(z.nonnull))

#First Do Efron lfdr fit to examine fit to f(z)
library(locfdr)
fdr.model<-locfdr(qnorm(p.vals), bre = 120, df = 10, pct = 0, pct0 = 1/4, nulltype = 1, type =0, plot = 1, main = " ", sw = 0)

minfo <- sampler.prep(p.vals, num.bins=120, degree=10, interceptQ=T, overdispersionQ=F, sampler="jags")
jsim <- jags(data=minfo$Data, inits=minfo$Initialization.Function, minfo$Model.Parameters, n.iter=10000, n.chains=4, model.file=minfo$BUG.Model.File.Path)
jsim

posterior.f <- jsim$BUGSoutput$sims.list$lambda
x <- minfo$Bin.Midpoints
z <- minfo$z.Values

lfdr.info <- make.fdr.functions(z, x, pct0=0.25, posterior.f, credibility.level=0.95, interval.type="hpd", plotQ=T)
efdrf <- lfdr.info$fdr.means.func
med.fdrf <- lfdr.info$fdr.means.func
uci.fdrf <- lfdr.info$upper.fdr.ci.func
lci.fdrf <- lfdr.info$lower.fdr.ci.func
etdrf <- lfdr.info$tdr.means.func
med.tdrf <- lfdr.info$tdr.means.func
uci.tdrf <- lfdr.info$upper.tdr.ci.func
lci.tdrf <- lfdr.info$lower.tdr.ci.func

p0 <- lfdr.info$p0
d0 <- lfdr.info$delta0
s0 <- lfdr.info$sigma0

hist(p0)
mean(p0)

#-------------------------
# Test input Zte and output Platts
mln <- mean(log(pscs[,1]))
sln <- sd(log(pscs[,1]))
tinfo3 <- posterior.probs.for.unks2(training.dmat=Ztr, training.labels=lbltr, C.param=0.1, test.dmat = Zte,
                                    standardizeQ=zscs.info$standardization.flag, 
                                    mean.log.null.score=mln, sd.log.null.score=sln,
                                    pvalue.method = "integral",
                                    null.vec.training = zscs.info$score.null.training,
                                    distribution.name = zscs.info$fit.distribution.name,
                                    distribution.fit.info = zscs.info$fit.info.and.diagnostics,
                                    pp.point.est.func=efdrf,
                                    pp.uci.est.func=uci.fdrf,
                                    pp.lci.est.func=lci.fdrf)
tinfo3
plot(tinfo3[,4],typ="h")

z.unknowns <- tinfo3[,2]
z.unknowns
zoomed.post.prob.plot(c(z.null,z.nonnull), zbounds=c(NA, max(z.unknowns)), 
                      point.est.func=efdrf, 
                      upper.est.func=uci.fdrf, 
                      lower.est.func=lci.fdrf, prob.scale="percent", xlab=NULL, ylab="fdr (%)", main=NULL)
min(z.unknowns)
min(c(z.null,z.nonnull))
max(z.unknowns)
max(c(z.null,z.nonnull))

points(z.unknowns[which(100*efdrf(z.unknowns) <=5)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) <=5)],col="green")
points(z.unknowns[which(100*efdrf(z.unknowns) > 5 & 100*efdrf(z.unknowns) <= 50)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) > 5 & 100*efdrf(z.unknowns) <= 50)],col="yellow")
points(z.unknowns[which(100*efdrf(z.unknowns) > 50)], 100*efdrf(z.unknowns)[which(100*efdrf(z.unknowns) > 50)],col="red")
#Got these one wrong
wrong.idxs <- which((tinfo3[,1] == lblte)==FALSE)
data.frame(lblte[wrong.idxs],tinfo3[wrong.idxs,])
points(z.unknowns[wrong.idxs], 100*efdrf(z.unknowns)[wrong.idxs],col="black", pch=16)

length(wrong.idxs)/nrow(Zte) * 100
