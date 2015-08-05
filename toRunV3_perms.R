# R CMD BATCH "--args phenotype=\"$p\" "  toRunV3.R

args=(commandArgs(TRUE))
for(i in 1:length(args))
{
  eval(parse(text=args[[i]]))
}

rm(list=ls())
load("eqtl2015_finalbinnormCrosses.RData")

pkg <- c("rrBLUP","RCurl","plyr","qtl")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))

function.names<-c("pipeV3.R")
us<-paste("https://raw.githubusercontent.com/jtlovell/eQTL_functions/master/",function.names, sep="")
for(i in 1:length(function.names)){
  cat("loading",function.names[i],"\n")
  script <- getURL(us[i], ssl.verifypeer = FALSE)
  eval(parse(text = script))
}

transformExpression<-function(x) {
  x<-as.numeric(x)
  if(sum(x[!is.na(x)]==0)/length(x[!is.na(x)]) >= 0.05){
    log2(x+1) ######### If there are more than 5% 0's in the data, use the 0's
  }else{
    x[x==0]<-NA  ######### If there are less than 5% 0's in the data, toss the 0's
    log2(x+1)
  }
}


allphes<-phenames(cross)[grep("Pahal", phenames(cross))]
cross2<-transformPheno(cross, pheno.col=allphes, transf=transformExpression)
phe1<-pull.pheno(cross, pheno.col=goodphes)
phe2<-pull.pheno(cross2, pheno.col=goodphes)

cross3<-transformPheno(cross, pheno.col=allphes, transf=transformExpression2)
phe3<-pull.pheno(cross3, pheno.col=goodphes)


outliers<-sapply(phe1, function(x) getOutliers(as.numeric(x))$nOut)
outliers1<-sapply(phe2, function(x) getOutliers(as.numeric(x))$nOut)
outliers2<-sapply(phe3, function(x) getOutliers(as.numeric(x))$nOut)
out<-outliers[1,]+outliers[2,]
out1<-outliers1[1,]+outliers1[2,]
out2<-outliers2[1,]+outliers2[2,]
plot(jitter(out), jitter(out1), pch=".")
plot(jitter(out), jitter(out2), pch=".")

transformExpression
phe=normPhes[1:10]
nperm=10
ncluster=1
covar=trt
cross<-subset(cross, ind=!is.na(covar$trt))
covar=data.frame(trt=as.numeric(pull.pheno(cross, pheno.col="Treatment")))
set.seed(42)
system.time(addPerms<-scanone(cross=cross, pheno.col=phe, addcovar=covar, intcovar=NULL, perm.strata=covar[,1], n.perm=nperm, n.cluster=ncluster, verbose=T))
set.seed(42)
system.time(intPerms<-scanone(cross=cross, pheno.col=phe, addcovar=covar, intcovar=covar, perm.strata=covar[,1], n.perm=nperm, n.cluster=ncluster, verbose=T))

dat<-pull.pheno(cross, pheno.col=c("id","Treatment","Pahal.J02143"))
dat$Pahal.J02143<-transformExpression(dat$Pahal.J02143)
out<-outlierTest(lm(Pahal.J02143~Treatment, dat), cutoff=Inf, n.max=Inf)
df<-data.frame(pvals=out$p, id=names(out$p))
tp<-merge(phe2, df, by="id")
with(tp, plot(Pahal.J02143,-log10(pvals)))
packageVersion("qtl")
data(multitrait)
covar<-data.frame(test=sample(c(0,1), nind(multitrait), replace=T))
test<-scanone(multitrait, pheno.col=phenames(multitrait)[1:10], addcovar=covar, perm.strata=covar$test, n.perm=10)

covar2<-data.frame(test=sample(c(0,1, NA), nind(multitrait), replace=T))
test<-scanone(multitrait, pheno.col=phenames(multitrait)[1:10], addcovar=covar2, perm.strata=covar$test, n.perm=10)
test<-scanone(multitrait, pheno.col=phenames(multitrait)[1:10], addcovar=covar2, perm.strata=covar2$test, n.perm=10)

test<-full.pipe(cross=cross, phe=phenotype, nperm=10, verbose=T, ncluster=1, covar=trt)
save(test, file=paste("eqtlpipe3_",phenotype,".RData", sep=""))

test<-sapply(normPhes[20:35], function(x) paste("R CMD BATCH '--args phenotype=\"",x,"\"' toRunV3.R", sep=""))
write(test, file="testLauncher.txt")
