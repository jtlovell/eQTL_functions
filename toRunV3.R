#what is needed to run V3 pipe:
rm(list=ls())
load("eqtl2015_finalbinnormCrosses.RData")

pkg <- c("rrBLUP","doMC","doSNOW","cluster","reshape2","reshape","parallel","snow","RCurl","plyr","qtl","car","MASS","ggplot2")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))

function.names<-c("pipeV3.R")
us<-paste("https://raw.githubusercontent.com/jtlovell/eQTL_functions/master/",function.names, sep="")
for(i in 1:length(function.names)){
  cat("loading",function.names[i],"\n")
  script <- getURL(us[i], ssl.verifypeer = FALSE)
  eval(parse(text = script))
}

test<-full.pipe(cross=cross, phe=normPhes[30], nperm=10, verbose=T, ncluster=1, covar=trt)
