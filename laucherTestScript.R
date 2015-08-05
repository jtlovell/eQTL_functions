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

test<-full.pipe(cross=cross, phe=phenotype, nperm=10, verbose=T, ncluster=1, covar=trt)
save(test, file=paste("eqtlpipe3_",phenotype,".RData", sep=""))


sbatch normal.slurm test1.launcher

/tmp/slurmd/job5578090/slurm_script: line 20: /paramrun: No such file or directory

$T