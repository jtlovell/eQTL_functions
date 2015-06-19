binNormCat<-function(rawcounts, normCut=.2, binCut=.95, plotit=T){
  p0<-apply(rawcounts[,-which(colnames(rawcounts)=="id")], 2, function(x) sum(x==0)/length(x))
  mnon0<-apply(rawcounts[,-which(colnames(rawcounts)=="id")], 2, function(x) mean(x[x!=0], na.omit=T))
  binaryPhes<-colnames(rawcounts)[-which(colnames(rawcounts)=="id")][p0>.2 & p0<.95]
  normalPhes<-colnames(rawcounts)[-which(colnames(rawcounts)=="id")][p0<=.2]
  missingPhes<-colnames(rawcounts)[-which(colnames(rawcounts)=="id")][p0>=.95]
  hiPhes<-colnames(rawcounts)[-which(colnames(rawcounts)=="id")][mnon0>(mean(mnon0)*10)]
  
  length(normalPhes)
  length(binaryPhes)
  if(plotit){
    plot(c(0,1),c(min(log10(mnon0+1)),max(log10(mnon0+1))), type="n", bty="n",
         ylab="log10 mean counts (excluding 0's", xlab="proportion of 0 counts")
    points(p0[-which(colnames(rawcounts)=="id")][p0<=.2],log10(mnon0[-which(colnames(rawcounts)=="id")][p0<=.2]+1), col="black")
    points(p0[-which(colnames(rawcounts)=="id")][p0>.2],log10(mnon0[-which(colnames(rawcounts)=="id")][p0>.2]+1), col="red")
    points(p0[-which(colnames(rawcounts)=="id")][p0>.95],log10(mnon0[-which(colnames(rawcounts)=="id")][p0>.95]+1), col="black", pch=4)
    legend("topright", c(paste("n normal genes =", length(normalPhes)),
                         paste("n binary genes =", length(binaryPhes)),
                         paste("n missing genes =", length(missingPhes))),
           pch=c(1,1,4),col=c("black","red","black"), adj=0)
    title("distribution of raw counts for eQTL phenotypes") 
  }
  return(
    list(
      binaryPhes=binaryPhes,
      normalPhes=normalPhes,
      missingPhes=missingPhes,
      hiPhes=hiPhes))
}