linearBp2cM<-function(gff, markers, cm, start, plotit=T,...){
  d<-data.frame(markers, cm, start)
  d<-d[order(d$start),]
  if(plotit)  plot(1,1, ylim=range(d$cm), xlim=range(d$start), type="n",...)
  out<-list()
  for(i in 2:length(d$start)){
    
    p<-gff[gff$start>=d$start[i-1] & gff$start<d$start[i],]
    
    d1<-d[(i-1):i,]
    if(plotit) lines(d1[,c("start","cm")])
    mod<-lm(cm~start, data=d1)
    pred<-as.numeric(predict(mod, newdata=p))
    out[[i]]<-data.frame(predictedCM=pred, geneID=p$geneID)
  }
  return(ldply(out, data.frame))
}