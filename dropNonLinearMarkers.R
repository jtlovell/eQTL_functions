dropNonLinearMarkers<-function(markers, bp, cm, plotit=T,...){
  dat<-data.frame(markers, bp, cm, stringsAsFactors=F)
  dat<-dat[order(dat$bp),]
  cm=dat$cm
  bp=dat$bp
  if(plotit){
    plot(bp,cm, type="n",...)
    points(bp[diff(bp)/diff(cm)>=0],cm[diff(bp)/diff(cm)>=0])
    points(bp[diff(bp)/diff(cm)<0],cm[diff(bp)/diff(cm)<0], col="red", pch=19)
  }
  if(min(diff(bp)/diff(cm))<0){
    while(min(diff(bp)/diff(cm))<0){
      bad<-which(diff(cm)==min(diff(cm)))
      dat<-dat[-c(bad, bad+1),]
      cm=dat$cm
      bp=dat$bp
      if(plotit){
        plot(bp,cm, type="n",...)
        points(bp[diff(bp)/diff(cm)>=0],cm[diff(bp)/diff(cm)>=0])
        points(bp[diff(bp)/diff(cm)<0],cm[diff(bp)/diff(cm)<0], col="red", pch=19)
      }
    }
  }
  cat("dropping markers: ",markers[!markers %in% dat$markers], "\n")
  return(dat$markers)
}