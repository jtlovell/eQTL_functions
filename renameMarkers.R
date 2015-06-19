renameMarkers<-function(cross, crossfile, oldnames, newnames, outputName){
  cat("renaming markers and writing to a new file:", outputName, "\n")
  m<-markernames(cross)
  c<-read.csv(crossfile, header=T, quote="", na.strings="NA")
  n<-data.frame(oldnames,newnames)
  for(x in oldnames){
    if(x %in% m){
      colnames(c)[which(colnames(c)==x)]<-as.character(n$newnames[n$oldnames==x])
    }
  }
  if(is.na(c$id[1])) c$id[1:2]<-""
  write.csv(c,file=outputName, row.names=F, quote=F)
}