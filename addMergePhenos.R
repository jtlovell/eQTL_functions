addMergePhenos<-function(cross, dat){
  id<-as.numeric(as.character(getid(cross)))
  dat$id<-as.numeric(as.character(dat$id))
  dat<-dat[dat$id %in% id,]
  nodat<-id[!id %in% dat$id]
  nadf<-dat[1:length(nodat),]
  nadf[, -which(colnames(nadf)=="id")]<-NA
  nadf$id<-nodat
  dat<-rbind(dat,nadf)
  dat$id<-as.numeric(dat$id)
  dat<-dat[order(dat$id),]
  cross$pheno <- cbind(cross$pheno, dat)
  return(cross)
}