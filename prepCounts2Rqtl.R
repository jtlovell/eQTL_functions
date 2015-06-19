prepCounts2Rqtl<-function(dat){
  dat<-dat[,-which(colnames(dat)=="FH058_H311.1_591")]  #drop 58, #2.
  dat<-dat[,-which(colnames(dat) %in% c("v11.id","chr","type","start", "end", "strand"))]
  rownames(dat)<-dat$id
  dat$id<-NULL
  dat<-data.frame(t(dat))
  dat$id<-as.numeric(as.character(sapply(row.names(dat),function(x) substr(x, 3,5))))
  return(dat)
}