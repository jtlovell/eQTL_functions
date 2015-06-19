quantnorm.no0<-function(x) {
  y<-as.numeric(x)
  y[y==0]<-NA
  n=sum(!is.na(y),na.rm=T)
  y=rank(y)/(n+1)
  x=qnorm(y)
  x[is.infinite(x)]=NA
  x[x=="NaN"]<-NA
  x
}