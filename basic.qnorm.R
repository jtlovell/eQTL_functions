quantnorm<-function(x) {
  n=sum(!is.na(x),na.rm=T)
  x=rank(x)/(n+1)
  x=qnorm(x)
  x[is.infinite(x)]=NA
  x
}