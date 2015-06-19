binary.perms<-function(cross, nperms, nclusters, series=NULL, trt){
  if(is.null(series)) series<-c(seq(from=.15, to=.25, by=.01),seq(from=.26, to=.35, by=.03),.4,.45,.5)
  binperms.out<-data.frame()
  for ( i in series){
    print(i)
    prop<-as.character(i)
    zeros<-rep(0,i*nind(cross))
    ones<-rep(1,nind(cross)-length(zeros))
    bin<-sample(c(zeros,ones), replace=F)
    bin<-bin+rnorm(1, mean=0, sd=0.0001)
    perm.bin<-data.frame(bin); colnames(perm.bin)<-"perm.bin"
    cross.binperm<-add.phenos(cross.null,perm.bin)
    cross.binperm<-calc.genoprob(cross.binperm, step=1,error.prob=0.001, map.function="kosambi")
    pens<-fourperms(cross=cross.binperm,phename="perm.bin", np=nperms, alpha=.05, nc=nclusters, covar=trt)[[1]]
    pens.out<-data.frame(rep(i,4),paste("pen",1:4,sep=""),pens); colnames(pens.out)<-c("prop.0","penalty.type","penalty")
    binperms.out<-rbind(binperms.out,pens.out)
  }
  par(mfrow=c(2,2))
  spline.out<-list()
  for (i in c("pen1", "pen2", "pen3", "pen4")){
    sub<-binperms.out[binperms.out$penalty.type==i,]
    modNew<-smooth.spline(sub$prop.0, sub$penalty, df=6)
    spline.out[[i]]<-modNew
  }
  return(list(
    penalties=binperms.out,
    splines=spline.out))
}
