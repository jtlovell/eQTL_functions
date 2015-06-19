fourperms<-function(cross, phename, alpha=0.05, np=100, nc=1, fit3=F, covar, pens4ctQTL=F){
  cis.trt<-data.frame(covar,
                      sample(c(0,1,2),nind(cross), replace=T))
  colnames(cis.trt)<-c("trt","cis")
  cis.trt3<-data.frame(covar,
                       sample(c(0,1,2),nind(cross), replace=T),
                       sample(c(0,1,2),nind(cross), replace=T))
  colnames(cis.trt3)<-c("trt","cis","trans")
  
  set.seed(42)
  cat("running perm set: \n 1..")
  s1perm.1<-scanone(cross, phe=phename, addcovar=cis.trt, method="hk",
                    n.perm=np, perm.strata=cis.trt$trt, n.cluster=nc,verbose=F)
  set.seed(42)
  cat("2..")
  s1perm.2<-scanone(cross, phe=phename, addcovar=cis.trt, intcovar=cis.trt$trt,  method="hk",
                    n.perm=np, perm.strata=cis.trt$trt, n.cluster=nc, verbose=F)
  set.seed(42)
  cat("3..")
  s1perm.3<-scanone(cross, phe=phename, addcovar=cis.trt, intcovar=cis.trt$cis,  method="hk",
                    n.perm=np, perm.strata=cis.trt$trt, n.cluster=nc, verbose=F)
  set.seed(42)
  cat("4..")
  s1perm.4<-scanone(cross, phe=phename, addcovar=cis.trt, intcovar=cis.trt,  method="hk",
                    n.perm=np, perm.strata=cis.trt$trt, n.cluster=nc, verbose=F)
  pens.1<-summary(s1perm.1, alpha=alpha)[1]
  pens.2<-summary(s1perm.2, alpha=alpha)[1]
  pens.3<-summary(s1perm.3, alpha=alpha)[1]
  pens.4<-summary(s1perm.4, alpha=alpha)[1]
  if(fit3){
    set.seed(42)
    cat("5..")
    s1perm.5<-scanone(cross, phe=phename, addcovar=cis.trt3, method="hk",
                      n.perm=np, perm.strata=cis.trt3$trt, n.cluster=nc, verbose=F)
    set.seed(42)
    cat("6..")
    s1perm.6<-scanone(cross, phe=phename, addcovar=cis.trt3, intcovar=cis.trt3$trt,  method="hk",
                      n.perm=np, perm.strata=cis.trt3$trt, n.cluster=nc, verbose=F)
    set.seed(42)
    cat("7..")
    s1perm.7<-scanone(cross, phe=phename, addcovar=cis.trt3, intcovar=cis.trt3$cis,  method="hk",
                      n.perm=np, perm.strata=cis.trt3$trt, n.cluster=nc, verbose=F)
    set.seed(42)
    cat("8..\n")
    s1perm.8<-scanone(cross, phe=phename, addcovar=cis.trt3, intcovar=cis.trt3[,c("cis","trt")],  method="hk",
                      n.perm=np, perm.strata=cis.trt3$trt, n.cluster=nc, verbose=F)
    pens.5<-summary(s1perm.5, alpha=alpha)[1]
    pens.6<-summary(s1perm.6, alpha=alpha)[1]
    pens.7<-summary(s1perm.7, alpha=alpha)[1]
    pens.8<-summary(s1perm.8, alpha=alpha)[1]
    if(pens4ctQTL){
      
    }
    return(c(pens.1,pens.2,pens.3,pens.4,pens.5,pens.6,pens.7,pens.8))
  }else{
    return(list(
      penalties = c(pens.1,pens.2,pens.3,pens.4),
      permutations= list(s1perm.1,s1perm.2,s1perm.3,s1perm.4))
    )
  }  
}