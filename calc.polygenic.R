calc.polygenic<-function(cross, chrs, poss, phe, trt, best.form){
  g<-pull.genoprob(cross, include.pos.info=T, rotate=T, omit.first.prob=F)
  mars<-find.marker(cross,chrs, poss)
  to.vc<-data.frame(cbind(pull.pheno(cross,pheno.col=phe),trt))
  colnames(to.vc)[1]<-c("y")
  for(i in mars){
    gs<-g[g$marker ==i,]
    gs.out<-as.data.frame(as.numeric(unlist(apply(gs[,-c(1:4)],2,function(x) which(x==max(x))[1]))))
    colnames(gs.out)<-i
    to.vc<-cbind(to.vc,gs.out)
  }
  colnames(to.vc)[3:length(colnames(to.vc))]<-paste("Q",seq(from=1, to=length(colnames(to.vc))-2), sep="")
  snps<-as.data.frame(pull.geno(cross))
  to.vc$SNP<-snps
  to.vc<<-to.vc
  v1<-varComp(as.formula(best.form) , to.vc, ~ibs(SNP))
  pg<-as.numeric(v1$varComps)
  err<-as.numeric(v1$sigma2)
  prop.pg<-pg/(pg+err)
  return(c(pg,err,prop.pg))
}