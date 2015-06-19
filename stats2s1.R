stats2s1<-function(cross, phe, trt.covar, stats){
  s<-stats[stats$phenotype==phe,]
  c<-as.numeric(s[s$category=="cis","chr"])
  p<-as.numeric(s[s$category=="cis","pos"])

  g<-pull.genoprob(cross.norm, chr=stats[stats$phenotype==phe & stats$category=="cis","chr"],
                   include.pos.info=T, rotate=T)
  g<-g[g$chr==c,]; g<-g[which(abs(g$pos-p)==min(abs(g$pos-p))),-c(1:4)]
  cis<-as.numeric(apply(g,2,function(x) ifelse(x[1]>.9,1,ifelse(x[2]>.9,2,ifelse(x[3]>.9,3,NA)))))
  cis<-as.data.frame(cis)
  rownames(cis)<-colnames(g); colnames(cis)<-"cis"
  trt<-trt.covar
  rownames(trt)<-colnames(g)
  ac<-cbind(trt,cis)
  if(any(grepl("epi", s$category)) & any(grepl("trans.trt.int", s$category))){
    ic<-ac
  }else{
    if(any(grepl("epi", s$category))){
      ic<-cis
    }else{
      if(any(grepl("trans.trt.int", s$category))){
        ic<-trt
      }else{
        ic<-NULL
      }
    }
  }
  s1<-scanone(cross.norm, pheno.col=phe, method="hk", addcovar=ac, intcovar=ic)
  colnames(s1)[which(colnames(s1)=="lod")]<-phe
  return(s1)
}