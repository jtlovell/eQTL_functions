find.ciseqtl<-function(cross,phe,pos.gene,chr.gene, model.type,threshold){
  r<-refineqtl(cross, 
               pheno.col=phe,
               chr=chr.gene, 
               pos=pos.gene,
               model=model.type, 
               method="hk",
               keeplodprofile=T, 
               verbose=F)
  b<-bayesint(r,qtl.index=1, 
              prob=.95, 
              expandtomarkers=F) 
  cq<-ifelse(pos.gene > b$pos[1] & b$pos[2] > pos.gene, 
             "yes",
             "no")
  q<-ifelse(max(b$lod>threshold),"yes","no")
  out<-data.frame(phe,chr.gene, pos.gene,b$pos[1],b$pos[2], q ,cq)
  colnames(out)<-c("phenotype","chr","pos.gene","pos.lo.qtl","pos.hi.qtl","is.qtl","is.cis")
  return(out)
}