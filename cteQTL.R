cteQTL<-function(cross, chromosome, position, phe, pens=NULL, forms.in=NULL,trt, wiggle=1, refine.qtl=FALSE, ci.method="bayes", estCI=TRUE){
  #for a standard cis-trans eqtl, we search exhaustively through mode space
  #   the cis eQTL ("Q1") and treatment covariate ("trt") are always included. 
  if(is.null(forms.in)){
    forms<-c("y ~ Q1 + trt",
             "y ~ Q1 + Q1*trt + trt",
             "y ~ Q1 + Q2 + trt",
             "y ~ Q1 + Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + Q2*trt + trt")
  }else{
    forms<-forms.in
  }
  #get name for phe, so that we index by character not number
  if(is.factor(phe)){
    phe<-as.character(phe)
  }
  if(is.numeric(phe)){
    phe<-phenames(cross)[phe]
  }
  mods<-list()
  fits<-list()
  scans<-list()
  cis.pos<-position
  #if our position is not 100%, allow the position to move
  if(wiggle>0 | length(wiggle>1)){
    s1.1a<-scanone(cross, pheno.col=phe, method="hk",  addcovar=trt)
    scans[[1]]<-s1.1a
    s1.2a<-scanone(cross, pheno.col=phe, method="hk",  addcovar=trt, intcovar=trt)
    scans[[2]]<-s1.2a
    s1.1<-as.data.frame(s1.1a)
    s1.2<-as.data.frame(s1.2a)
    if(length(wiggle>1)){
      s1.wiggle<-cbind(s1.1[s1.1$chr==chromosome & s1.1$pos<max(wiggle) & s1.1$pos>min(wiggle),],
                       s1.2$lod[s1.1$chr==chromosome & s1.1$pos<max(wiggle) & s1.1$pos>min(wiggle)])
    }else{
      s1.wiggle<-cbind(s1.1[s1.1$chr==chromosome & s1.1$pos<position+wiggle & s1.1$pos>position-wiggle,],
                       s1.2$lod[s1.1$chr==chromosome & s1.1$pos<position+wiggle & s1.1$pos>position-wiggle])
    }
    
    wig.sum<-s1.wiggle[,3]+s1.wiggle[,4]
    position.new<-s1.wiggle$pos[which(wig.sum==max(wig.sum))[1]]
    wiggle.move<-abs(position.new-position)
    position<-position.new
  }
  # the first two models need to be built manually because of a single QTL
  trt<-data.frame(trt=trt)
  cis.eqtl<-makeqtl(cross=cross, chr=chromosome, pos=position, what="prob")
  mods[[1]]<-cis.eqtl
  mods[[2]]<-cis.eqtl
  fit.1<-fitqtl(cross, qtl=cis.eqtl, 
                formula=forms[1], pheno.col=phe, 
                covar=trt, method="hk", dropone=T)
  
  fit.2<-fitqtl(cross, qtl=cis.eqtl, 
                formula=forms[2], pheno.col=phe, 
                covar=trt, method="hk", dropone=T)
  
  lod.1<-data.frame(fit.1$result.drop)
  fits[[1]]<-lod.1
  lod.nullmod<-lod.1[cis.eqtl$name,"LOD"]
  lod.1<-sum(lod.1[which(rownames(lod.1)!=colnames(trt)),"LOD"])-lod.nullmod
  lod.2<-data.frame(fit.2$result.drop)
  fits[[2]]<-lod.2
  lod.2<-sum(lod.2[which(rownames(lod.2)!=colnames(trt)),"LOD"])-lod.nullmod
  lod.all<-c(lod.1,lod.2)
  # fit the remaining 7 models
  for (i in 3:length(forms)){
    form.in<-forms[i]
    scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                   method="hk", covar=trt, pheno.col=phe)
    scans[[i]]<-scan
    mod <- addtoqtl(cross, cis.eqtl, max(scan)$chr, max(scan)$pos)
    if(length(unique(mod$chr))==1 & abs(diff(mod$pos)) < 35){
      scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                     chr = chrnames(cross)[-which(chrnames(cross)==cis.eqtl$chr)],
                     method="hk", covar=trt, pheno.col=phe)
      scans[[i]]<-scan
      mod <- addtoqtl(cross, cis.eqtl, max(scan)$chr, max(scan)$pos)
    }
    fit <- fitqtl(cross, qtl=mod, 
                  formula=form.in, pheno.col=phe, 
                  covar=trt, method="hk", dropone=T, get.ests=T)
    lod.out<-data.frame(fit$result.drop)
    fits[[i]]<-lod.out
    lodi<-sum(lod.out[which(rownames(lod.out)!=colnames(trt)),"LOD"])-lod.nullmod
    lod.all<-c(lod.all,lodi)
    mods[[i]]<-mod
  }
  #deterime which mode fit is best by taking the best pLOD score from the QTL and QTL*trt interactions
  plods<-lod.all-pens
  best<-which(plods == max(plods))
  best.form<-forms[best]
  best.lod<-as.numeric(lod.all[best])
  best.mod<-mods[[best]]
  best.fit<-fits[[best]]
  best.scan<-scans[[best]]
  cis.name<-cis.eqtl$name
  #refine if desired
  if(refine.qtl){
    best.mod<-refineqtl(cross, pheno.col=phe, qtl=best.mod, covar=trt, method="hk", model="normal",
                        formula=best.form, keeplodprofile=T, verbose=F)
    fitref <- fitqtl(cross, qtl=best.mod, 
                     formula=best.form, pheno.col=phe, 
                     covar=trt, method="hk", dropone=T, get.ests=T)
    best.fit<-data.frame(fitref$result.drop)
    cis.name<-best.mod$name[1]
    
    if(estCI){
      cis<-data.frame()
      nqtls<-nqtl(best.mod)
      for (j in 1:nqtls){
        if(ci.method=="drop"){
          ciout<-lodint(best.mod,qtl.index=j, expandtomarkers=F)
        }else{
          ciout<-bayesint(best.mod,qtl.index=j, expandtomarkers=F)
        }
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      cis$term.id<-best.mod$name
    }
  }
  
  
  #make the output object
  
  all.out<-data.frame(rownames(best.fit)); colnames(all.out)[1]<-"term.id"
  all.out$cisPostion<-cis.pos
  
  mgsub <- function(pattern, replacement, x, ...) {
    result <- x
    for (i in 1:length(pattern)) {
      result <- gsub(pattern[i], replacement[i], result, ...)
    }
    result
  }
  
  cis.ints<-c(paste(":",cis.name, sep=""),paste(cis.name,":", sep=""))
  qnames<-rownames(best.fit)
  qnames<-mgsub(c(":trt","trt:",cis.ints),rep("",4),qnames)
  chr.out<-sapply(qnames, function (x) strsplit(x,"@")[[1]][1]) 
  chr.out[chr.out=="trt"]<-NA
  pos.out<-sapply(qnames, function (x) strsplit(x,"@")[[1]][2]) 
  pos.out[pos.out=="trt"]<-NA
  
  all.out$chr<-chr.out
  all.out$pos<-pos.out
  
  #add category information
  if(grepl("Q2", best.form)){
    trans.name<-best.mod$name[!best.mod$name==cis.name]
    category<-rownames(best.fit)
    category[intersect(grep("trt", category), grep(paste(trans.name,collapse="|"), category))]<-"trans.trt.int"
    category[intersect(grep("trt", category), grep(":", category))]<-"cis.trt.int"
    category[grep(":", category)]<-"epi"
    category[grep(paste(trans.name,collapse="|"), category)]<-"trans"
    category[category==cis.name]<-"cis"
    
  }else{
    category<-rownames(best.fit)
    category[intersect(grep("trt", category), grep(":", category))]<-"cis.trt.int"
    category[grep("@", category)]<-"cis"
  }
  all.out$category<-category
  
  #add phenotype name
  all.out$phenotype<-phe
  
  #add estimates
  fit<-fitqtl(cross, qtl=best.mod, 
              formula=best.form, pheno.col=phe, 
              covar=trt, method="hk", dropone=T, get.ests=T)
  ests.all<-summary(fit)$ests
  rows.dom<-intersect(grep("d",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE))
  rows.add<-intersect(intersect(grep("a",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE)), 
                      grep(colnames(trt),rownames(ests.all),invert=TRUE))
  rows.cov<-which(rownames(ests.all) == "trt")
  covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
  colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
  rownames(covar.ests_out)<-"trt"
  if(grepl("Q2", best.form)){
    dom.ests_out<-data.frame(ests.all[rows.dom,1:3]);colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
    add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
    qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
  }else{
    dom.ests_out<-data.frame(t(data.frame(ests.all[rows.dom,1:3])));colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
    add.ests_out<-data.frame(t(data.frame(ests.all[rows.add,1:3]))); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
    qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
    row.names(qtl.ests_out)<-rownames(ests.all)[rows.add]
    qtl.ests_out<-data.frame(qtl.ests_out)
    rownames(covar.ests_out)<-"trt"
  }
  
  ests.out<-rbind(covar.ests_out,qtl.ests_out)
  ests.out$term.id<-gsub("a","",rownames(ests.out))
  
  all.out<-merge(all.out, ests.out, by="term.id", all.x=T)
  if(estCI){
    all.out<-merge(all.out, cis, by="term.id", all.x=T)
  }
  
  lod.df<-data.frame(best.fit)
  lod.df$term.id<-rownames(lod.df)
  all.out<-merge(all.out, lod.df, by="term.id", all.x=T)
  
  all.out
  
}