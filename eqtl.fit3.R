eqtl.fit3<-function(cross, ct.obj, trt, pens=NULL, q3forms=NULL, refine.qtl=FALSE){
  if(is.null(pens)){pens<-pens.all<-c(0.000000,  1.301030,  3.967144,  5.268174,  6.144188,  6.922426,  8.223456,  9.576449, 10.877479)}
  if(is.null(q3forms)){q3forms<-c("+ Q3", "+ Q3 + Q3*trt", "+ Q3 + Q1*Q3", "+ Q3 + Q3*trt + Q3*Q1") }
  
  chr<-ct.obj$chr
  pos<-ct.obj$pos
  best.form<-ct.obj$formula
  plods<-ct.obj$plods
  plods.init<-plods
  phe<-ct.obj$stats$phenotype[1]
  lod.nullmod<-ct.obj$ciseqtl.lod
  cis.eqtl<-as.character(ct.obj$stats[ct.obj$stats$category=="cis","term.id"])
  
  best.mod<-makeqtl(cross, chr=chr, pos=pos, what="prob")
  lod.all<-vector()
  mods<-list()
  fits<-list()
  forms.out<-vector()
  for (i in 1:length(q3forms)){
    form.in3<-paste(best.form,q3forms[i])
    forms.out<-c(forms.out, form.in3)
    scan3 <- addqtl(cross, qtl=best.mod, formula=form.in3, 
                    method="hk", covar=trt, pheno.col=phe)
    mod3 <- addtoqtl(cross, best.mod, max(scan3)$chr, max(scan3)$pos)
    if(length(unique(mod3$chr))==2){
      chr.dup<-mod3$chr[(duplicated(mod3$chr))]
      diff.pos<-abs(diff(mod3$pos[mod3$chr==chr.dup]))
      if(diff.pos < 50 & !is.na(diff.pos)){
        scan3 <- addqtl(cross, qtl=best.mod, formula=form.in3,
                        chr = chrnames(cross)[-which(chrnames(cross) %in% chr)],
                        method="hk", covar=trt, pheno.col=phe)
        mod3 <- addtoqtl(cross, best.mod, max(scan3)$chr, max(scan3)$pos)
      }
    }else{
      if(length(unique(mod3$chr))==1){
        scan3 <- addqtl(cross, qtl=best.mod, formula=form.in3,
                        chr = chrnames(cross)[-which(chrnames(cross) == unique(mod3$chr))],
                        method="hk", covar=trt, pheno.col=phe)
        mod3 <- addtoqtl(cross, best.mod, max(scan3)$chr, max(scan3)$pos)
      }
    }
    fit3 <- fitqtl(cross, qtl=mod3, 
                   formula=form.in3, pheno.col=phe, 
                   covar=trt, method="hk", dropone=T, get.ests=T)
    lod.out3<-data.frame(fit3$result.drop)
    lodi<-sum(lod.out3[which(rownames(lod.out3)!=colnames(trt)),"LOD"])-lod.nullmod
    lod.all<-c(lod.all,lodi)
    mods[[i]]<-mod3
    fits[[i]]<-data.frame(fit3$result.drop)
  }
  best<-which(plods==max(plods))
  max.2qtl.plod<-max(plods)  
  
  pens3<-c((pens[best] + pens[3]),  
           (pens[best] + pens[4]),
           (pens[best] + pens[6]),
           (pens[best] + pens[7]))
  plods3<-lod.all-pens3
  best3<-which(plods3==max(plods3))
  
  #output
  if(max(plods3)>max(plods)){
    plods<-plods3
    best<-best3
    best.form<-forms.out[best]
    best.lod<-as.numeric(lod.all[best])
    best.mod<-mods[[best]]
    best.fit<-fits[[best]]
    cis.name<-cis.eqtl
    if(refine.qtl){
      best.mod<-refineqtl(cross, pheno.col=phe, qtl=best.mod, covar=trt, method="hk", model="normal",
                          formula=best.form, keeplodprofile=F, verbose=F)
      fitref <- fitqtl(cross, qtl=best.mod, 
                       formula=best.form, pheno.col=phe, 
                       covar=trt, method="hk", dropone=T, get.ests=T)
      best.fit<-data.frame(fitref$result.drop)
      cis.name<-best.mod$name[1]
    }
    #make the output object
    all.out<-data.frame(rownames(best.fit)); colnames(all.out)[1]<-"term.id"
    mgsub <- function(pattern, replacement, x, ...) {
      if (length(pattern)!=length(replacement)) {
        stop("pattern and replacement do not have the same length.")
      }
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
    trans.name<-best.mod$name[!best.mod$name==cis.name]
    category<-rownames(best.fit)
    category[intersect(grep("trt", category), grep(paste(trans.name,collapse="|"), category))]<-"trans.trt.int"
    category[intersect(grep("trt", category), grep(":", category))]<-"cis.trt.int"
    category[grep(":", category)]<-"epi"
    category[grep(paste(trans.name,collapse="|"), category)]<-"trans"
    category[category==cis.name]<-"cis"
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
                        grep(colnames(covar),rownames(ests.all),invert=TRUE))
    rows.cov<-which(rownames(ests.all) == "trt")
    covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
    colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
    rownames(covar.ests_out)<-"trt"
    dom.ests_out<-data.frame(ests.all[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
    add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
    qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
    ests.out<-rbind(covar.ests_out,qtl.ests_out)
    ests.out$term.id<-gsub("a","",rownames(ests.out))
    
    all.out<-merge(all.out, ests.out, by="term.id", all.x=T)
    
    lod.df<-data.frame(best.fit)
    lod.df$term.id<-rownames(lod.df)
    all.out<-merge(all.out, lod.df, by="term.id", all.x=T)
    plods<-data.frame(t(data.frame(plods3)))
    colnames(plods)<-c(paste("lod",10:13, sep="."))
    plods.out<-cbind(plods.init, plods)
    return(list(formula=best.form,
                plods=plods.out,
                ciseqtl.lod=lod.nullmod,
                chr=best.mod$chr,
                pos=best.mod$pos,
                stats=all.out,
                cis.position.move=NA))
  }else{
    return(ct.obj)
  }
}