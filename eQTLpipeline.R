scanFormula<-function(cross, chromosome, position, phe, pens=NULL, forms.in=NULL,trt, wiggle=1, verbose=F, saveFULL=F){
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
    if(length(wiggle)>1){
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
  
  fits[[1]]<-data.frame(fit.1$result.drop)
  fits[[2]]<-data.frame(fit.2$result.drop)
  lod.1<-fit.1$result.full["Model","LOD"]
  lod.2<-fit.2$result.full["Model","LOD"]
  if(verbose){
    cat(forms[1], "...  max scan LOD / model LOD :", max(s1.1$lod), "/",fit.1$result.full["Model","LOD"], "\n")
    cat(forms[2], "...  max scan LOD / model LOD :", max(s1.2$lod), "/",fit.2$result.full["Model","LOD"], "\n")
  }
  lod.all<-c(lod.1,lod.2)
  transLowCI<-c(NA,NA)
  transHighCI<-c(NA,NA)
  # fit the remaining 7 models
  for (i in 3:length(forms)){
    form.in<-forms[i]
    scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, model="normal",
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
                  covar=trt, method="hk", dropone=T, get.ests=F)
    fits[[i]]<-data.frame(fit$result.drop)
    lodi<-fit$result.full["Model","LOD"]
    lod.all<-c(lod.all,lodi)
    mods[[i]]<-mod
    
    ciout<-lodint(scan, chr=mod$chr[2])
    transLowCI<-c(transLowCI, ciout[1,2])
    transHighCI<-c(transHighCI,ciout[3,2])
    
    if(verbose){
      cat(forms[i], "...  max scan LOD / model LOD :", max(scan$lod), "/",fit$result.full["Model","LOD"], "\n")
    }
  }
  if(saveFULL){
    return(list(models=mod, fitqtlResults=fit, scansData=scans, modelLOD=lod.all))
  }else{
    out<-sapply(mods, function(x) {
      if(length(x$chr)==1){
        chr<-as.numeric(as.character(c(x$chr,NA)))
        pos<-c(x$pos,NA)
      } else{
        chr<-as.numeric(as.character(x$chr))
        pos<-x$pos
      }
      c(chr[1], chr[2],pos[1], pos[2])
    })
    dat<-data.frame(t(out))
    colnames(dat)<-c("chrCis","chrTrans", "posCis", "posTrans")
    dat$modelLOD<-lod.all
    out<-data.frame(phenotype=phe, penalties=pens, formula = forms, dat, pLOD=dat$modelLOD - pens, transLowCI,transHighCI)
    return(out)
  }
}

bestpLOD<-function(dat){
  chr=dat$chrCis[dat$pLOD==max(dat$pLOD)]
  pos<-dat$posCis[dat$pLOD==max(dat$pLOD)]
  if(!is.na(dat$chrTrans[dat$pLOD==max(dat$pLOD)])){
    chr<-c(chr, dat$chrTrans[dat$pLOD==max(dat$pLOD)])
    pos<-c(pos, dat$posTrans[dat$pLOD==max(dat$pLOD)])
  }
  formula<-as.character(dat$formula[dat$pLOD==max(dat$pLOD)])
  list(chr=chr, pos=pos, formula=formula)
}

makeModel<-function(cross,dat){
  chr=dat[["chr"]]
  pos=dat[["pos"]]
  if(length(chr)==2){
    mod<-makeqtl(cross, chr=chr, pos=pos, qtl.name=c("cis","trans"), what = "prob")
  }else{
    mod<-makeqtl(cross, chr=chr, pos=pos, qtl.name="cis", what = "prob")
  }
  mod
}

cteQTLStats<-function(cross, model, formula, phe, covar=NULL, scanFormulaOutput=NULL){
  err<-"good"
  tryCatch(fit<-fitqtl(cross,
                       pheno.col=phe,
                       qtl=model,
                       formula=formula,
                       get.ests=T,dropone=T,covar=covar,
                       method="hk"),
           error=function(e) { 
             fit<-fitqtl(cross,
                                  pheno.col=phe,
                                  qtl=model,
                                  formula=formula,
                                  get.ests=F,dropone=T,covar=covar,
                                  method="hk")
           },
           finally = function(e) { 
             fit<-fitqtl(cross,
                         pheno.col=phe,
                         qtl=model,
                         formula=formula,
                         get.ests=F,dropone=T,covar=covar,
                         method="hk")
           }
           ) 
  nqtls<-nqtl(model)
  nterms<-sum(countqtlterms(formula, ignore.covar=F)[c(1,4)])
  ncovar<-length(covar)
  d1<-data.frame(fit$result.drop)
  d1$term<-row.names(d1)
  if(!is.null(summary(fit)$ests)){
    #parse estimated effects data
    ests<-summary(fit)$ests
    ests<-ests[!(grepl("cisa:trans",rownames(ests)) | grepl("cisd:trans",rownames(ests))),]
    estsA<-ests[c(grep("a$", rownames(ests)), grep("a:trt", rownames(ests))),]
    estsD<-ests[c(grep("d$", rownames(ests)), grep("d:trt", rownames(ests))),]
    estsT<-ests[c("Intercept","trt"),]
    if(length(estsA)==3){
      estsA<-data.frame(t(data.frame(estsA)))
      estsD<-data.frame(t(data.frame(estsD)))
      rownames(estsA)<-"cis"
      rownames(estsD)<-"cis"
    }else{
      estsA<-data.frame(estsA)
      estsD<-data.frame(estsD)
      rownames(estsA)<-gsub("cisa", "cis", rownames(estsA))
      rownames(estsA)<-gsub("transa", "trans", rownames(estsA))
      rownames(estsD)<-gsub("cisd", "cis", rownames(estsD))
      rownames(estsD)<-gsub("transd", "trans", rownames(estsD))
    }
    
    estsT<-data.frame(estsT)
    colnames(estsA)<-paste(colnames(estsA), "add",sep="_")
    colnames(estsD)<-paste(colnames(estsD), "dom",sep="_")
    colnames(estsT)<-paste(colnames(estsT), "add",sep="_")
    estsA$term<-as.character(rownames(estsA))
    estsD$term<-as.character(rownames(estsD))
    estsT$term<-c("Intercept","trt")
    
    estsAD<-merge(estsA, estsD, by="term")
    
    estsT$est_dom<-NA; estsT$SE_dom<-NA; estsT$t_dom<-NA
    ests.out<-rbind(estsT, estsAD)
    out<-merge(d1,ests.out, by="term", all=T)
    out$phenotype<-phe
    out$formula=formula
    out$modelLOD<-NA; out$pLOD<-NA; out$transLowCI<-NA; out$transHighCI<-NA
  }else{
    out<-d1
    for(i in c("est_add", "SE_add", "t_add", "est_dom", "SE_dom", "t_dom")) out[,i]<-NA
    out$modelLOD<-NA; out$pLOD<-NA; out$transLowCI<-NA; out$transHighCI<-NA
  }
  if(!is.null(scanFormulaOutput)){
    dat<-scanFormulaOutput[scanFormulaOutput$formula==formula,]
    out$transLowCI[grep("trans", out$term)]<-dat$transLowCI
    out$transHighCI[grep("trans", out$term)]<-dat$transHighCI
    out$cisChr<-dat$chrCis
    out$cisPos<-dat$posCis
    out$transChr<-dat$chrTrans
    out$transPos<-dat$posTrans
    out$modelLOD<-dat$modelLOD
    out$pLOD<-dat$pLOD
  }
  out[,c(which(colnames(out) %in% c("phenotype","formula","modelLOD","pLOD")),
         which(!colnames(out) %in% c("phenotype","formula","modelLOD","pLOD")))]
}

add1<-function(cross, formula, covar, model, phe){
  chr<-model$chr
  pos<-model$pos
  scan<-addqtl(cross=cross, qtl=model, covar=trt, formula=form, method="hk", model="normal", pheno.col=phe)
  max(scan)
}

eQTLcistransPipe<-function(cross, gene, gff, trt, geno.probs, snps, ...){
  pos=gff$predictedCM[gff$geneID==gene]
  chr=gff$lg[gff$geneID==gene]
  
  if(verbose) cat("1. exhaustively searching for the best cis/trans eQTL model\n")
  scanbypLOD<-scanFormula(cross=cross, chromosome=chr, position=pos, phe=gene, trt=trt, wiggle=5, pens=pens)
  
  if(verbose) cat("2. compiling statistics\n")
  toppLOD<-bestpLOD(scanbypLOD)
  models<-makeModel(cross=cross, dat=toppLOD)
  result<-cteQTLStats(cross=cross, model=models, formula=toppLOD[["formula"]], phe=gene,covar=trt,
                      scanFormulaOutput=scanbypLOD)
  
  if(verbose) cat("3. calculating the extent of polygenic inheritance\n")
  chr=as.numeric(models[["chr"]])
  pos=as.numeric(models[["pos"]])
  pg<-calc.polygenic(cross=cross, geno.probs=geno.probs, snps=snps,
                     chrs=chr, 
                     poss=pos, 
                     best.form=models[["formula"]], 
                     trt=trt, phe=gene)
  
  if(grep("Q2", toppLOD[["formula"]])){
    if(verbose) cat("4. searching for a 2nd trans eQTL\n")
    trans2<-add1(cross=cross, formula=models[["formula"]], phe=gene)
  }else{
    trans2<-NULL
  }
  return(list(scanbypLOD=scanbypLOD, stats=result, polygenic=pg, trans2=trans2))
}
