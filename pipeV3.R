# eQTL pipeline V3.1

genoprob2marker<-function(cross, genotypes=c("AA","AB","BB"), threshold=0.9){
  dat<-pull.genoprob(cross,omit.first.prob=F, rotate=T,include.pos.info=T)
  out<-data.frame(apply(dat[,5:ncol(dat)],2,function(x) {
    tapply(x, dat$marker, function(y) {
      ifelse(max(y)<threshold,
             NA,
             genotypes[which(y==max(y))]
      )
    })
  }))
  out$marker<-as.character(row.names(out))
  map<-dat[!duplicated(dat$marker),c(1,3,4)]
  merge(map,out)
}

calcPolygenic<-function(data, pheno=phenoMatrix, kinshipMatrix, covar="trt"){
  v3<-kin.blup(data=data, geno="id", pheno=pheno, covariate=covar, K=kinshipMatrix)
  c(v3$Vg, v3$Ve)
}

compileScanone<-function(cross, s1Add, s1Int, addPerms, intPerms, phe){
  summaryAdd<-summary(s1Add,  perms=addPerms, pvalues=TRUE)
  maxAdd<-summaryAdd[summaryAdd$lod==max(summaryAdd$lod),]
  summaryInt<-summary(s1Int,  perms=intPerms, pvalues=TRUE)
  maxInt<-summaryInt[summaryInt$lod==max(summaryInt$lod),]
  diffPerms<-intPerms-addPerms
  s1Diff<-s1Int-s1Add
  summaryDiff<-summary(s1Diff,  perms=diffPerms, pvalues=TRUE)
  maxDiff<-summaryDiff[summaryDiff$lod==max(summaryDiff$lod),]
  if(maxDiff$pval < 0.05){
    bestForm<-"y ~ Q1 + Q1*trt + trt"
    bestMod<-makeqtl(cross, chr=maxInt$chr, pos=maxInt$pos, what="prob")
  }else{
    bestForm<-"y ~ Q1 + trt"
    bestMod<-makeqtl(cross, chr=maxAdd$chr, pos=maxAdd$pos, what="prob")
  }    
  colnames(maxAdd)<-paste("scan_", colnames(maxAdd), sep="")
  colnames(maxInt)<-paste("scan_", colnames(maxInt), sep="")
  colnames(maxDiff)<-paste("diff_", colnames(maxDiff), sep="")
  maxDiffNA<-maxDiff; maxDiffNA[1,]<-NA
  out<-rbind(data.frame(maxAdd, maxDiffNA), data.frame(maxInt, maxDiff))
  out$formula=c("y ~ Q1 + trt", "y ~ Q1 + Q1*trt + trt")
  out$phe=phe
  return(list(maxScanData=out, bestForm=bestForm,bestMod=bestMod))
}

add1QTL<-function(cross, intial.model, new.formula, perm2intepi, perm2int, perm2epi, perm2, phe, ncluster){
  baseScan<-addqtl(cross, qtl=intial.model, formula=new.formula[nchar(new.formula)==min(nchar(new.formula))],
                   model="normal", method="hk", covar=trt, pheno.col=phe)
  out<-mclapply(new.formula,  mc.preschedule = F, mc.cores=ncluster, function(x) {
    if(grepl("Q2[*]trt", x) & grepl("Q1[*]Q2", x)){
      perms<-perm2intepi
      is.base=FALSE
    }else{
      if(grepl("Q2[*]trt", x)){
        perms<-perm2int
        is.base=FALSE
      }else{
        if(grepl("Q1[*]Q2", x)){
          perms<-perm2epi
          is.base=FALSE
        }else{
          perms<-perm2
          is.base=TRUE
        }
      }
    }
    permsdiff<-perms-perm2
    scan <- addqtl(cross, qtl=intial.model, formula=x, model="normal",
                   method="hk", covar=trt, pheno.col=phe)
    if(is.base){
      diff<-scan
    }else{
      diff<-scan-baseScan
    }
    
    scan <- scan[scan$chr != intial.model$chr | abs(scan$pos-intial.model$pos)>30,]
    diff <- diff[diff$chr != intial.model$chr | abs(diff$pos-intial.model$pos)>30,]
    maxScan<-summary(scan, perms=perms, pvalues=T)
    maxDiff<-summary(diff, perms=permsdiff, pvalues=T)
    maxScan<-data.frame(maxScan[maxScan$lod==max(maxScan$lod),])
    maxDiff<-data.frame(maxDiff[maxDiff$lod==max(maxDiff$lod),])
    colnames(maxScan)<-paste("scan_", colnames(maxScan), sep="")
    colnames(maxDiff)<-paste("diff_", colnames(maxDiff), sep="")
    print(permsdiff)

    print(summary(scan))
    print(summary(diff))
    
    print(maxScan)
    print(maxDiff)
    if(is.base) maxDiff[1,]<-NA
    out<-data.frame(formula=x,phe=phe,maxScan, maxDiff)
    return(out)
  })
  ldply(out, data.frame)
}

full.pipe<-function(cross, phe, covar, verbose=T, nperm=10, ncluster=1){
  require(rrBLUP)
  require(qtl)
  require(parallel)
  #part 1: initial permutations
  if(verbose) cat("running",nperm,"scanone permutations\n")
  set.seed(42)
  addPerms<-scanone(cross=cross, pheno.col=phe, addcovar=covar, intcovar=NULL, perm.strata=covar[,1], n.perm=nperm, n.cluster=ncluster, verbose=F)
  set.seed(42)
  intPerms<-scanone(cross=cross, pheno.col=phe, addcovar=covar, intcovar=covar, perm.strata=covar[,1], n.perm=nperm, n.cluster=ncluster, verbose=F)
  
  #part 2: functions - genoprob2marker
  if(verbose) cat("generating kinship matrix\n")
  genoprob<-calc.genoprob(cross, step=2,  error.prob=0.001, map.function="kosambi")
  grid<-reduce2grid(genoprob)
  datGrid<-genoprob2marker(grid, genotypes=c(-1,0,1))
  covar.name=colnames(covar)
  phenoMatrix<-data.frame(pull.pheno(cross, c("id",phe)),covar)
  datGrid2<-datGrid[,-c(1:3)]
  n<-gsub("X","",colnames(datGrid2))
  k<-A.mat(t(datGrid2))
  colnames(k)<-n; rownames(k)<-n
  
  #part 3: functions - calcPolygenic
  if(verbose) cat("running variance component analysis for polygenic inheritance\n")
  Vs<-calcPolygenic(data=phenoMatrix, pheno=phe, kinshipMatrix=k, covar="trt")
  VsOut_noQTL<-data.frame(t(unlist(Vs))); colnames(VsOut_noQTL)<-c("Vg","Ve"); VsOut_noQTL$phe<-phe
  
  # part 4: run scanone for initial QTL
  if(verbose) cat("running initial scanones\n")
  s1Add<-scanone(cross=cross, pheno.col=phe, method="hk", addcovar=trt, intcovar=NULL)
  s1Int<-scanone(cross=cross, pheno.col=phe, method="hk", addcovar=trt, intcovar=trt)
  
  # part 5: compiling scanone results
  if(verbose) cat("compiling scanone results\n")
  gp<-genoprob2marker(cross)
  comp<-compileScanone(cross=cross, s1Add=s1Add, s1Int=s1Int, addPerms=addPerms, intPerms=intPerms, phe=phe)
  bestForm=comp[["bestForm"]]
  bestMod=comp[["bestMod"]]
  s1Stats=comp[["maxScanData"]]
  
  # part 6: running second polygenic inheritance model (w/ Q1 in model)
  if(verbose) cat("running second polygenic inheritance model (w/ Q1 in model)\n")
  covar2<-data.frame(marker1=as.numeric(as.character(gp[gp$chr==bestMod$chr & gp$pos==bestMod$pos , -c(1:3)])))
  phenoMatrix2<-cbind(covar2,phenoMatrix)
  Vs2<-calcPolygenic(data=phenoMatrix2, pheno=phe, kinshipMatrix=k, covar=c("trt","marker1"))
  VsOut_1QTL<-data.frame(t(unlist(Vs2))); colnames(VsOut_1QTL)<-c("Vg","Ve"); VsOut_1QTL$phe<-phe
  
  covar2<-data.frame(marker=as.numeric(as.character(gp[gp$chr==bestMod$chr & gp$pos==bestMod$pos , -c(1:3)])), trt=covar$trt)
  marker<-data.frame(marker=as.numeric(as.character(gp[gp$chr==bestMod$chr & gp$pos==bestMod$pos , -c(1:3)])))
  # part 7: run additional permutations
  if(verbose) cat("running permutations for a second eQTL ... 1st set")
  set.seed(42)
  perm2<-scanone(cross=cross, pheno.col=phe, addcovar=covar2, intcovar=NULL, perm.strata=covar2$trt, n.perm=nperm, n.cluster=ncluster, verbose=F)
  if(verbose) cat("... 2nd set")
  set.seed(42)
  perm2int<-scanone(cross=cross, pheno.col=phe, addcovar=covar2, intcovar=covar, perm.strata=covar$trt, n.perm=nperm, n.cluster=ncluster, verbose=F)
  set.seed(42)
  if(verbose) cat("... 3rd set")
  perm2epi<-scanone(cross=cross, pheno.col=phe, addcovar=covar2, intcovar=marker, perm.strata=covar$trt, n.perm=nperm, n.cluster=ncluster, verbose=F)
  set.seed(42)
  if(verbose) cat("... 4th set\n")
  perm2intepi<-scanone(cross=cross, pheno.col=phe, addcovar=covar2, intcovar=covar2, perm.strata=covar$trt, n.perm=nperm, n.cluster=ncluster, verbose=F)
  
  # part 8: scanning for second eQTL
  if(verbose) cat("scanning for a second eQTL\n")
  if(bestForm=="y ~ Q1 + Q1*trt + trt"){
    forms<-c("y ~ Q1 + Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q1*trt + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + Q2*trt + trt")
  }else{
    forms<-c("y ~ Q1 + Q2 + trt",
             "y ~ Q1 + Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q2*trt + trt")
  }
  scans<-add1QTL(cross=cross, intial.model=bestMod, new.formula=forms, 
                 perm2intepi = perm2intepi, perm2int = perm2int, perm2epi = perm2epi, perm2 = perm2, 
                 phe = phe, ncluster=ifelse(ncluster>4,4,ncluster))
  bestQ2index<-ifelse(sum(scans$scan_pval==min(scans$scan_pval))==1, 
                      which(scans$scan_pval==min(scans$scan_pval)),
                      which(scans$scan_lod==max(scans$scan_lod)))
  bestQ2<-scans[bestQ2index,c("scan_chr","scan_pos")]
  
  # part 9: running third polygenic inheritance model (w/ Q1 + Q2 in model)
  if(verbose) cat("running third polygenic inheritance model (w/ Q1 in model)\n")
  covar3<-data.frame(marker2=as.numeric(as.character(gp[gp$chr==as.character(bestQ2$scan_chr) & gp$pos==bestQ2$scan_pos , -c(1:3)])))
  phenoMatrix3<-cbind(covar3,phenoMatrix2)
  Vs3<-calcPolygenic(data=phenoMatrix3, pheno=phe, kinshipMatrix=k, covar=c("trt","marker1","marker2"))
  VsOut_2QTL<-data.frame(t(unlist(Vs3))); colnames(VsOut_2QTL)<-c("Vg","Ve"); VsOut_2QTL$phe<-phe
  
  
  return(list(
    scanoneStats=rbind(s1Stats,scans),
    polygenStats=rbind(VsOut_noQTL,VsOut_1QTL,VsOut_2QTL),
    perms=list(addPerms=addPerms, intPerms=intPerms, perm2=perm2,perm2int=perm2int,perm2epi=perm2epi,perm2intepi=perm2intepi)
  ))
}
