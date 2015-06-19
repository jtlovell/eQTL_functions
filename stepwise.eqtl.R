stepwise.eqtl<-function(cross, covar,
                           phe, pos.gene, chr.gene, 
                           model.type="normal", method="hk",
                           penalties=c(3,3,3),
                           nqtl=4, binary.tolerance=.0001){

  #Part 1: set up starting model for stepwise
  init.model<-makeqtl(cross, chr=chr.gene, pos=pos.gene, what="prob")
  #Part 2: Stepwise model selection
  mod<-stepwiseqtl(cross, pheno.col=phe, 
                   covar=covar, qtl=init.model,
                   method=method,model="normal",scan.pairs=F,
                   additive.only=F, max.qtl=nqtl,
                   penalties=penalties, keeplodprofile=T, verbose=F,tol=binary.tolerance)
  stepout<-mod
  #Part 3: output if NULL QTL from stepwise
  if(nqtl(stepout)==0){
    output<-list(NULL,NULL)
  }else{
    #Part 4: fit QTL model
    nqtls<-nqtl(stepout)
    nterms<-sum(countqtlterms(formula(stepout), ignore.covar=F)[c(1,4)])
    ncovar<-length(covar)
    n.epi<-nterms-(ncovar + nqtl)
    fit<-fitqtl(cross,
                pheno.col=phe,
                qtl=stepout,
                formula=formula(stepout),
                get.ests=T,dropone=T,covar=covar,
                method="hk", model=model.type)
    
    #part 5 output estimates
    ests.all<-summary(fit)$ests
    rows.dom<-intersect(grep("d",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE))
    rows.add<-intersect(intersect(grep("a",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE)), 
                        grep(colnames(covar),rownames(ests.all),invert=TRUE))
    rows.cov<-which(rownames(ests.all) %in% colnames (covar))
    if(nqtls==1){
      covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
      colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
      dom.ests_out<-data.frame(t(ests.all[rows.dom,1:3])); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
      add.ests_out<-data.frame(t(ests.all[rows.add,1:3])); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
      qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
      ests.out<-rbind(covar.ests_out,qtl.ests_out)
    }else{
      covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
      colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
      rownames(covar.ests_out)<-colnames(covar)
      dom.ests_out<-data.frame(ests.all[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
      add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
      qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
      ests.out<-rbind(covar.ests_out,qtl.ests_out)
    }
    #add nas if epistasis
    if(n.epi>0){
      for (i in 1:n.epi) ests.out<-rbind(ests.out,NA)     
    }
    
    #part 6 get confidence intervals
    cis<-data.frame()
    #calculate confidence intervals for each qtl
    for (j in 1:nqtls){
      ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=1.5)
      lowmarker<-rownames(ciout)[1]
      highmarker<-rownames(ciout)[3]
      lowposition<-ciout[1,2]
      highposition<-ciout[3,2]
      cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
    }
    #add nas if epistasis
    if(n.epi>0){ for (i in 1:n.epi) cis<-rbind(cis,NA)}
    if(ncovar>0){ for(i in 1:ncovar) cis<-rbind(rep(NA,4),cis)}
    
    #part 7 full model stats out
    fit.out<-data.frame(fit$result.drop)
    pos.out<-sapply(rownames(fit.out), function(x) strsplit(x,"@")[[1]][2])
    pos.out[grep(":",pos.out)]<-NA
    fit.out$pos<-pos.out
    chr.out<-sapply(rownames(fit.out), function(x) strsplit(x,"@")[[1]][1])
    chr.out[grep(":",pos.out)]<-NA
    fit.out$chr<-chr.out
    fit.out$phenotype<-phe

    #part 8 check for covariate interactions
    covint<-addcovarint(cross, pheno.col=phe, 
                        covar=covar, icovar=colnames(covar),
                        qtl=mod, formula=formula(mod),
                        method=method, model=model.type,tol=binary.tolerance)
    covint<-data.frame(covint)[,c("LOD","X.var","Pvalue.Chi2.")]  
    if(n.epi>0){ for (i in 1:n.epi) covint<-rbind(covint,NA)}
    if(ncovar>0){ for(i in 1:ncovar) covint<-rbind(NA,covint)}
    colnames(covint)<-c("covint.LOD","covint.X.var","covint.Pva")
    
    all.out<-cbind(fit.out,cis,ests.out,covint)
    return(list(stepout,all.out))
  }
}