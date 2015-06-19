stepwiseStats<-function(cross, model.in, phe, covar=NULL, ci.method="drop", drop=1.5, prob=.95, plot=FALSE, printout=TRUE){
  
  if(class(cross)[1]=="riself" | class(cross)[1]=="bc"){
    statcols<-c("phenotype", "chromosome", "position",  "df", "type3SS", "LOD", "perc.var", "Fstat", "P.chi2", "P.F",
                "effect.estimate", "effect.SE", "effect.t",
                "lowCImarker", "hiCImarker", "lowCIpos", "hiCIpos")
  }else{
    if(class(cross)[1]=="f2"){
      statcols<-c("phenotype","chromosome","position","df","type3SS","LOD","perc.var","Fstat","P.chi2","P.F",
                  "est.dom", "SE.dom", "t.dom", "est.add", "SE.add", "t.add",
                  "lowCImarker", "hiCImarker","lowCIpos", "hiCIpos")
    }else{
      stop("stepwiseStats is only implemented for f2, bc and ril experimental designs")
    }
  }

  stepout<-model.in
  
  if(nqtl(stepout)>0){
    
    
    if(plot){
      plotLodProfile(stepout,main=paste(phe,"formula: ", formula(stepout)))
    }
    
    #fit the qtl model
    fit<-fitqtl(cross,
                pheno.col=phe,
                qtl=stepout,
                formula=formula(stepout),
                get.ests=T,dropone=T,covar=covar,
                method="hk")
    
    #determine the number of terms/qtl in the model
    nqtls<-nqtl(stepout)
    nterms<-sum(countqtlterms(formula(stepout), ignore.covar=F)[c(1,4)])
    ncovar<-length(covar)
    
    #extract information if the qtl model has multiple qtls/covariates
    if(nterms>1){
      ests_out<-summary(fit)$ests[2:(nterms+1),1:3]
      cis<-data.frame()
      #calculate confidence intervals for each qtl
      for (j in 1:nqtls){
        attr(stepout, "lodprofile")
        if(ci.method=="drop"){ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=drop)
        }else{
          if(ci.method=="bayes"){ciout<-bayesint(stepout,qtl.index=j, expandtomarkers=F, prob=prob)}
        }
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      chrs<-stepout$chr
      pos<-stepout$pos
      #if covariates, add in NAs into the dataframe
      if(ncovar>0) {
        covar.ci<-(rep(NA,4))
        cis<-rbind(covar.ci, cis)
        chrs<-c(names(covar), stepout$chr)
        pos<- c(names(covar), stepout$pos)
      }
      #if epistasis, add NAs into the dataframe
      nepi<-as.numeric(countqtlterms(formula(stepout))[4])
      if(nepi>0) {
        epi.ci<-data.frame()
        for (j in 1:nepi){
          epi.ci.out<-(rep(NA,4))
          cis<-rbind(cis,epi.ci.out)
        }
        chrs<-c(chrs,rep(NA,nepi))
        pos<- c(pos,rep(NA,nepi))
      }
      
      if(class(cross)[1] %in% c("riself","bc")){
        ests_out<-summary(fit)$ests[2:(nterms+1),1:3]
        stats<-data.frame(rep(phe,nterms),chrs,pos,fit$result.drop[1:nterms,],ests_out,cis)
        colnames(stats)<-statcols
      }else{
        ests.all<-summary(fit)$ests
        rows.dom<-grep("d",rownames(ests.all))
        rows.add<-grep("a",rownames(ests.all))
        dom.ests_out<-ests.all[rows.dom,1:3]
        add.ests_out<-ests.all[rows.add,1:3]
        ests.out<-data.frame(dom.ests_out,add.ests_out)
        ests.out<-ests.out[1:nqtls,]
        if(nepi>0) {
          epi.ests<-ests.out[1:nepi,]
          epi.ests[!is.na(epi.ests)]<-NA
          ests.out<-rbind(ests.out, epi.ests)
        }        
        colnames(ests.out)<-c("est.dom","SE.dom","t.dom","est.add","SE.add","t.add")
        stats<-data.frame(rep(phe,nterms),chrs,pos,fit$result.drop[1:nterms,],ests.out,cis)
        colnames(stats)<-statcols
      }
    }else{
      cis<-data.frame()
      for (j in 1:nqtls){
        if(ci.method=="drop"){ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=drop)
        }else{
          if(ci.method=="bayes"){ciout<-bayesint(stepout,qtl.index=j, expandtomarkers=F, prob=prob)}
        }
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      if(class(cross)[1]=="riself" | class(cross)[1]=="bc"){
        ests_out<-summary(fit)$ests[2:(nqtls+1),1:3]
        stats<-data.frame(c(phe,stepout$chr[1],stepout$pos[1],
                            as.numeric(fit$result.full[1,c(1,2,4,5,3,6,7)]),as.numeric(ests_out),cis[1,]))
      }else{
        ests_out<-summary(fit)$ests
        rows.dom<-grep("d",rownames(ests_out))
        rows.add<-grep("a",rownames(ests_out))
        if(nterms==1){
          dom.ests_out<-ests_out[rows.dom,1:3]; names(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
          add.ests_out<-ests_out[rows.add,1:3]; names(add.ests_out)<-c("est.add","SE.add","t.add")
          ests.out<-data.frame(t(data.frame(c(dom.ests_out,add.ests_out))))
        }else{
          if(nterms>nqtls){
            rows.epi<-grep(":",rownames(ests_out))
            rows.dom<-rows.dom[!rows.dom %in% rows.epi]
            rows.add<-rows.add[!rows.add %in% rows.epi]
            dom.ests_out<-data.frame(ests_out[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
            add.ests_out<-data.frame(ests_out[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
            ests.out<-cbind(dom.ests_out,add.ests_out)
            #add rows of nas to data frame for epis
            n.epi<-nterms-nqtls
            for (i in 1:n.epi) ests.out<-rbind(ests.out,NA)     
          }else{
            dom.ests_out<-data.frame(ests_out[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
            add.ests_out<-data.frame(ests_out[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
            ests.out<-cbind(dom.ests_out,add.ests_out)
          }
        }
        stats<-data.frame(c(phe,stepout$chr[1],stepout$pos[1],
                            as.numeric(fit$result.full[1,c(1,2,4,5,3,6,7)]),
                            ests.out,
                            cis[1,]))
      }
      colnames(stats)<-statcols; rownames(stats)<-stepout$name
    }
    stats$id<-rownames(stats)
    stats<-stats[,c(which(colnames(stats)=="id"),1:(length(colnames(stats))-1))]
    return(stats)
  }else{
    cat(" ******* phenotype",phe,"has a null qtl model","\n")
    stats<-NULL
  }
}
