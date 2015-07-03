#function to call alleles from ASE RNA-seq data
#requires the name of the gene and two files with the ASE data with colnames as genes
genobyMclust<-function(gene.name, hal.file, fil.file, ids=c("id","type") ,
                       logtrans=TRUE, n.clust=3,p.tol=1e-10, plot=FALSE){
  geno.in<-merge(hal.file[,c(ids,gene.name)],fil.file[,c(ids,gene.name)], by=ids)
  colnames(geno.in)[3:4]<-c("hal","fil")
  if(logtrans){
    geno.in$hal<-log10(geno.in$hal+1)
    geno.in$fil<-log10(geno.in$fil+1)
  }
  clust<-Mclust(geno.in[,c("hal","fil")],G=n.clust, prior = priorControl())
  geno.in$group<-clust$classification
  geno.in$uncertainty<-clust$uncertainty
  geno.in$group[geno.in$uncertainty>=p.tol]<-NA
  F1.group<-unique(geno.in$group[geno.in$type=="F1"])
  if(length(F1.group)>1){
    tab<-table(geno.in$group[geno.in$type=="F1"])
    tab.names<-names(tab)
    F1.group<-tab.names[which(as.numeric(tab)==max(as.numeric(tab)))]
    cat(gene.name,":",sum(as.numeric(tab)[as.numeric(tab)!=max(as.numeric(tab))]), "samples miscategorized for F1 \n")
  }
  hal.group<-names(table(geno.in$group[geno.in$type=="HA"]))
  if(length(hal.group)>1){
    tab<-table(geno.in$group[geno.in$type=="HA"])
    tab.names<-names(tab)
    F1.group<-tab.names[which(as.numeric(tab)==max(as.numeric(tab)))]
    cat(gene.name,":",sum(as.numeric(tab)[as.numeric(tab)!=max(as.numeric(tab))]), "samples miscategorized for Hal \n")
  }
  fil.group<-names(table(geno.in$group[geno.in$type=="FI"]))
  if(length(fil.group)>1){
    tab<-table(geno.in$group[geno.in$type=="FI"])
    tab.names<-names(tab)
    F1.group<-tab.names[which(as.numeric(tab)==max(as.numeric(tab)))]
    cat(gene.name,":",sum(as.numeric(tab)[as.numeric(tab)!=max(as.numeric(tab))]), "samples miscategorized for Fil \n")
  }
  geno.in$allele<-ifelse(is.na(geno.in$group),NA,
                         ifelse(geno.in$group==fil.group,"fil",
                                ifelse(geno.in$group==hal.group,"hal",
                                       ifelse(geno.in$group==F1.group,"het",NA))))
  if(plot){
    p1<-ggplot(geno.in, aes(x=hal, y=fil,col=as.factor(type)))+
      geom_point()+
      scale_y_continuous(name="FIL-allele-specific expression (log10)")+
      scale_x_continuous(name="HAL-allele-specific expression (log10)")+
      theme_bw()+
      ggtitle("expression categorization by experimental unit")
    p2<-ggplot(geno.in, aes(x=hal, y=fil,col=as.factor(allele)))+
      geom_point()+
      scale_y_continuous(name="FIL-allele-specific expression (log10)")+
      scale_x_continuous(name="HAL-allele-specific expression (log10)")+
      theme_bw()+
      ggtitle("expression categorization by  genotype")
    print(multiplot(p1, p2, cols=1))
  }
  
  geno.out<-geno.in[,c("id","allele")]
  colnames(geno.out)[2]<-gene.name
  return(geno.out)
}