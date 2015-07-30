parseLP<-function(cross, models, formulae,covar=NULL, allchr=F){
  if(is.null(names(models))) stop("model names must match phenotype names")
  if(!identical(names(models), names(formulae))) stop("model names must match formulae names")
  phenotypes<-names(models)
  refs<-lapply(phenotypes, function(x) {
    if(allchr & length(chrnames(cross))!= length(unique(models[[x]]$chr))){
      chrs<-chrnames(cross)[which(!chrnames(cross) %in% models[[x]]$chr)]
      poss<-rep(0, length(chrs))
      models[[x]]<-addtoqtl(cross=cross, qtl=models[[x]], chr=chrs, pos=poss,
                            qtl.name=paste("fake",chrs, sep="_"))
      formulae[[x]]<-paste(formulae[[x]], paste("fake",chrs, sep="_"), sep=" + ")
    }
    ref<-refineqtl(cross, qtl=models[[x]], formulae[[x]], pheno.col=x, covar=covar, 
                   verbose=F, method="hk")
    lp<-ldply(attr(ref, "lodprofile"), data.frame)
    colnames(lp)<-c("qtlname","chr","pos","lodProfile")
    lp$stLodProfile<-lp$lodProfile/max(lp$lodProfile)
    return(lp)
  })
  names(refs)<-phenotypes
  dat<-ldply(refs,data.frame)
  colnames(dat)[1]<-"phenotype"
  dat
}

plotMultiLP<-function(parseLPoutput, title="multipleLodProfile Plot",...){
  library(ggplot2)
  parseLPoutput$chr<-as.numeric(as.character(parseLPoutput$chr))
  ggplot(data=parseLPoutput, aes(x=pos, y=lodProfile, col=phenotype, group=interaction(qtlname, chr, phenotype)))+
    geom_line(...)+
    facet_grid(.~chr, scales="free", space="free")+
    theme_bw()+
    theme(
      panel.background = element_rect(fill = "ghostwhite")
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,strip.background = element_blank()
      ,axis.text.x = element_blank()
      ,axis.ticks.x = element_blank()
      ,axis.text.y = element_text(colour="black")
    ) +
    theme(strip.text.y = element_text(angle = 0,hjust=0))+
    ggtitle("chromosome")+
    scale_x_continuous("")+
    scale_y_continuous("LOD")
}

makeLpPlot<-function(cross, models, formulae, verbose=T, covar=NULL, allchr=F, title="multipleLodProfile Plot",...){
  if(verbose) cat("calculating LOD profiles via refineQTL\n")
  parsed<-parseLP(cross, models, formulae, covar, allchr)
  if(verbose) cat("generating the plot\n")
  print(plotMultiLP(parsed, title,...))
}