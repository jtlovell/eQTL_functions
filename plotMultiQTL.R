plotMultiQTL<-function(cross, stats=NULL, phes=NULL,chrs=NULL, peak=NULL, right=NULL, left=NULL, col=NULL, 
                       chr.subset=NULL, ylabelcex=NULL, rugsize=NULL, cex=NULL, pch=19,lty=1,lwd=1,
                       plotQTLdensity=TRUE, binwidth=1, adj.ylabsize=TRUE,
                       colbychr=TRUE, palette=rainbow, showConfidenceInterval=TRUE, showPointEstimate=TRUE,
                       outline=FALSE, background=TRUE,plotNullPheno=FALSE, setmargin=NULL, ...){

  if(is.null(stats)){
    stats<-data.frame(phes, chrs, peak, right, left)
    colnames(stats)<-c("phenotype","chromosome","position","lowCIpos","hiCIpos")
  }
  if(plotNullPheno){
    stats<-stats[stats$chromosome %in% chr.subset,]
  }
  nqtl<-dim(stats)[1]
  nphes<-length(unique(stats$phenotype))
  
  if(is.null(cex)){
    cex<-(1/(1+(nphes*.01)))
  }
  if(!is.integer(col) & !is.character(col)){
    col<-"black"
  }
  for(i in c(lty, lwd, pch)) {
    if(is.null(i)){
      assign(i,1)
    }
  }
  for(i in c("pch", "lty", "lwd", "cex", "col")){
    dat<-get(i)
    if(length(i)!=1 & length(i)!=nrow(stats)){
      cat("warning: if specified,",i,"must be a vector of length 1 or ", nrow(stats), "...\n setting", i, "to",dat[1])
      stats[,i]<-dat[1]
    }else{
      stats[,i]<-dat
    }
  }
  #simplify the dataset
  stats1<-stats[complete.cases(stats),c("phenotype","chromosome","position","lowCIpos","hiCIpos",
                                        "col","pch","lty","lwd","cex")]
  for(i in c("position","lowCIpos","hiCIpos","pch","lty","lwd","cex")) stats1[,i]<-as.numeric(as.character(stats1[,i]))
  
  #convert colors to chromosomes if specified
  if(colbychr){
    pal<-palette(nchr(cross))
    stats1$col<-as.character(unlist(sapply(stats1$chromosome, function(x) pal[which(chrnames(cross)==x)])))
  }
  
  #prepare the rug
  if(is.null(rugsize)){
    rugsize<-(1/(nphes^2))+.01
  }
  if(adj.ylabsize){
    ylab.adj<-(1/(nphes))+.1
    par(mar=c(5,4,4,2)+.1)
  }else{
    ylab.adj<-.5
    par(mar=c(5,a,4,2)+.1)
  }
  #set plotting window
  if(is.null(setmargin)){
    a<-(sqrt(2*max(sapply(as.character(unique(stats$phenotype)),nchar)))*(ylab.adj^2))+4
    par(mar=c(5,a,4,2)+.1)
  }else{
    par(mar=setmargin)
  }
  #prepare the map
  if(!is.null(chr.subset)){
    map<-pull.map(cross, as.table=T, chr=chr.subset)
    totlen<-sum(chrlen(cross)[which(chrnames(cross) %in% chr.subset)])
    chrn<-chrnames(cross)[which(chrnames(cross) %in% chr.subset)]
    gapsize=totlen/40
    gaps<-c(0,cumsum(rep(gapsize,length(chrn)-1)))
    chrslens<-c(0, cumsum(chrlen(cross)[which(chrnames(cross) %in% chr.subset)][-length(chrn)]))
    corr<-gaps+chrslens
    corr<-data.frame(chrn, chrslens,corr, stringsAsFactors=F)
  }else{
    map<-pull.map(cross, as.table=T)
    totlen<-sum(chrlen(cross))
    chrn<-chrnames(cross)
    gapsize=totlen/40
    gaps<-c(0,cumsum(rep(gapsize,length(chrn)-1)))
    chrslens<-c(0, cumsum(chrlen(cross)[-length(chrn)]))
    corr<-gaps+chrslens
    corr<-data.frame(chrn, chrslens,corr, stringsAsFactors=F)
  }
 

  
  #convert positions to plotting positions
  dat<-stats1
  map2<-map
  for (i in corr$chrn){
    for(j in c("position", "lowCIpos", "hiCIpos")){
      dat[dat$chromosome==i,j]<-dat[dat$chromosome==i,j]+corr$corr[corr$chrn==i]
    }
    map2[map2$chr==i,"pos"]<-map2[map2$chr==i,"pos"]+corr$corr[corr$chrn==i]
  }
  
  #add space at the top of the plotting window if needed for density plots
  if(plotQTLdensity){
    plot(0,0, ylim=c(0,(nphes+(.1*nphes)+1)), xlim=c(0,max(map2$pos)), type="n", bty="n", yaxt="n",xaxt="n",ylab="", ...)
  }else{
    plot(0,0, ylim=c(0,nphes), xlim=c(0,max(map2$pos)), type="n", bty="n", yaxt="n",xaxt="n",ylab="",...)
  }
  
  #add the x axis
  if(outline) box()
  tloc<-ddply(map2,"chr",summarize,mean=mean(pos))$mean
  axis(side=1, at=tloc, labels=chrn)
  
  phes<-as.character(unique(stats$phenotype))
  
  #add the y axis
  
  for(i in 1:nphes){
    p<-phes[i]
    if(is.null(chr.subset)){
      tp<-dat[dat$phenotype==p,]
    }else{
      tp<-dat[dat$phenotype==p & dat$chromosome %in% chr.subset,]
    }
    if(is.null(ylabelcex)){
      axis(side=2, at=i, labels=p, las=2, cex.axis=ylab.adj)
    }else{
      axis(side=2, at=i, labels=p, las=2, cex.axis=ylabelcex)
    }
    if(nrow(tp)==0){
      next
    }else{
      if(showConfidenceInterval){
        segments(tp$lowCIpos,i,tp$hiCIpos,i,col=tp$col, lty=tp$lty, lwd=tp$lwd)
      }
      if(showPointEstimate){
        points(tp$position,rep(i,length(tp$position)), col=tp$col,cex=tp$cex, pch=tp$pch)
      }
    }    
  }

  rug(map2$pos, ticksize=rugsize)
  if(background){
    if(colbychr){
      min.pos<-tapply(map2$pos, map2$chr, min)
      max.pos<-tapply(map2$pos, map2$chr, max)
      if(is.null(chr.subset)){
        pal2<-palette(nchr(cross), alpha=.05)
      }else{
        pal2<-palette(nchr(cross), alpha=.05)[which(chrnames(cross) %in% chr.subset)]
      }
      
      rect(min.pos,rep(0,nchr(cross)),max.pos, rep(nphes), col=pal2, border = pal2)
    }else{
      min.pos<-tapply(map2$pos, map2$chr, min)
      max.pos<-tapply(map2$pos, map2$chr, max)
      if(is.null(chr.subset)){
        rect(min.pos,rep(0,nchr(cross)),max.pos, rep(nphes), col=rgb(0,0,0,.05), border = rgb(0,0,0,.05))
      }else{
        rect(min.pos,rep(0,length(chr.subset)),max.pos, rep(nphes), col=rgb(0,0,0,.05), border = rgb(0,0,0,.05))
      }
    }
  }
  if(plotQTLdensity){
    if(!is.null(chr.subset)){
      denpos<-density(dat$position[dat$chromosome %in% chr.subset],bw=binwidth)$x
      den<-density(dat$position[dat$chromosome %in% chr.subset],bw=binwidth)$y
    }else{
      denpos<-density(dat$position,bw=binwidth)$x
      den<-density(dat$position,bw=binwidth)$y
    }
    den<-den/max(den)
    den<-den*(.1*nphes)+1
    den<-den+nphes
    lines(denpos,den)
  }
}
  