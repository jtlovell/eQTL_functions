bp2cm<-function(cross,geneID=gff$A.id,chrID=gff$chr, bpPos=gff$start, markerBP=NULL, markerID=NULL, splitByCentromere=TRUE){
  #create an annotation file with the necessary columns
  if(!is.null(markerBP) & !is.null(markerID)){
    gff<-data.frame(geneID,markerID=NULL,chrID,bpPos)
    markers<-data.frame(gene.id=NULL,markerID, chrID=NULL, bpPos=markerBP)
    gff<-rbind(gff,markers)
  }else{
    if(is.null(markerID)) markerID<-geneID
    gff<-data.frame(geneID,markerID,chrID,bpPos)
  }
  
  require(splines)
  require(car)
  
  par(mfrow=c(2,1))
  all.dat<-data.frame()
  
  for(i in chrnames(cross)){
    m<-pull.map(cross, chr=i, as.table=T)
    colnames(m)[colnames(m)=="chr"]<-"lg"
    m$markerID<-rownames(m)
    
    g<-merge(m,gff, by="markerID")
    
    bestChr<-as.data.frame(table(as.character(g[,"chrID"])))
    bestChr<-as.character(bestChr[,1][bestChr$Freq==max(bestChr$Freq)])
    
    if(length(unique(g$chrID))!=1){
      cat("markers on lg", i, "are on multiple physical chromosomes...\n",
          "running predicted positions only for markers found on physical chromosome", bestChr)
      g<-g[g$chrID==bestChr,]
    }
    g<-g[order(g$bpPos),]
    plot(g$bpPos,g$pos, main=paste("mapping vs. physical position of markers on chr",i), cex=.2, pch=19, ylab="mapping position (cM)", xlab="physical postion (bp)")
    if(splitByCentromere){
      diffs<-diff(g$bpPos)
      ot<-outlierTest(lm(diff(g$bpPos)~diff(g$pos)),cutoff=Inf, n.max=10)
      cmIndex<-as.numeric(names(ot[[2]][ot[[2]]<0.1]))
      centromere<-g$pos[c(min(cmIndex),(max(cmIndex)+1))]
      rect(ybottom=centromere[1], ytop=centromere[2],xleft=0,xright=max(g$bpPos), col=rgb(1,0,0,.2))
      splits<-list(top=c(0,(min(centromere))),
                   centromere=range(centromere),
                   bottom=c((max(centromere)), max(g$pos)))
      out.all<-data.frame()
      for(j in 1:3){
        g1<-g[g$pos>=splits[[j]][1] & g$pos<=splits[[j]][2],]
        bps<-c(min(g1$bpPos), max(g1$bpPos))
        if(j==2){
          gff1<-gff[gff$chr==bestChr & gff$bpPos>=bps[1] & gff$bpPos<=bps[2],]
        }else{
          gff1<-gff[gff$chr==bestChr & gff$bpPos>bps[1] & gff$bpPos<bps[2],]
        }        
        geneID<-gff1$geneID
        bpPos<-data.frame(bpPos=gff1$bpPos)
        if(nrow(g1)<4){
          mod<-with(g1, lm(pos~bpPos))
          abline(mod)
          out<-data.frame(predict(mod,bpPos),bpPos,geneID)
          colnames(out)[1]<-"imputed.cm"
        }else{
          mod<-smooth.spline(g1$bpPos,g1$pos, df=sqrt(nrow(g1)))
          lines(mod)
          out<-data.frame(predict(mod,x=bpPos),geneID)
          colnames(out)[1:2]<-c("bpPos","imputed.cm")
        }
        max.marker<-max(g1$pos)
        min.marker<-min(g1$pos)
        out$imputed.cm[out$imputed.cm>max.marker]<-max.marker
        out$imputed.cm[out$imputed.cm<min.marker]<-min.marker
        out.dat<-merge(gff1,out, by=c("geneID","bpPos"))
        out.all<-rbind(out.all,out.dat) 
      }
      with(out.all, plot(bpPos,imputed.cm, main=paste("imputed positions for each gene on chr", i), cex=.2, pch=19))
    }else{
      g1<-g
      gff1<-gff[gff$chr==bestChr,]
      geneID<-gff1$geneID
      bpPos<-data.frame(bpPos=gff$bpPos[gff$chr==bestChr])
      mod<-smooth.spline(g1$bpPos,g1$pos, df=sqrt(nrow(g1)))
      lines(mod)
      out<-data.frame(predict(mod,x=bpPos),geneID)
      colnames(out)[1:2]<-c("bpPos","imputed.cm")
      out.dat<-merge(gff1,out, by=c("geneID","bpPos"))
      out.all<-rbind(out.all,out.dat) 
      with(out.all, plot(bpPos,imputed.cm, main=paste("imputed positions for each gene on chr", i)))
    }
    all.dat<-rbind(all.dat,out.all)
  }
  return(all.dat)
}