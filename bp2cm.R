bp2cm<-function(cross,geneID=gff$A.id,chrID=gff$chr, bpPos=gff$start, markerBP=NULL, markerID=NULL, splitByCentromere=TRUE,minCI=1){
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
    #plot(g$bpPos,g$pos, main=paste("mapping vs. physical position of markers on chr",i), cex=.2, pch=19, ylab="mapping position (cM)", xlab="physical postion (bp)")
    if(splitByCentromere){
      diffs<-diff(g$bpPos)
      ot<-outlierTest(lm(diff(g$bpPos)~diff(g$pos)),cutoff=Inf, n.max=10)
      cmIndex<-as.numeric(names(ot[[2]][ot[[2]]<0.1]))
      centromere<-g$pos[c(min(cmIndex),(max(cmIndex)+1))]
      #rect(ybottom=centromere[1], ytop=centromere[2],xleft=0,xright=max(g$bpPos), col=rgb(1,0,0,.2))
      centromere.range<-range(centromere)
      if(min(centromere)<5){
        centromere.range[1]<-0
        arms<-2:3
      }else{
        arms<-1:3
      }
      splits<-list(top=c(0,(min(centromere))),
                   centromere=centromere.range,
                   bottom=c((max(centromere)), max(g$pos)))
      out.all<-data.frame()
      par(mfrow=c(3,1))
      
      for(j in arms){
        g1<-g[g$pos>=splits[[j]][1] & g$pos<=splits[[j]][2],]
        bps<-c(min(g1$bpPos), max(g1$bpPos))
        if(j==2){
          gff1<-gff[gff$chr==bestChr & gff$bpPos>=bps[1] & gff$bpPos<=bps[2],]
        }else{
          gff1<-gff[gff$chr==bestChr & gff$bpPos>bps[1] & gff$bpPos<bps[2],]
        }        
        geneID<-gff1$geneID
        bpPos<-data.frame(bpPos=gff1$bpPos)
        df1<-ifelse(floor(nrow(g1)/6)==0,1,floor(nrow(g1)/6))
        mod <- lm(pos ~ ns(bpPos, df=df1) , data=g1)
        plot(g1$bpPos,g1$pos)
        pred<-predict(mod, newdata=data.frame(bpPos=gff1$bpPos), interval = "prediction", level=.95)
        pred<-cbind(data.frame(bpPos=gff1$bpPos), pred)
        colnames(pred)[2]<-"pos"
        out<-pred[order(pred$bpPos),]
        for(k in c("pos","lwr","upr")) lines(out[,"bpPos"],out[,k], col="red")
        colnames(out)[1:2]<-c("bpPos","imputed.cm")
        out$geneID<-geneID
        
        max.marker<-max(g1$pos)
        min.marker<-min(g1$pos)
        out$imputed.cm[out$imputed.cm>max.marker]<-max.marker
        out$imputed.cm[out$imputed.cm<min.marker]<-min.marker
        out$lwr[out$lwr<0]<-0
        if(j==3){
          out$upr[out$upr>max(g1$pos)]<-max(g1$pos)
        }
        out$arm<-ifelse(j==1,"proximate",ifelse(j==2, "centromere", "distal"))
        out$chromosome<-i
        out.all<-rbind(out.all,out) 
      }
      par(mfrow=c(1,1))
      plot(g$bpPos,g$pos, main=paste("mapping vs. physical position of markers on chr",i), cex=.2, pch=19, ylab="mapping position (cM)", xlab="physical postion (bp)")
      if(min(centromere)<5){
        plwr.max=0
        pupr.max=min(out.all$upr[out.all$arm=="distal"])
        p.max=0
      }else{
        pupr.max<-max(out.all$upr[out.all$arm=="proximate"])
        plwr.max<-max(out.all$lwr[out.all$arm=="proximate"])
        p.max<-max(out.all$bpPos[out.all$arm=="proximate"])
      }

      qupr.min<-min(out.all$upr[out.all$arm=="distal"])
      qlwr.min<-min(out.all$lwr[out.all$arm=="distal"])
      q.min<-min(out.all$bpPos[out.all$arm=="distal"])
      lwr<-data.frame(lwr=c(plwr.max,qlwr.min),bp=c(p.max,q.min))
      upr<-data.frame(upr=c(pupr.max,qupr.min),bp=c(p.max,q.min))
      
      upr.mod<-lm(upr~bp, data=upr)
      lwr.mod<-lm(lwr~bp, data=lwr)
      upr.pred<-predict(upr.mod, data.frame(bp=out.all$bpPos[out.all$arm=="centromere"]))
      lwr.pred<-predict(lwr.mod, data.frame(bp=out.all$bpPos[out.all$arm=="centromere"]))
      out.all$upr[out.all$arm=="centromere"]<-upr.pred
      out.all$lwr[out.all$arm=="centromere"]<-lwr.pred
      if(!is.na(minCI)){
        out.all$upr[out.all$upr-out.all$imputed.cm<1]<-out.all$imputed.cm[out.all$upr-out.all$imputed.cm<1]+minCI
        out.all$lwr[out.all$imputed.cm-out.all$lwr<1]<-out.all$imputed.cm[out.all$imputed.cm-out.all$lwr<1]-minCI
      }
      out.all<-out.all[order(out.all$bpPos),]
      #with(out.all, plot(bpPos,imputed.cm, type="p", cex=.1))
      with(out.all, lines(bpPos,lwr, type="l",col="blue"))
      with(out.all, lines(bpPos,upr, type="l",col="red"))
      #with(out.all, plot(bpPos,imputed.cm, main=paste("imputed positions for each gene on chr", i), cex=.2, pch=19))
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