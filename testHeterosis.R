testHeterosis<-function(p1vf1,p2vf1,p1vp2, genoID=c("f","f1","h")){
  dat1<-data.frame(p1vf1)
  id1<-paste(genoID[1],"v",genoID[2], sep="")
  colnames(dat1)<-paste(colnames(dat1),id1, sep="_")
  dat1$gene<-rownames(dat1)
  
  id2<-paste(genoID[3],"v",genoID[2], sep="")
  dat2<-data.frame(p2vf1)
  colnames(dat2)<-paste(colnames(dat2),id2, sep="_")
  
  id3<-paste(genoID[1],"v",genoID[3], sep="")
  dat3<-data.frame(p1vp2)
  colnames(dat3)<-paste(colnames(dat3),id3, sep="_")
  
  dat<-cbind(dat1, dat2, dat3)
  
  d<-dat[,c("gene",colnames(dat)[c(grep("log2FoldChange", colnames(dat)),grep("stat", colnames(dat)), grep("padj", colnames(dat)))])]
  for(i in colnames(d)[grep("padj", colnames(d))]) d[is.na(d[,i]),]<-1
  
  d$hp<-ifelse(d$log2FoldChange_fvh>=0, genoID[1], genoID[3])
  d1<-d[d$hp == genoID[1],]
  d2<-d[d$hp == genoID[3],]
  d3<-rbind(d1,d2)
  lfcID1<-paste("log2FoldChange",id1, sep="_")
  lfcID2<-paste("log2FoldChange",id2, sep="_")
  lfcID3<-paste("log2FoldChange",id3, sep="_")
  statID1<-paste("stat",id1, sep="_")
  statID2<-paste("stat",id2, sep="_")
  statID3<-paste("stat",id3, sep="_")
  pID1<-paste("padj",id1, sep="_")
  pID2<-paste("padj",id2, sep="_")
  pID3<-paste("padj",id3, sep="_")
  
  d3$hpvF1_lfc<-c(d1[,lfcID1], d2[,lfcID2])
  d3$hpvF1_p<-c(d1[,pID1], d2[,pID2])
  d3$hpvF1_stat<-c(d1[,statID1], d2[,statID2])
  
  d3$lpvF1_lfc<-c(d1[,lfcID2], d2[,lfcID1])
  d3$lpvF1_p<-c(d1[,pID2], d2[,pID1])
  d3$lpvF1_stat<-c(d1[,statID2], d2[,statID1])
  
  d3$hpvlp_lfc<-c(d1[,lfcID3], (-1*d2[,lfcID3]))
  d3$hpvlp_p<-c(d1[,pID3], d2[,pID3])
  d3$hpvlp_stat<-c(d1[,statID3], (-1*d2[,statID3]))
  
  d3$f1pos<-with(d3, ifelse(hpvF1_lfc>0 & lpvF1_lfc<0, "mid", 
                            ifelse(hpvF1_lfc<=0,"high","low")))
  
  d3$heterosis_cat<-with(d3, ifelse(hpvF1_p>0.05 & lpvF1_p>0.05 & hpvlp_p>0.05,"noDE",
                                    ifelse(f1pos=="mid" & hpvF1_p<=0.05 & lpvF1_p<=0.05,"additive",
                                           ifelse(f1pos %in% c("mid","low") & hpvF1_p<=0.05 & lpvF1_p>0.05,"below mid-parent heterosis",
                                                  ifelse(f1pos %in% c("mid","high") & hpvF1_p>0.05 & lpvF1_p<=0.05,"above mid-parent heterosis",
                                                         ifelse(f1pos=="high" & hpvF1_p<=0.05,"above high-parent heterosis",
                                                                ifelse(f1pos=="low" & lpvF1_p<=0.05,"below low-parent heterosis","ambiguous")))))))
  d3
}