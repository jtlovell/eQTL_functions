ctPens<-function(fourPens, nineForms=NULL, cistrtPen=NULL, epiPen=NULL){
  if(is.null(cistrtPen) & is.null(epiPen)){
    cis.trt.pen<-qchisq(0.95,df=2)/2/log(10)
    epi.pen<-qchisq(0.95,df=4)/2/log(10)
  }
  if(is.null(nineForms)){
    forms<-c("y ~ Q1 + trt",
             "y ~ Q1 + Q1*trt + trt",
             "y ~ Q1 + Q2 + trt",
             "y ~ Q1 + Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + Q2*trt + trt")
  }
  pens.all<-c(0,  # plod = 0
              cis.trt.pen, #plod = lod (Q1 *trt)- chi2.2
              fourPens[1], #plod = lod (Q2) - pen1
              fourPens[1] + cis.trt.pen, #plod= lod(Q2) +lod(Q1*trt) - (pen1 + chi2.2)
              fourPens[2] + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen2 + chi2.2)
              fourPens[3] + epi.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen3 + chi2.4) 
              fourPens[3] + epi.pen + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2) 
              fourPens[4] + epi.pen + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2) 
              fourPens[4] + epi.pen + cis.trt.pen + cis.trt.pen) #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2 + chi2.2) -c("y ~ Q1 + trt",
  data.frame(forms=forms, penalties=pens.all, stringsAsFactors=F)
}