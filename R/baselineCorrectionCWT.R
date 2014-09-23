"baselineCorrectionCWT" <-
function(x,peakWidth,threshold=0.5,lambda=100,differences=1) {

  # fitting an rough background with proper value of lambda
  x.w=rep(1,length(x))
  peakIndex=peakWidth$peakIndex
  for (i in 1:length(peakIndex)){
    x.w[peakWidth[[paste(peakIndex[i])]]]=0
  }
  backgr = WhittakerSmooth(x,x.w,lambda,differences)
  
  #assigning the corresponding value of the spectrum to the part of the rough background,which is larger than the spectrum. and fitting another background without value lager than orignal spectrum based on the new rough background
  
  backgr.final=backgr
  backgr.final[x<=backgr]=x[x<=backgr]
  backgr.final= WhittakerSmooth(backgr.final,rep(1,length(backgr.final)),differences)
  x= backgr.final
  
  #refine the new rough background
  
  background = NULL
  peakIndex=peakWidth$peakIndex
  signal_baseline= x[1:(peakWidth[[paste(peakIndex[1])]][1])-1]
  signal_baseline.w= rep(1,length(signal_baseline))
  baseline_baseline =   WhittakerSmooth(signal_baseline,signal_baseline.w,differences)
  background = c(background,baseline_baseline)
  for (i in 1:(length(peakIndex)-1)){

    peakWidth.1=peakWidth[[paste(peakIndex[i])]]
    peakWidth.2=peakWidth[[paste(peakIndex[i+1])]]

    if(length(intersect(peakWidth.1, peakWidth.2))==0)
    {
      signal_peak=x[peakWidth.1]
      signal_peak.w=c(1,rep(0,(length(peakWidth.1)-2)),1)
      if((abs(signal_peak[length(signal_peak)]-signal_peak[1])/mean(signal_peak-min(signal_peak)))>=threshold){
        if(signal_peak[length(signal_peak)]>signal_peak[1]){
          signal_peak.w=as.numeric(!(signal_peak>signal_peak[length(signal_peak)]))
        }else{
          signal_peak.w=as.numeric(!(signal_peak>signal_peak[1]))
        }
        signal_peak.w[1]=1
        signal_peak.w[length(signal_peak.w)]=1
      }
      peak_baseline= WhittakerSmooth(signal_peak,signal_peak.w,differences)
      background = c(background,peak_baseline)

      if(peakWidth.2[1]-peakWidth.1[length(peakWidth.1)]==1){
        # two peak just together
      }else{
        signal_baseline=x[(peakWidth.1[length(peakWidth.1)]+1):(peakWidth.2[1]-1)]
        if(length(signal_baseline)<=3){
          baseline_baseline=signal_baseline
        }else{
          signal_baseline.w=c(1,rep(1,(length(signal_baseline)-2)),1)
	    baseline_baseline= WhittakerSmooth(signal_baseline,signal_baseline.w,differences)
        }
        background = c(background,baseline_baseline)
      }

    }else{
      if(peakWidth.1[1]>peakWidth.2[1]){
         # the last peak contain the preivous peak
         peakWidth[[paste(peakIndex[i+1])]]=(peakWidth.1[1]):(peakWidth.2[length(peakWidth.2)])
      }else{
        if(length(peakWidth.1)>=length(peakWidth.2)){
          peakWidth.11=setdiff(peakWidth.1,intersect(peakWidth.1, peakWidth.2))
        }else{
          peakWidth.11=peakWidth.1
          peakWidth[[paste(peakIndex[i+1])]]=setdiff(peakWidth.2,intersect(peakWidth.1, peakWidth.2))
        }

        if(length(peakWidth.11)<=2){
          peakWidth[[paste(peakIndex[i+1])]]=(peakWidth.1[1]):(peakWidth.2[length(peakWidth.2)])
        } else{
          signal_peak=x[peakWidth.11]
          signal_peak.w=c(1,rep(0,(length(peakWidth.11)-2)),1)
          if((abs(signal_peak[length(signal_peak)]-signal_peak[1])/mean(signal_peak-min(signal_peak)))>=threshold){
            if(signal_peak[length(signal_peak)]>signal_peak[1]){
              signal_peak.w=as.numeric(!(signal_peak>signal_peak[length(signal_peak)]))
            }else{
              signal_peak.w=as.numeric(!(signal_peak>signal_peak[1]))
            }
            signal_peak.w[1]=1
            signal_peak.w[length(signal_peak.w)]=1
          }
          peak_baseline= WhittakerSmooth(signal_peak,signal_peak.w,differences)
          background = c(background,peak_baseline)
        }

      }
    }
  }

  peakWidth.end= peakWidth[[paste(peakIndex[length(peakIndex)])]]
  signal_peak=x[peakWidth.end]
  signal_peak.w=c(1,rep(0,(length(peakWidth.end)-2)),1)
  if((abs(signal_peak[length(signal_peak)]-signal_peak[1])/mean(signal_peak-min(signal_peak)))>=threshold){
    if(signal_peak[length(signal_peak)]>signal_peak[1]){
      signal_peak.w=as.numeric(!(signal_peak>signal_peak[length(signal_peak)]))
    }else{
      signal_peak.w=as.numeric(!(signal_peak>signal_peak[1]))
    }
  }
  peak_baseline= WhittakerSmooth(signal_peak,signal_peak.w,1)

  signal_baseline=x[(peakWidth.end[length(peakWidth.end)]+1):length(x)]
  signal_baseline.w=rep(1,length(signal_baseline))
  baseline_baseline= WhittakerSmooth(signal_baseline,signal_baseline.w,1)
  background = c(background,peak_baseline)
  background = c(background,baseline_baseline)

  return(background)
}