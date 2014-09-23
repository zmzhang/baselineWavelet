"widthEstimationCWT" <-
function(x,majorPeakInfo) {

  wCoefs_haar <- cwt(x, 1:max(majorPeakInfo$peakScale), wavelet='haar')
  peakIndex <- majorPeakInfo$peakIndex
  peakScale <- majorPeakInfo$peakScale[findInterval(majorPeakInfo$peakIndex,majorPeakInfo$allPeakIndex)]
  
  peakWidth <- list()
  peakWidth[["peakIndex"]]= majorPeakInfo$peakIndex
  
  for(i in 1:length(peakIndex)){
    peakIndex.i=peakIndex[i]
    peakScale.i=peakScale[i]
    wCoefs_haar.i=wCoefs_haar[,peakScale.i]
    wCoefs_haar.i.abs=abs(wCoefs_haar.i)
    localmax=localMaximum(-wCoefs_haar.i.abs,winSize=5)      
#    localmax=localmax & abs(wCoefs_haar.i)<(mean(wCoefs_haar.i.abs[localmax==1])+0.5*sd(wCoefs_haar.i.abs[localmax==1]))
    localmax=as.numeric(localmax)
    localmax[peakIndex]=0
    localmax[(peakIndex.i-peakScale.i/2+1):(peakIndex.i+peakScale.i/2-1)]=0
        
    Lef =0
    Rig =0
    
    peakScale.i.3=2*peakScale.i
    
     
    if(i==1){
       maxIndexL=max(c((peakIndex.i-peakScale.i.3),1))
    }else{
       maxIndexL=max(c((peakIndex.i-peakScale.i.3),peakIndex[i-1]))
    }
    
    if(i==length(peakIndex)){
       minIndexR=min(c((peakIndex.i+peakScale.i.3),length(localmax)))
    } else{
       minIndexR=min(c((peakIndex.i+peakScale.i.3),peakIndex[i+1]))
    }
    ignoreL=1:maxIndexL
    ignoreR=minIndexR:length(localmax)        
    localmax[c(ignoreL,ignoreR)]=0
    localmax[c(peakIndex.i,(peakIndex.i-(peakScale.i/2)):(peakIndex.i+(peakScale.i/2)))]=0
    bi=which(localmax==1)
        
    biLeft=bi[bi<peakIndex.i]
    useL= maxIndexL:peakIndex.i
    minIndexLeft=useL[which(min(x[useL])==x[useL])]
    
    if(length(biLeft)==0){
      Lef=minIndexLeft
    }else{
      minbaselineIndexLeft=biLeft[which(min(x[biLeft])==x[biLeft])]
      if(minIndexLeft>=(peakIndex.i-peakScale.i/2+1)){
        Lef=minbaselineIndexLeft
      }else{
        Lef=max(c(minIndexLeft,minbaselineIndexLeft))
      }   
    }
    
    biRight=bi[bi>peakIndex.i]
    useR= peakIndex.i:minIndexR
    minIndexRight=useR[which(min(x[useR])==x[useR])]
      
    if(length(biRight)==0){
      Rig= minIndexRight
    }else{
      minbaselineIndexRight=biRight[which(min(x[biRight])==x[biRight])] 
      if(minIndexRight<=(peakIndex.i+peakScale.i/2-1)){
        Rig=minbaselineIndexRight
      }else{
        Rig=min(c(minIndexRight,minbaselineIndexRight))
      }      
    }
        
    peakWidth[[paste(peakIndex.i)]]=Lef:Rig
    
  }
      
      
      
#      lmindex=1:length(localmax)
#      plot(wCoefs_haar.i,type='l',ylim=c(min(c(x,wCoefs_haar.i)),max(c(x,wCoefs_haar.i))))
#      lines(x,col='green')
#      points(peakIndex,x[peakIndex],col='blue')   
#      points(bi,x[bi])
         
#      points(lmindex[c(Lef,Rig)],x[c(Lef,Rig)],cex='15',col='red')
#      points(lmindex[localmax==1],wCoefs_haar.i[localmax==1],cex='15',col='red')
#      points(lmindex[localmax==1],x[localmax==1],cex='15',col='red')
      
      


  return(peakWidth)		
}

