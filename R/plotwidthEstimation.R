"plotwidthEstimation" <-
function(x,peakWidth) {
  lmindex=1:length(x)
  peakIndex=peakWidth$peakIndex
  LR = NULL
#  plot(x,type='l')                       
  points(peakIndex,x[peakIndex])
  for (i in 1:length(peakIndex)){
    peakWidth.i=peakWidth[[paste(peakIndex[i])]]
    LR=c(LR,peakWidth.i[c(1,length(peakWidth.i))])
    text(peakWidth.i[c(1,length(peakWidth.i))],x[peakWidth.i[c(1,length(peakWidth.i))]]-100,paste(i),cex = .8)    
  }
  points(LR,x[LR])
  
  return ("successful")
}