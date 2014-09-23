"peakDetectionCWT" <-
function(ms, scales=c(1, seq(2,30,2),seq(32, 64, 4)), SNR.Th=3, peakThr=NULL, nearbyPeak=TRUE, peakScaleRange=5, amp.Th=0.05, minNoiseLevel=amp.Th/SNR.Th, ridgeLength=24, tuneIn=FALSE, ...) {

	if (minNoiseLevel > 1)  names(minNoiseLevel) <- 'fixed' 
	## Perform Continuous Wavelet Transform
	wCoefs <- cwt(ms, scales=scales, wavelet='mexh')

	## Attach the raw data as the zero level of decomposition
	wCoefs <- cbind(as.vector(ms), wCoefs)
	colnames(wCoefs) <- c(0, scales)

	##-----------------------------------------
	## Identify the local maximum by using a slide window
	## The size of slide window changes over different levels, with the coarse level have bigger window size
	if (is.null(amp.Th)) {
		amp.Th <- 0 
	} else {
		if (is.null(names(amp.Th))) {
			amp.Th <- max(wCoefs) * amp.Th 
		} else if (names(amp.Th) != 'fixed') {
			amp.Th <- max(wCoefs) * amp.Th
		}
	} 
	localMax <- getLocalMaximumCWT(wCoefs, amp.Th=amp.Th)
	colnames(localMax) <- colnames(wCoefs)

	## In order to fastern the calculation, we can filter some local maxima with small amplitude
	## In this case a baseline estimation was performed.
	if (!is.null(peakThr)) {
		otherPar <- list(...)
		if ('fl' %in% names(otherPar)) {
			filterLength <- otherPar$fl
			otherPar <- otherPar[- which(names(otherPar) == 'fl')]
		} else {
			filterLength <- 1000
		}
		if ('forder' %in% names(otherPar)) {
			fOrder <- otherPar$forder
			otherPar <- otherPar[- which(names(otherPar) == 'forder')]
		} else {
			fOrder <- 2
		}
		## Baseline estimation using Savitzky Golay Filter  
		## this part was added by Steffen Neumann
		sg <- sav.gol(ms, fl=filterLength,forder=fOrder)
		localMax[(ms - sg) < peakThr,] <- 0
	}

	##-----------------------------------------
	## Indentify the ridges from coarse level to more detailed levels
	
#comment by zzm for the simulate data	
# the skip paramter should understand
	ridgeList <- getRidge(localMax, gapTh=3, skip=2)


	##-----------------------------------------
	## Indentify the major peaks and their nearby peaks 
	majorPeakInfo <- identifyMajorPeaks(ms, ridgeList, wCoefs, SNR.Th=SNR.Th, peakScaleRange=peakScaleRange, 
			nearbyPeak=nearbyPeak, minNoiseLevel=minNoiseLevel, ridgeLength=ridgeLength, ...)

	if (tuneIn) {
		refinedPeakInfo <- tuneInPeakInfo(ms, majorPeakInfo)
		return(list(majorPeakInfo=refinedPeakInfo, ridgeList=ridgeList, localMax=localMax, wCoefs=wCoefs[, -1], oldPeakInfo=majorPeakInfo))		
	} else {
		return(list(majorPeakInfo=majorPeakInfo, ridgeList=ridgeList, localMax=localMax, wCoefs=wCoefs[, -1]))		
	}
}

