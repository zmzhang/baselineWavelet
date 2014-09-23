"identifyMajorPeaks" <-
function(ms, ridgeList, wCoefs, scales=as.numeric(colnames(wCoefs)), SNR.Th=1, peakScaleRange=5, 
		ridgeLength=10, nearbyPeak=FALSE, nearbyWinSize=100, winSize.noise=500, SNR.method='quantile', minNoiseLevel=0.001) {

	if (is.null(scales)) {
		scales <- 1:ncol(wCoefs)
		colnames(wCoefs) <- scales
	} else if (is.character(scales)) {
		scales <- as.numeric(scales)
	}
	if (ridgeLength > max(scales))  ridgeLength <- max(scales)
		
	if (length(peakScaleRange) == 1) {
		peakScaleRange <- scales[scales >= peakScaleRange]
	} else {
		peakScaleRange <- scales[scales >= peakScaleRange[1] & scales <= peakScaleRange[2]]
	}

	## Limit the minNoiseLevel to avoid the case of very low noise level, e.g., smoothed spectrum
	if (minNoiseLevel >= 1)  names(minNoiseLevel) <- 'fixed'
	if (is.null(minNoiseLevel)) {
		minNoiseLevel <- 0
	} else {	# By default the threshold is the ratio of the maximum coefficient
		if (is.null(names(minNoiseLevel))) {
			minNoiseLevel <- max(wCoefs) * minNoiseLevel 
		} else if (names(minNoiseLevel) != 'fixed') {		
			minNoiseLevel <- max(wCoefs) * minNoiseLevel
		}
	}
		
	## Get the peak values
	#mzInd <- as.numeric(names(ridgeList))
	ridgeLen <- sapply(ridgeList, length)
	ridgeName <- names(ridgeList)
	ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName, '_'))), nrow=2)
	ridgeLevel <- ridgeInfo[1,]
	#mzInd <- sapply(ridgeList, function(x) x[1])
	notnull <- sapply(ridgeList, function(x) {!is.null(x[1])})  # fixed by Steffen Neumann
	mzInd <- sapply(ridgeList[notnull], function(x) {x[1]})     # fixed by Steffen Neumann
	#mzInd <- ridgeInfo[2,]

	## Reorder them by m/z index
	ord <- order(mzInd)	
	ridgeName <- ridgeName[ord]
	ridgeLen <- ridgeLen[ord]
	ridgeLevel <- ridgeLevel[ord]
	ridgeList <- ridgeList[ord]
	mzInd <- mzInd[ord]

	peakScale <- NULL
	peakCenterInd <- NULL
	peakValue <- NULL
	#ridgeValue <- NULL
	## Get the ridge values within the provided peakScaleRange
	for (i in 1:length(ridgeList)) {
		ridge.i <- ridgeList[[i]]
		level.i <- ridgeLevel[i]
		levels.i <- level.i:(level.i + ridgeLen[i] - 1)
		scales.i <- scales[levels.i]
		# Only keep the scales within the peakScaleRange
		selInd.i <- which(scales.i %in% peakScaleRange)
		if (length(selInd.i) == 0) {
			peakScale <- c(peakScale, scales.i[1])
			peakCenterInd <- c(peakCenterInd, ridge.i[1])
			peakValue <- c(peakValue, 0)
			next
		} 
		
		levels.i <- levels.i[selInd.i]
		scales.i <- scales.i[selInd.i]
		ridge.i <- ridge.i[selInd.i]
		if (scales.i[1] == 0) {
			ind.i <- cbind(ridge.i[-1], levels.i[-1])
		} else {
			ind.i <- cbind(ridge.i, levels.i)
		}
		ridgeValue.i <- wCoefs[ind.i]
		maxInd.i <- which.max(ridgeValue.i)
		peakScale <- c(peakScale, scales.i[maxInd.i])
		peakCenterInd <- c(peakCenterInd, ridge.i[maxInd.i])
		peakValue <- c(peakValue, ridgeValue.i[maxInd.i])
		#ridgeValue <- c(ridgeValue, list(ridgeValue))
	}
	#names(ridgeValue) <- names(ridgeList)		
	#ridgeLen <- get.ridgeLength(ridgeValue, Th=0.5)
	
	## Compute SNR of each peak
	noise <- abs(wCoefs[,'1'])
	peakSNR <- NULL	
	nMz <- nrow(wCoefs)		# The length of ms signal
	
	for (k in 1:length(ridgeList)) {
		ind.k <- mzInd[k]
		start.k <- ifelse(ind.k - winSize.noise < 1, 1, ind.k - winSize.noise)
		end.k <- ifelse(ind.k + winSize.noise > nMz, nMz, ind.k + winSize.noise)
		ms.int <- ms[start.k:end.k] ## m/z intensity values in ind.k + /- winSize.noise (Added by Steffen Neumann)
		noiseLevel.k <- switch(SNR.method,
			quantile = quantile(noise[start.k:end.k], probs=0.95), 
			sd = sd(noise[start.k:end.k]),
			mad = mad(noise[start.k:end.k], center=0),
			data.mean = mean(ms.int),		# (data.mean and data.mean.quant were added by Steffen Neumann)
			data.mean.quant = mean(ms.int[ ms.int < quantile(ms.int,probs=.95) ])  )
		## Limit the minNoiseLevel to avoid the case of very low noise level, e.g., smoothed spectrum
		if (noiseLevel.k < minNoiseLevel) noiseLevel.k <- minNoiseLevel
		peakSNR <- c(peakSNR, peakValue[k]/noiseLevel.k)
	}

	## Rule 1: ridge length should larger than a certain threshold
	#selInd1 <- (scales[ridgeLen] >= ridgeLength)
	selInd1 <- (scales[ridgeLevel + ridgeLen - 1] >= ridgeLength)
	
	## In the case of nearbyPeak mode, it will include the nearby peaks within a certain range
	if (0) {
		selInd1 <- which(selInd1)
		index <- 1:length(mzInd)
		nearbyWinSize <- 150
		tempInd <- NULL
		for (ind.i in selInd1) {
			tempInd <- c(tempInd, index[mzInd >= mzInd[ind.i] - nearbyWinSize & mzInd <= mzInd[ind.i] + nearbyWinSize])
		}
		selInd1 <- (index %in% tempInd)
	}
	
	## Rule 2: Based on the peak SNR
	selInd2 <- (peakSNR > SNR.Th)
	
	## Because of the boundary effects,
	## remove the peaks (half of the nearbyWinSize) at both ends of the signal profile if exists
	selInd3 <- !(mzInd %in% c(1:(nearbyWinSize/2), (nrow(wCoefs) - (nearbyWinSize/2) + 1):nrow(wCoefs)))

	## combine SNR and peak length rule and other rules
	selInd <- (selInd1 & selInd2 & selInd3)
#  selInd <- (selInd1 & selInd2)
#selInd <- (selInd1)
	
	names(peakSNR) <- names(peakScale) <- names(peakCenterInd) <- names(peakValue) <- names(mzInd) <- ridgeName

	return(list(peakIndex=mzInd[selInd], peakValue=peakValue, peakCenterIndex=peakCenterInd, peakSNR=peakSNR, peakScale=peakScale, potentialPeakIndex=mzInd[selInd1 & selInd3], allPeakIndex=mzInd))
}

