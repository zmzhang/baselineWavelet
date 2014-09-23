"getLocalMaximumCWT" <-
function(wCoefs, minWinSize=5, amp.Th=0) {

	localMax <- NULL
	scales <- as.numeric(colnames(wCoefs))
	
	for (i in 1:length(scales)) {
		scale.i <- scales[i]
		winSize.i <- scale.i * 2 + 1
		if (winSize.i < minWinSize) {
			winSize.i <- minWinSize
		} 
		temp <- localMaximum(wCoefs[,i], winSize.i)
		localMax <- cbind(localMax, temp)
	}
	# Set the values less than peak threshold as 0
	localMax[wCoefs < amp.Th] <- 0
	colnames(localMax) <- colnames(wCoefs)
	rownames(localMax) <- rownames(wCoefs)
	return(localMax)
}