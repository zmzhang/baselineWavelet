"plotRidgeList" <-
function(ridgeList, wCoefs=NULL, range=NULL, colorMap='RYB', main=NULL, pch='.', cex=3, ...) {
	if (colorMap == 'RYB') {
		rgb.palette <- colorRampPalette(c("red", "orange", "blue"),
                                space = "rgb")
		colorMap <- rgb.palette(255)
		#colorMap <- heat.colors(255)
	}
	ridgeLen <- sapply(ridgeList, length)
	ridgeName <- names(ridgeList)
	ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName, '_'))), nrow=2)
	ridgeLevel <- ridgeInfo[1,]
	mzInd <- ridgeInfo[2,]
	
	scales <- attr(ridgeList, 'scales')
	if (is.null(scales)) scales <- 0:max(ridgeLen)
	## Check range parameter
	if (!is.null(range)) {
		range <- sort(range)
		if (length(range) == 1) range <- c(range, max(mzInd))
		if (range[1] < 1) range[1] <- 1
		if (range[2] > max(mzInd)) range[2] <- max(mzInd)
	} else {
		range <- c(1, max(mzInd))
	}
	
	selInd <- which(mzInd >= range[1] & mzInd <= range[2])
	selMzInd <- mzInd[selInd]
	selRidgeList <- ridgeList[selInd]
	selRidgeLen <- ridgeLen[selInd]
	selRidgeLevel <- ridgeLevel[selInd]

	plot(seq(range[1], range[2], length=length(scales)), scales, type='n', xlab='index', ylab='CWT coefficient scale', main=main)
	if (!is.null(wCoefs)) {
		if (length(scales) != ncol(wCoefs))  stop('The length of "scales" does not match with "wCoefs"!')
		#if (!is.null(skip)) {
		#	wCoefs <- wCoefs[, -skip]
		#}
		range <- round(range)
		maxSelV <- max(wCoefs[range[1]:range[2],])
	}

	for (i in 1:length(selMzInd)) {
		ridge.i <- selRidgeList[[i]]
		level.i <- selRidgeLevel[i]
		levels.i <- level.i:(level.i + selRidgeLen[i] - 1)
		scales.i <- scales[levels.i]

		if (is.null(wCoefs)) {
			col <- 'blue'
		} else {
			coefInd.i <- cbind(ridge.i, levels.i)
			ridgeValue.i <- wCoefs[coefInd.i]
			col <- colorMap[1+log(1+abs(ridgeValue.i))*255/log(1+maxSelV)]	
		}
		points(ridge.i, scales.i, pch=pch, cex=cex, col=col, ...) 
	}
}

