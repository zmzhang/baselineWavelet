"localMaximum" <-
function (x, winSize = 5) {
	len <- length(x)
	rNum <- ceiling(len/winSize)
	
	## Transform the vector as a matrix with column length equals winSize
	##		and find the maximum position at each row.
	y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow=winSize)
	y.maxInd <- apply(y, 2, which.max)
	## Only keep the maximum value larger than the boundary values
	selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
	
	## keep the result
	localMax <- rep(0, len)
	localMax[(selInd-1) * winSize + y.maxInd[selInd]] <- 1
	
	## Shift the vector with winSize/2 and do the same operation
	shift <- floor(winSize/2)
	rNum <- ceiling((len + shift)/winSize)	
	y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow=winSize)
	y.maxInd <- apply(y, 2, which.max)
	## Only keep the maximum value larger than the boundary values
	selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
	localMax[(selInd-1) * winSize + y.maxInd[selInd] - shift] <- 1
	
	## Check whether there is some local maxima have in between distance less than winSize
	maxInd <- which(localMax > 0)
	selInd <- which(diff(maxInd) < winSize)
	if (length(selInd) > 0) {
		selMaxInd1 <- maxInd[selInd]
		selMaxInd2 <- maxInd[selInd + 1]
		temp <- x[selMaxInd1] - x[selMaxInd2]
		localMax[selMaxInd1[temp <= 0]] <- 0
		localMax[selMaxInd2[temp > 0]] <- 0
	}
	
	return(localMax)
}

