"cwt" <-
function(ms, scales=1, wavelet='mexh') {
	## Check for the wavelet format
	if (wavelet == 'mexh') {
		psi_xval <- seq(-8, 8, length=1024)
		psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *exp(-psi_xval^2/2)
		#plot(psi_xval, psi)
	} else if (wavelet=='haar') {
    psi_xval <- seq(0,1,length=1024)
    psi <- c(0,rep(1,511),rep(-1,511),0)
	}else if (is.matrix(wavelet)) {
		if (nrow(wavelet) == 2) {
			psi_xval <- wavelet[1,]
			psi <- wavelet[2,]
		} else if (ncol(wavelet) == 2) {
			psi_xval <- wavelet[,1]
			psi <- wavelet[,2]
		} else {
			stop('Unsupported wavelet format!')
		}
	} else {
		stop('Unsupported wavelet!')
	}
		
    oldLen <- length(ms)
	## To increase the computation effeciency of FFT, extend it as the power of 2
	## because of a strange signal length 21577 makes the FFT very slow!
	#ms <- extend.nBase(ms, nLevel=1, base=2)
	ms <- extendNBase(ms, nLevel=NULL, base=2)
	len <- length(ms)
    nbscales <- length(scales)
    wCoefs <- NULL

    psi_xval <- psi_xval - psi_xval[1]
    dxval <- psi_xval[2]
    xmax  <- psi_xval[length(psi_xval)]
    for (i in 1:length(scales)) {
		scale.i <- scales[i]
		f <- rep(0, len)
        j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
        if (length(j) == 1)		j <- c(1, 1)
		lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
		if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
		wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
		## Shift the position with half wavelet width
		wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
		wCoefs <- cbind(wCoefs, wCoefs.i)
    }
	if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
	colnames(wCoefs) <- scales
	wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
	return(wCoefs)
}

