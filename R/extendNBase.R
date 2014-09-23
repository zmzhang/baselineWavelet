"extendNBase" <-
function(x, nLevel=1, base=2, ...) {
	if (!is.matrix(x)) x <- matrix(x, ncol=1)	
	
	nR <- nrow(x)
	if (is.null(nLevel)) {
		nR1 <- nextn(nR, base)		
	} else {
		nR1 <- ceiling(nR / base^nLevel) * base^nLevel		
	}
	if (nR != nR1) {
		x <- extendLength(x, addLength=nR1-nR, ...)
	}

	return(x)
}

