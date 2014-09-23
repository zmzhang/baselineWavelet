"extendLength" <-
function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both')) {
	if (is.null(addLength)) stop('Please provide the length to be added!')
	if (!is.matrix(x)) x <- matrix(x, ncol=1)	
	method <- match.arg(method)
	direction <- match.arg(direction)
	
	nR <- nrow(x)
	nR1 <- nR + addLength
	if (direction == 'both') {
		left <- right <- addLength
	} else if (direction == 'right') {
		left <- 0
		right <- addLength
	} else if (direction == 'left') {
		left <- addLength
		right <- 0
	}
	
	if (right > 0) {
		x <- switch(method,
			reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
			open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
			circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
	}

	if (left > 0) {
		x <- switch(method,
			reflection =rbind(x[addLength:1, , drop=FALSE], x),
			open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
			circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
	}
	if (ncol(x) == 1)  x <- as.vector(x)
	
	return(x)
}

