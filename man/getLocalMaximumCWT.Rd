\name{getLocalMaximumCWT}
\alias{getLocalMaximumCWT}
\title{Identify the local maximum of each column in 2-D CWT coefficients matrix}
\description{
Identify the local maximum of each column in 2-D CWT coefficients matrix by using a slide window.
The size of slide window linearly changes from the coarse scale (bigger window size) to detail scale.
The scale of CWT increases with the column index.
}
\usage{
getLocalMaximumCWT(wCoefs, minWinSize= 5, amp.Th = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wCoefs}{2-D CWT coefficients, each column corresponding to CWT coefficient at one scale. The column name is the scale.}
  \item{minWinSize}{The minimum slide window size used.}
  \item{amp.Th}{ The minimum peak amplitude. }
}

\value{
return a matrix with same dimension as CWT coefficient matrix, wCoefs. The local maxima are marked as 1, others are 0.
}

\author{ Pan Du }

\seealso{ \code{\link{localMaximum}} }
\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')
	
	localMax <- getLocalMaximumCWT(wCoefs)
}
\keyword{ methods }


