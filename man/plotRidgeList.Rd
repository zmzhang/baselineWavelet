\name{plotRidgeList}
\alias{plotRidgeList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the ridge list}
\description{
  Plot the ridge list returned by \code{\link{getRidge}}
}
\usage{
plotRidgeList(ridgeList, wCoefs = NULL, range = NULL, colorMap = "RYB", main = NULL, pch = ".", cex = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ridgeList}{ returned by \code{\link{getRidge}} }
  \item{wCoefs}{ 2-D CWT coefficients }
  \item{range}{ plot range of m/z index }
  \item{colorMap}{ colorMap to plot the points of local maximum }
  \item{main}{parameter of \code{\link{plot}}}
  \item{pch}{ parameter of \code{\link{plot}}}
  \item{cex}{ parameter of \code{\link{plot}} }
  \item{\dots}{  other parameters of \code{\link{points}}  }
}


\author{ Pan Du}

\seealso{\code{\link{getRidge}}}
\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')
	
	localMax <- getLocalMaximumCWT(wCoefs)
	ridgeList <- getRidge(localMax)
	plotRidgeList(ridgeList)
}
\keyword{hplot}
