\name{getRidge}
\alias{getRidge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify ridges based on the local maximum matrix}
\description{
Identify ridges by connecting the local maximum of 2-D CWT coefficients from the coarse scale to detail scale. The local maximum matrix is returned  from \code{\link{getLocalMaximumCWT}}
}
\usage{
getRidge(localMax, iInit = ncol(localMax), step = -1, iFinal = 1, minWinSize= 5, gapTh = 3, skip = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{localMax}{ The local maximum matrix is returned  from \code{\link{getLocalMaximumCWT}} with 1 represents maximum, others are 0. }
  \item{iInit}{ The start column to search ridge. By default, it starts from the coarsest scale level.}
  \item{step}{ Search step. -1 by default, which means searching from coarse scale to detail scale column by column. }
  \item{iFinal}{ The final column index of search ridge. }
  \item{minWinSize}{ The minimum slide window size used. }
  \item{gapTh}{ The gap allowed during searching for ridge. 3 by default. }
  \item{skip}{ The column to be skipped during search. }
}

\value{
 Return a list of ridge. As some ridges may end at the scale larger than 1, in order to keep the uniqueness of the ridge names, we combined the smallest scale of the ridge and m/z index of the peak at that scale together to name the ridges. For example the ridge name "1\_653" means the peak ridge ends at the CWT scale 1 with m/z index 653 at scale 1. 
}
\references{ 
Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics, 22, 2059-2065.
}
\author{ Pan Du, Simon Lin }
\seealso{ \code{\link{getLocalMaximumCWT}}, \code{\link{identifyMajorPeaks}}}

\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')
	
	localMax <- getLocalMaximumCWT(wCoefs)
	ridgeList <- getRidge(localMax)
	plotRidgeList(ridgeList)
}
\keyword{methods}
