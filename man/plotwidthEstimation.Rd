\name{plotwidthEstimation}
\alias{plotwidthEstimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the estimation of peak width based on the CWT}
\description{
  Plot the the estimation of peak returned by \code{\link{widthEstimationCWT}}
}
\usage{
plotwidthEstimation(x,peakWidth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ raman spectrum }
  \item{peakWidth}{ returned by \code{\link{widthEstimationCWT}} }
}

\value{
Whether the drawing is successful or not.
}


\author{Yizeng Liang ,Zhang Zhimin}

\seealso{\code{\link{widthEstimationCWT}}}
\examples{
x=c(40*dnorm(seq(-5,5,by=0.1),sd=0.2),60*dnorm(seq(-5,5,by=0.1),sd=0.5),30*dnorm(seq(-5,5,by=0.1),sd=0.7),100*rep(0.001,100))
x=x+40*rnorm(length(x))*0.01
b=1:length(x)
yc=30*sin(1:length(x)/100)
yl=30*(1:length(x)/100)
xl=x+yl
xc=x+yc
x=xc  
scales <-seq(1, 30, 1)
wCoefs <- cwt(x, scales=scales, wavelet='mexh')
image(1:nrow(wCoefs), scales, wCoefs, col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
box()
localMax <- getLocalMaximumCWT(wCoefs)
ridgeList <- getRidge(localMax, gapTh=3, skip=2)
majorPeakInfo=identifyMajorPeaks(x, ridgeList, wCoefs,SNR.Th=3,ridgeLength=5)
peakWidth=widthEstimationCWT(x,majorPeakInfo)  
plot(xc,type='l', main="pure signal with curved background and random noise",xlab=expression("Wavenumber / cm"^-1),ylab="Raman Intensity")
plotwidthEstimation(x,peakWidth)   
}
\keyword{plotwidthEstimation}
