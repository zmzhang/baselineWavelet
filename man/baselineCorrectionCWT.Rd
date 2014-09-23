\name{baselineCorrectionCWT}
\alias{baselineCorrectionCWT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{background fitting via penalized least squares algorithm}
\description{
   intelligent background-correction algorithm is developed, which simulates manual background-correction procedure of an expert intelligently. It basically consists of works of three aspects: 1) the accurate peak position detection of Raman spectrum by continuous wavelet transform (CWT) with the Mexican Hat wavelet as the mother wavelet; 2) peak widths estimation by enhanced signal-to-noise ratio (SNR) derivative calculation based on CWT but with the Haar wavelet as the mother wavelet; and 3) background fitting using penalized least squares with binary masks.
}
\usage{
baselineCorrectionCWT(x,peakWidth,threshold=0.5,lambda=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Raman spectrum }
  \item{peakWidth}{ returned by \code{\link{widthEstimationCWT}} }
  \item{threshold}{ the user define peak shape threshold }
  \item{lambda}{lambda is an adjustable parameter, it can be adjusted by user. The larger lambda is,  the smoother fitted background will be }
  \item{differences}{ an integer indicating the order of the difference of penalties of Whittaker Smoother method, see \code{\link{WhittakerSmooth}}}
}

\value{
  the final background
}

\author{Yizeng Liang ,Zhang Zhimin}

\seealso{\code{\link{WhittakerSmooth}}}
\examples{
data(raman)
x=m[7,]
scales <-seq(1, 70, 1)
wCoefs <- cwt(x, scales=scales, wavelet='mexh')
image(1:nrow(wCoefs), scales, wCoefs, col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
box()
localMax <- getLocalMaximumCWT(wCoefs)
ridgeList <- getRidge(localMax, gapTh=3, skip=2)
plotRidgeList(ridgeList)
majorPeakInfo = identifyMajorPeaks(x, ridgeList, wCoefs, SNR.Th=1,ridgeLength=5)
peakWidth=widthEstimationCWT(x,majorPeakInfo)                                                
backgr = baselineCorrectionCWT(x,peakWidth,lambda=1000,differences=1)
corrected=x-backgr 
plot(xa,x,type='l',ylim=c(min(c(x,corrected)),max(c(x,corrected))),main="The background-correction result of Raman Spectra",xlab=expression("Wavenumber / cm"^-1),ylab="Raman Intensity/Arbitr. Units")
points(xa[majorPeakInfo$peakIndex],x[majorPeakInfo$peakIndex])
lines(xa,backgr,lty=5)
lines(xa,corrected)
}
\keyword{widthEstimationCWT}
