\name{identifyMajorPeaks}
\alias{identifyMajorPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identify peaks based on the ridges in 2-D CWT coefficient matrix }
\description{
Indentify the peaks based on the ridge list (returned by \code{\link{getRidge}}) in 2-D CWT coefficient matrix and estimated Signal to Noise Ratio (SNR) }
\usage{
identifyMajorPeaks(ms, ridgeList, wCoefs, scales = as.numeric(colnames(wCoefs)), SNR.Th = 3, peakScaleRange = 5, ridgeLength = 32, nearbyPeak = FALSE, nearbyWinSize = 100, winSize.noise = 500, SNR.method = "quantile", minNoiseLevel = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ the mass spectrometry spectrum }  
  \item{ridgeList}{returned by \code{\link{getRidge}}}
  \item{wCoefs}{ 2-D CWT coefficients }
  \item{scales}{ scales of CWT, by default it is the colnames of wCoefs }
  \item{SNR.Th}{ threshold of SNR }
  \item{peakScaleRange}{ the CWT scale range of the peak. }
  \item{ridgeLength}{ the maximum ridge scale of the major peaks. }
  \item{nearbyPeak}{ determine whether to include the small peaks close to large major peaks }
  \item{nearbyWinSize}{ the window size to determine the nearby peaks. Only effective when nearbyPeak is true.}
  \item{winSize.noise}{ the local window size to estimate the noise level. }
  \item{SNR.method}{ method to estimate noise level. Currently, only 95 percentage quantile is supported. }
  \item{minNoiseLevel}{ the minimum noise level used in calculating SNR, i.e., if the estimated noise level is less than "minNoiseLevel", it will use "minNoiseLevel" instead. If the noise level is less than 0.5, it will be treated as the ratio to the maximum amplitude of the spectrum. }
}

\details{
The determination of the peaks is based on three rules:
Rule 1: The maximum ridge scale of the peak should larger than a certain threshold
Rule 2: Based on the scale of the peak (corresponding to the maximum value of the peak ridge) should be within certain range
Rule 3: Based on the peak SNR 
}
\value{
Return a list with following elements:
  \item{peakIndex}{the m/z indexes of the identified peaks}
  \item{peakCenterIndex}{the m/z indexes of peak centers, which correspond to the maximum on the ridge. peakCenterIndex includes all the peaks, not just the identified major peaks.}
  \item{peakCenterValue}{the CWT coefficients (the maximum on the ridge) corresponding to peakCenterIndex}
  \item{peakSNR}{the SNR of the peak, which is the ratio of peakCenterValue  and noise level}
  \item{peakScale}{the estimated scale of the peak, which corresponds to the peakCenerIndex}
  \item{potentialPeakIndex}{the m/z indexes of all potential peaks, which satisfy all requirements of a peak without considering its SNR. Useful, if you want to change to a lower SNR threshold later.}
  \item{allPeakIndex}{the m/z indexes of all the peaks, whose order is the same as peakCenterIndex, peakCenterValue, peakSNR and peakScale.}

All of these return elements have peak names, which are the same as the corresponding peak ridges. see \code{\link{getRidge}} for details.
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M. (2006) Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, Bioinformatics, 22, 2059-2065.
}
\author{Pan Du, Simon Lin}

\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS, scales=scales, wavelet='mexh')
	
	localMax <- getLocalMaximumCWT(wCoefs)
	ridgeList <- getRidge(localMax)
	
	SNR.Th <- 3
	majorPeakInfo <- identifyMajorPeaks(exampleMS, ridgeList, wCoefs, SNR.Th=SNR.Th)
	## Plot the identified peaks
	peakIndex <- majorPeakInfo$peakIndex
}
\keyword{methods}
