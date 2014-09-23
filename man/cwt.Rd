\name{cwt}
\alias{cwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Continuous Wavelet Transform (CWT) }
\description{
  CWT(Continuous Wavelet Transform) with Mexican Hat wavelet (by default) to match the peaks in Mass Spectrometry spectrum
}
\usage{
cwt(ms, scales = 1, wavelet = "mexh")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{ Mass Spectrometry spectrum (a vector of MS intensities)}
  \item{scales}{ a vector represents the scales at which to perform CWT. }
  \item{wavelet}{The wavelet base, Mexican Hat by default. User can provide wavelet Psi(x) as a form of two row matrix. The first 
	row is the x value, and the second row is Psi(x) corresponding to x.  }
}

\value{
The return is the 2-D CWT coefficient matrix, with column names as the scale. Each column is the CWT coefficients at that scale.
}

\author{Pan Du, Simon Lin}

\examples{
	data(exampleMS)
	scales <- seq(1, 64, 3)
	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')

	## Plot the 2-D CWT coefficients as image (It may take a while!)
	xTickInterval <- 1000
	image(5000:11000, scales, wCoefs, col=terrain.colors(256), axes=FALSE, xlab='m/z index', ylab='CWT coefficient scale', main='CWT coefficients')
	axis(1, at=seq(5000, 11000, by=xTickInterval))
	axis(2, at=c(1, seq(10, 64, by=10)))
	box()
}
\keyword{ methods }
