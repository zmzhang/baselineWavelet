\name{getRidgeLength}
\alias{getRidgeLength}
\title{ Estimate the length of the ridge}
\description{
  Estimate the length of the ridge line, which is composed of local maxima at adjacent CWT scales. The ridge line is cut off at the end point, whose amplitude divided by the maximum ridge amplitude is larger than the cutoff amplitude ratio threshold (0.5 by default).
}
\usage{
getRidgeLength(ridgeList, Th = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ridgeList}{ a list of identified ridges }
  \item{Th}{ the cutoff amplitude ratio (the amplitude divided by the maximum amplitude of the ridge) threshold of the ridge line end.  }
}

\value{
a vector of estimated ridge length
}

\author{Pan Du}

\keyword{methods}
