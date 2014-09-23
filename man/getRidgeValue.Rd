\name{getRidgeValue}
\alias{getRidgeValue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the CWT coefficient values corresponding to the peak ridge}
\description{
Get the CWT coefficient values corresponding to the peak ridge
}
\usage{
getRidgeValue(ridgeList, wCoefs, skip = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ridgeList}{ a list of ridge lines }
  \item{wCoefs}{ 2-D CWT coefficients}
  \item{skip}{ the CWT scale level to be skipped, by default the 0 scale level (raw spectrum) is skipped. }
}

\value{
	A list of ridge values corresponding to the input ridgeList.
}

\author{Pan Du}


\keyword{methods}
