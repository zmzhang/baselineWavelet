\name{localMaximum}
\alias{localMaximum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify local maximum within a slide window.}
\description{
 Find local maximum by transform the vector as matrix, then get the the maximum of each column. This operation is performed twice with vecctor shifted half of the winSize.}
\usage{
localMaximum(x, winSize = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector represents a signal profile}
  \item{winSize}{ the slide window size, 5 by default.}
}
\details{
Instead of find the local maximum by a slide window, which slide all possible positions, we find local maximum by transform the vector as matrix, then get the the maximum of each column. This operation is performed twice with vecctor shifted half of the winSize. The main purpose of this is to increase the efficiency of the algorithm.}
\value{
 Return a vector with the same length of the input x. The position of local maximum is set as 1, 0 else where.
}

\author{ Pan Du }
\seealso{ \code{\link{getLocalMaximumCWT}}}
\examples{
x <- rnorm(200)
lmax <- localMaximum(x, 5)
maxInd <- which(lmax > 0)
plot(x, type='l')
points(maxInd, x[maxInd], col='red')
}
\keyword{methods}
