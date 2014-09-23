\name{extendNBase}
\alias{extendNBase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extend the row number of a matrix as the exponential of base N}
\description{
   Extend the data as the exponential of base N by increasing row number.
}
\usage{
extendNBase(x, nLevel=1, base=2, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data matrix }
  \item{nLevel}{ the level of DWT decomposition. Basically, it is equivalent to changing the 'base' as base\^nLevel}
  \item{base}{ the base, 2 by default }
  \item{\dots}{ other parameters of used by \code{\link{extendLength}} }
}
\details{
The method 'open' is padding the the matrix with the last row.   
}
\value{
Return a extended matrix
}

\author{ Pan Du }

\seealso{ \code{\link{extendLength}}}
\examples{
	# a = matrix(rnorm(9), 3)
	# extendNBase(a)	## not exposed function
}
\keyword{ methods }% at least one, from doc/KEYWORDS
