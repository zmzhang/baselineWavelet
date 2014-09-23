\name{extendLength}
\alias{extendLength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend the length of a signal or matrix}
\description{
Extend the length of a signal or matrix by row
}
\usage{
extendLength(x, addLength = NULL, method = c("reflection", "open", "circular"), direction = c("right", "left", "both"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector or matrix with column with each column as a signal }
  \item{addLength}{ the length to be extended }
  \item{method}{ three methods available, c("reflection", "open", "circular"). By default, it is "reflection". }
  \item{direction}{ three options available: c("right", "left", "both") }
}

\value{
return the extended vector or matrix. 
}

\author{ Pan Du }

\seealso{ \code{\link{extendNBase}}}
\examples{
	# a = matrix(rnorm(9), 3)
	# extendLength(a, 3, direction='right') 	## not exposed function
}
\keyword{methods}
