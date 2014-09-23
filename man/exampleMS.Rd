\name{exampleMS}
\alias{exampleMS}
\docType{data}
\title{SELDI-TOF spectra with known polypeptide positions}   
\description{
The CAMDA 2006 dataset (CAMDA, 2006, http://camda.duke.edu) of All-in-1 Protein Standard II (Ciphergen Cat. C100-0007) was the reference dataset. The MS spectra were measured on Ciphergen NP20 chips. There are seven polypeptides in the sample with the m/z values of 7034, 12 230, 16 951, 29 023,46 671, 66 433 and 147 300.
}
\usage{data(exampleMS)}
\format{
  matrix with 1 objects and 37656 variables:
  \describe{
    \item{\code{exampleMS}}{matrix with 1 rows and 37656 columns}
  }
}
\details{
The data set can be used for background-correct.
}
\examples{
data(exampleMS)
names(exampleMS)
}
\keyword{datasets}
