\name{raman}
\alias{raman}
\docType{data}
\title{Raman spectra of Prednisone Acetate Tablets}   
\description{
Prednisone Acetate Tablets (PAT) were collected from ten different pharmaceutical factories, which were analyzed by BWTEK BTC112 Raman spectrometer. Measurements were carried out using laser of 784.7nm wavelength for excitation. A thermoelectric cooled (TEC) linear charge-coupled device (CCD) arrays with 2048 elements was used for signal detection, which provided high stability and low dark counts. All the tablets are recorded in the same integration times of 5000ms to obtain comparable spectra. Then these Raman spectra are corrected to test the proposed background-correction algorithm. 
}
\usage{data(raman)}
\format{
  two matrix with 10 objects and 1715 variables:
  \describe{
    \item{\code{m}}{matrix with 10 rows and 1715 columns}
    \item{\code{xa}}{matrix with 1 rows and 1715 columns}
  }
}
\details{
The data set can be used for background-correct.
}
\examples{
data(raman)
}
\keyword{datasets}
