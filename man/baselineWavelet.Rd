\name{baselineWavelet-package}
\alias{baselineWavelet-package}
\alias{baselineWavelet}
\docType{package}
\title{
Intelligent background-correction algorithm for highly fluorescent sample in Raman spectroscopy
}
\description{
an intelligent background-correction algorithm is developed, which simulates manual background-correction procedure of an expert intelligently. It basically consists of works of three aspects: 1) the accurate peak position detection of Raman spectrum by continuous wavelet transform (CWT) with the Mexican Hat wavelet as the mother wavelet; 2) peak widths estimation by enhanced signal-to-noise ratio (SNR) derivative calculation based on CWT but with the Haar wavelet as the mother wavelet; and 3) background fitting using penalized least squares with binary masks.\cr

What's new in newer version:  \cr
1: from version 2.0 to 3.0: Rewirte the WhittakerSmooth function, don't use the cholskey decomposition any more.\cr
2: from version 1.0 - 2.0: Two functions, say baselineCorrectionCWT() and WhittakerSmoother(), in the baselineWavelet package were modified (add a parameter) so that one could easily perform first, second or even higher differences penalties by adjusting the parameter for the purpose. \cr
}
\details{
\tabular{ll}{
Package: \tab baselineWavelet\cr
Type: \tab Package\cr
Version: \tab 3.0.0\cr
Date: \tab 2009-09-03\cr
License: \tab GPL (>= 2)\cr
}
The package can be used to not only an intelligent and valid tool for background correction in Raman spectroscopy, but also a useful tool in other spectral signals and chromatogram.
}
\author{ yizeng liang<yizeng_liang@263.net>, zhimin zhang <zhangzhimin.csu@gmail.com>} 
\references{
[1]	D. L. A. de Faria, M. A. de Souza, Journal of Raman Spectroscopy 1999, 30, 169.\cr
[2]	R. P. Millen, M. L. A. Temperini, D. L. A. de Faria, D. N. Batchelder, Journal of Raman Spectroscopy 1999, 30, 1027.\cr
[3]	M. J. Pelletier, R. Altkorn, Applied Spectroscopy 2000, 54, 1837.\cr
[4]	A. O'Grady, A. C. Dennis, D. Denvir, J. J. McGarvey, S. E. J. Bell, Analytical Chemistry 2001, 73, 2058.\cr
[5]	D. Rohleder, W. Kiefer, W. Petrich, Analyst 2004, 129, 906.\cr
[6]	D. M. Zhang, M. F. Mrozek, Y. Xie, D. Ben-Amotz, Applied Spectroscopy 2004, 58, 929.\cr
[7]	A. M. Macdonald, P. Wyeth, Journal of Raman Spectroscopy 2006, 37, 830.\cr
[8]	Y. Oshima, Y. Komachi, C. Furihata, H. Tashiro, H. Sato, Applied Spectroscopy 2006, 60, 964.\cr
[9]	I. Osticioli, A. Zoppi, E. M. Castellucci, Journal of Raman Spectroscopy 2006, 37, 974.\cr
[10]	I. Osticioli, A. Zoppi, E. M. Castellucci, Applied Spectroscopy 2007, 61, 839.\cr
[11]	T. Hasegawa, J. Nishijo, J. Umemura, Chemical Physics Letters 2000, 317, 642.\cr
[12]	D. M. Zhang, D. Ben-Amotz, Applied Spectroscopy 2000, 54, 1379.\cr
[13]	A. Jirasek, G. Schulze, M. M. L. Yu, W. Blades, R. F. B. Turner, Applied Spectroscopy 2004, 58, 1488.\cr
[14]	V. Mazet, C. Carteret, D. Brie, M. Idier, B. Humbert, Chemometrics and Intelligent Laboratory Systems 2005, 76, 121.\cr
[15]	P. M. Ramos, I. Ruisanchez, Journal of Raman Spectroscopy 2005, 36, 848.\cr
[16]	C. Camerlingo, F. Zenone, G. M. Gaeta, R. Riccio, M. Lepore, Measurement Science & Technology 2006, 17, 298.\cr
[17]	M. N. Leger, A. G. Ryder, Applied Spectroscopy 2006, 60, 182.\cr                                                  
[18]	M. Clupek, P. Matejka, K. Volka, Journal of Raman Spectroscopy 2007, 38, 1174.\cr
[19]	Y. G. Hu, T. Jiang, A. G. Shen, W. Li, X. P. Wang, J. M. Hu, Chemometrics and Intelligent Laboratory Systems 2007, 85, 94.\cr
[20]	J. Zhao, H. Lui, D. I. McLean, H. Zeng, Applied Spectroscopy 2007, 61, 1225.\cr
[21]	C. Yang, Z. Y. He, W. C. Yu, Bmc Bioinformatics 2009, 10.\cr
[22]	P. Du, W. A. Kibbe, S. M. Lin, Bioinformatics 2006, 22, 2059.\cr
[23]	L. Nie, S. G. Wu, X. Q. Lin, L. Z. Zheng, L. Rui, Journal of Chemical Information and Computer Sciences 2002, 42, 274.\cr
[24]	X. G. Shao, C. X. Ma, Chemometrics and Intelligent Laboratory Systems 2003, 69, 157.\cr
[25]	X. G. Shao, A. K. M. Leung, F. T. Chau, Accounts of Chemical Research 2003, 36, 276.\cr
[26]	X. G. Shao, C. Y. Pang, Q. D. Su, Fresenius Journal of Analytical Chemistry 2000, 367, 525.\cr
[27]	E. T. Whittaker, Proceedings of the Edinburgh Mathematical Society 1922, 41, 63.\cr
[28]	P. H. C. Eilers, Analytical Chemistry 2003, 75, 3631.\cr
[29]	J. Carlos Cobas, M. A. Bernstein, M. Mart-Pastor, P. G. Tahoces, Journal of Magnetic Resonance 2006, 183, 145.\cr
}                                                                                                                      
\keyword{ package }
