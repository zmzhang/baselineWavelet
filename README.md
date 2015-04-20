# baselinewavelet

An intelligence background-correction algorithm for highly fluorescent sample in Raman spectroscopy has been developed with peak detection and width estimation by CWT wavelet and background fitting by penalized least squares. The programming language is R(http://www.r-project.org/).

## Installation

### Intall using devtools from github
library(devtools); 
install_github("zmzhang/baselineWavelet")

### Install from Local zip

Install the downloaded packages from local zip or tar.gz file.

To start running this algorithm, load the baselineWavelet package through "library(baselineWavelet)" in the R commandline windows, try "?baselineWavelet" in the R commandline windows to open the documents.

## Correction example
This is a correction example:

![Correction Example](/images/logo.jpg)

## Contact
For any questions, please contact:
Yi-Zeng Liang: yizeng_liang@263.net
or:
Zhi-Min Zhang: zhangzhimin.csu@gmail.com

## How to cite:
Z.M. Zhang, S. Chen, Y.Z. Liang, et al., An intelligent background-correction algorithm for highly fluorescent samples in Raman spectroscopy. Journal of Raman Spectroscopy 41 (6), 659 (2010).

[Download pdf and endnote citation here](http://www3.interscience.wiley.com/journal/122630376/abstract)

## Release note of baselineWavelet
What's news:
1. What New of baselineWavelet 4.0.2:
Running smoothly in R 3.0 and above.
2. What New of baselineWavelet 4.0.0:
By taking the advantage of sparse matrix in R package "Matrix", we implemented the sparse version of whittaker smoother and baselineWavelet alogrithm. Now the speed of baselineWavelet 4.0 is faster than baselineWavelet 3.0 by 100 times or more. And it was built with R 2.12.2.
3. From version 2.0 to 3.0: Rewirte the WhittakerSmooth? function, don't use the cholskey decomposition any more.
4. From version 1.0 - 2.0: Two functions, say baselineCorrectionCWT() and WhittakerSmoother?(), in the baselineWavelet package were modified (add a parameter) so that one could easily perform first, second or even higher differences penalties by adjusting the parameter for the purpose.


