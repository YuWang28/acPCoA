
<!-- README.md is generated from README.Rmd. Please edit that file -->

# acPCoA

<!-- badges: start -->
<!-- badges: end -->

AC-PCoA is a method proposed by Yu Wang etc., which reduces the data
dimension while extracting the information from different distance
measures using principal coordinate analysis (PCoA), and adjusts the
confounding factors across multiple data sets  
by minimizing the associations between the lower dimensional
representations and the confounding variables. Application of the
proposed method is further extended to the scenario of classification
and prediction.

## Installation

To install “acPCoA”, first you need to install the P package “acPCA”
from <https://github.com/linzx06/AC-PCA>

Then, you can install the released version of “acPCoA” from github with

``` r
#install.packages("devtools")
#library(devtools)
#install_github("YuWang28/acPCoA")
```

## Example

This is a basic example which shows you how to implement acPCoA for
visualization after confounding factor adjustment:

``` r
library(acPCoA)
library(ggplot2)
X <- data_mbqc_groupA$DistMat.BC;
Y <- data_mbqc_groupA$ConfounderMat;
result_acPCoA <- acPCoA(DistanceMatrix=X, ConfounderMatrix=Y, nPC=2, lambda=seq(0, 20, 0.05), kernel="linear")
ggplot(as.data.frame(result_acPCoA$Xv),aes(x=V1,y=V2,color=data_mbqc_groupA$Specimen))+geom_point()
```

<img src="man/figures/README-example-1.png" width="100%" />
