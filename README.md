
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

You can install the released version of acPCoA from github with

``` r
library(devtools)
#> Loading required package: usethis
install_github("YuWang28/acPCoA")
#> Downloading GitHub repo YuWang28/acPCoA@HEAD
#> Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
#>   download from 'https://api.github.com/repos/YuWang28/acPCoA/tarball/HEAD' failed
```

## Example

This is a basic example which shows you how to implement acPCoA for
visualization after confounding factor adjustment:

``` r
library(acPCoA)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
