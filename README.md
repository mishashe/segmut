
<!-- README.md is generated from README.Rmd. Please edit that file -->

# segmut

<!-- badges: start -->
<!-- badges: end -->

The goal of segmut is to segment genome with mutations to segments with
significantly different mutations density.

## Installation

You can install the development version of segmut from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mishashe/segmut")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(segmut)
library(DEoptim)
#> Loading required package: parallel
#> 
#> DEoptim package
#> Differential Evolution algorithm in R
#> Authors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich
muts <- sort(c(sample(1:3000,400),sample(3001:8000,250),sample(8001:10000,400)))
L <- 10000
# res <- getBreaks(muts = muts, L = L, Kmin=0)
# plot(muts,rep(0,length(muts)),pch=".")
# breaks <- sort(c(0,res$optim$bestmem,L))
# colors <- brewer.pal(name="Paired", n=length(breaks)-1)
# for (i in 1:(length(breaks)-1))
# {
#   lines(c(breaks[i],breaks[i+1]),c(-0.05,-0.05),col=colors[i], lwd=5)
# }
```
