segmut
================
2022-11-29

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
# Sys.unsetenv("GITHUB_PAT")
devtools::install_github("mishashe/segmut")
```

To load nessesary libraries:

``` r
library(segmut)
library(DEoptim, quietly=T) # to find optimal break points using differential evolution algorithm
#> 
#> DEoptim package
#> Differential Evolution algorithm in R
#> Authors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich
library(RColorBrewer, quietly=T) # to plot results
```

## Example with known number of breaks n=2

This is a basic example which shows you how to find positions of n=2
breaks:

To generate example of genome of length L with vector of mutation
locations muts:

``` r
L <- 10000
muts <- sort(c(sample(1:3000,400),sample(3001:8000,250),sample(8001:10000,400)))
```

To find optimal breaks given n=2 number of breaks

``` r
res <- getBreaks(muts = muts, L = L, Kmin=0, n=2)
```

To plot the results

``` r
breaks <- sort(c(0,res$optim$bestmem,L))
colors <- brewer.pal(name="Paired", n=length(breaks)-1)
par(mar=c(2,0,0,0))
plot(muts,rep(0,length(muts)),pch=".", cex = 1.5,ylim=c(-0.06,0.01),ylab="",xlab="", axes=F)
axis(side=1, at=c(0,3000,8000,L))
for (i in 1:(length(breaks)-1))
{
  lines(c(breaks[i],breaks[i+1]),c(-0.05,-0.05),col=colors[i], lwd=5)
}
```

<img src="man/figures/README-plot results-1.png" width="100%" />
