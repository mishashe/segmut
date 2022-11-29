segmut
================
2022-11-29

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package description

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
library(doParallel, quietly=T)
ncores <- 4 # set to desired number of cores for parallel computation
registerDoParallel(cores = ncores)
library(rlist)
```

To set parameters

``` r
Kmin <- 300
```

## Example with known number of breaks `n=2`

This is a basic example which shows you how to find positions of `n=2`
breaks:

To generate example of genome of length `L` with vector of mutation
locations `muts`:

``` r
L <- 10000
muts <- sort(c(sample(1:3000,3000*0.1),sample(3001:8000,5000*0.18),sample(8001:10000,2000*0.12)))
```

To find optimal breaks locations given `n=2` number of breaks

``` r
res <- getBreaks(muts = muts, L = L, Kmin=Kmin, n=2)
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

## Example with unknown number of breaks

This is a basic example which shows you how to find optimal number of
breaks and their locations:

Gets $\log_{10}$ *p*-value vector for randomly shuffled mutations and
one break and estimates *p*-value of the null model as the 5% quantile
of this vector

``` r
p0 <- get_p0(muts,L=L,Kmin=Kmin,A=20)
```

Calculate p-value for `n=1`:

``` r
res <- getBreaks(muts = muts, L = L, Kmin=Kmin, n=1)
if (res$optim$bestval>p0) {print("There is no support even for a single break")} else
 {print("There is support for one break, let's test for more")}
#> [1] "There is support for one break, let's test for more"
```

To find optimal number of breaks

``` r
resList <- list(res,getBreaks(muts = muts, L = L, Kmin=Kmin, n=2))
imax <- length(resList)
while (resList[[imax]]$optim$bestval < resList[[imax-1]]$optim$bestval) 
{
  resList <- list.append(resList, getBreaks(muts = muts, L = L, Kmin=Kmin, n=length(resList)+1))
  imax <- length(resList)
  print(imax)
}
#> [1] 3
p.values <- sapply(1:length(resList),function(i){resList[[i]]$optim$bestval})
plot(1:length(resList),p.values)
```

<img src="man/figures/README-find optimal number of breaks-1.png" width="100%" />

``` r
n <- which.min(p.values)
res <- resList[[n]]
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

<img src="man/figures/README-plot best results-1.png" width="100%" />
