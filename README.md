segmut
================
2022-12-10

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
# install.packages("/home/misha/Documents/Development/segmut/", repos = NULL, type = "source")
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
library(stringr)
```

To set parameters

``` r
Kmin <- 100
```

To generate example of genome of length `L` with vector of mutation
locations `muts`:

``` r
L <- 4000
muts <- sort(c(sample(1:700,700*0.1),
               sample(701:2000,2300*0.2),
               sample(2001:3200,1200*0.15),
               sample(3201:4000,800*0.08)))
```

## Example with known number of breaks `n=3`

This is a basic example which shows you how to find positions of `n=3`
breaks:

To find optimal breaks locations given `n=3` number of breaks

``` r
res <- getBreaks(muts = muts, L = L, Kmin=Kmin, n=3)
```

To plot the results

``` r
breaks0L <- sort(c(0,res$optim$bestmem,L))
colors <- brewer.pal(name="Paired", n=length(breaks0L)-1)
par(mar=c(2,0,0,0))
plot(muts,rep(0,length(muts)),pch=".", cex = 1.5,ylim=c(-0.06,0.01),ylab="",xlab="", axes=F)
axis(side=1, at=c(0,700,2000,3200,L))
for (i in 1:(length(breaks0L)-1))
{
  lines(c(breaks0L[i],breaks0L[i+1]),c(-0.05,-0.05),col=colors[i], lwd=5)
}
```

<img src="man/figures/README-plot results-1.png" width="100%" />

## Example with unknown number of breaks

This is a basic example which shows you how to find optimal number of
breaks and their locations:

To find optimal number of breaks

``` r
breaks <- getNumberBreaks(muts,L=L,Kmin=Kmin,pThreshold=0.05)
#> [1] 0
#> [1] 1.396169e-09 5.180829e-01
#> [1] 0.471251747 0.064614271 0.001476743
#> [1] 0.8258191 0.9764116 0.6269692 0.4385365
```

To plot the results

``` r
breaks0L <- sort(c(0,breaks,L))
colors <- brewer.pal(name="Paired", n=length(breaks0L)-1)
par(mar=c(2,0,0,0))
plot(muts,rep(0,length(muts)),pch=".", cex = 1.5,ylim=c(-0.06,0.01),ylab="",xlab="", axes=F)
axis(side=1, at=c(0,700,2000,3200,L))
for (i in 1:(length(breaks0L)-1))
{
  lines(c(breaks0L[i],breaks0L[i+1]),c(-0.05,-0.05),col=colors[i], lwd=5)
}
```

<img src="man/figures/README-plot best results-1.png" width="100%" />

## Example with Escherichia coli (NZ_CP033020.1) vs. Salmonella enterica (NZ_AP026948) alignment

We align two bacterial genomes using nucmer and get all the alignment
blocks and the mutations

``` r
datadir <- "/home/misha/Documents/Development/segmut/data/"
file1 <- paste0(datadir,"Escherichia_coli_1.fasta")
file2 <- paste0(datadir,"Salmonella_enterica_1.fasta")
system(paste0("nucmer --mum --prefix=",datadir,"align ",file1," ",file2), intern = TRUE, wait=TRUE)
#> character(0)
alignment <- system(paste0("show-aligns -w 500000 ", datadir,"align.delta NZ_CP033020.1 NZ_AP026948.1"), intern = TRUE, wait=TRUE)
alignment <- alignment[str_detect(alignment,"\\^")]
divergence <- sum(str_count(alignment,"\\^"))/sum(nchar(alignment))
divergences <- (str_count(alignment,"\\^"))/(nchar(alignment))
```
