segmut
================
2022-12-21

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
devtools::install_github("mishashe/segmut")
```

To load necessary libraries:

``` r
library(segmut)
library(RColorBrewer, quietly=T) # to plot results
library(stringr)
library(tidyverse)
library(scales)
library(ggExtra)
library(stringi)
library(ggplot2)
library(plotrix)
library(Rcpp)
library(RcppArmadillo)
library(roxygen2)
options(warn=-1)
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

## Example with unknown number of breaks

This is a basic example which shows you how to find optimal number of
breaks and their locations:

To find optimal number of breaks

``` r
breaks0L <- improve(muts, L,c(0,L)) # this is the main function (improve) of the package
```

To plot the results

``` r
# breaks0L <- sort(c(0,res,L))
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

To plot the density results

``` r
ma <- stats::filter(1:L %in% muts, rep(1 / 200, 200), sides = 2,circular=FALSE)
colors <- brewer.pal(name="Paired", n=length(breaks0L)-1)
plot(ma)   
axis(side=1, at=c(0,700,2000,3200,L))
for (i in 1:(length(breaks0L)-1))
{
  lines(c(breaks0L[i],breaks0L[i+1]),c(length(muts)/L,length(muts)/L),col=colors[i], lwd=5)
}
```

<img src="man/figures/README-plot density results-1.png" width="100%" height="50%" />

## Example with Escherichia coli (NZ_CP033020.1) vs. Salmonella enterica (NZ_AP026948) alignment

### Alignment of the genomes

We align two bacterial genomes using nucmer and get all the alignment
blocks and the mutations

``` r
datadir <- "/home/misha/Documents/Development/segmut/data/"
file1 <- paste0(datadir,"Escherichia_coli_1.fasta")
file2 <- paste0(datadir,"Salmonella_enterica_1.fasta")
# system(paste0("nucmer --mum --prefix=",datadir,"align ",file1," ",file2), intern = TRUE, wait=TRUE)
alignment <- system(paste0("show-aligns -w 500000 ", datadir,"align.delta NZ_CP033020.1 NZ_AP026948.1"), intern = TRUE, wait=TRUE)
alignment <- substr(alignment,12,nchar(alignment)) # removing empty space in the front

# replace long indels by 1bp indel to make sure that (1 indel = 1 mutation)
Ls <- c()
number_muts <- c()
r <- c()
mutsS <- list()
for (i in seq(11,length(alignment),9))
{
  str1 <- str_split(alignment[[i-2]],"")[[1]]
  str2 <- str_split(alignment[[i-1]],"")[[1]]
  strmuts <- str_split(alignment[[i]],"")[[1]]
  Ind <- which(str1[1:(length(str1)-1)]=="." & str1[2:length(str1)]==".")
  while (length(Ind)>0)
  {
    str1 <- str1[-Ind]
    str2 <- str2[-Ind]
    strmuts <- strmuts[-Ind]
    Ind <- which(str1[1:(length(str1)-1)]=="." & str1[2:length(str1)]==".")
  }
  Ind <- which(str2[1:(length(str2)-1)]=="." & str2[2:length(str2)]==".")
  while (length(Ind)>0)
  {
    str1 <- str1[-Ind]
    str2 <- str2[-Ind]
    strmuts <- strmuts[-Ind]
    Ind <- which(str2[1:(length(str2)-1)]=="." & str2[2:length(str2)]==".")
  }
  Ls <- c(Ls,length(strmuts))
  number_muts <- c(number_muts,sum(strmuts=="^"))
  r <- c(r,diff(c(0,which(strmuts=="^"),length(strmuts)+1))-1)
  mutsS[[length(mutsS)+1]] <- which(strmuts=="^")
}
block_divergences <- number_muts/Ls
divergence <- sum(number_muts)/sum(Ls)
```

Plotting divergences and lengths of the alignment blocks (not segmented
ones):

``` r
p <- ggplot(data=data.frame(Ls=Ls,block_divergences=block_divergences), aes(Ls, block_divergences)) + 
          geom_point(size=0.2) +
          theme_bw() +
          scale_x_continuous(trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          # scale_y_continuous(trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          ylab("average block divergence") + xlab("block length") 
  ggExtra::ggMarginal(p, type = "histogram")
```

<img src="man/figures/README-blocks-1.png" width="50%" height="50%" />

In total there are 739 blocks with the total length of 1769587 bp.

### Segmentation of the genomes

Segmenting the blocks using segmut package:

``` r
Ks <- c()
taus <- c()
nmutsS <- c()
start_time <- Sys.time()
k <- 0
for (i in order(-Ls))
{
  k <- k+1
  muts <- mutsS[[i]]
  L <- Ls[[i]]
  breaks0L <- improve(muts, L,c(0,L))
  for (j in 1:(length(breaks0L)-1))
  {
    Ks <- c(Ks, breaks0L[j+1]-breaks0L[j]+1)
    nmutsS <- c(nmutsS, sum(muts<breaks0L[j+1] & muts>breaks0L[j]))
  }
}
taus <- nmutsS/Ks
end_time <- Sys.time()
print(end_time - start_time)
#> Time difference of 0.2218208 secs
```

In total there are 2697 segments.

Plotting divergences and lengths of the segments:

``` r
p <- ggplot(data=data.frame(Ks=Ks,taus=taus), aes(Ks, taus)) + 
          geom_point(size=0.2) +
          theme_bw() +
          scale_x_continuous(trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          # scale_y_continuous(trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          ylab("segment divergence") + xlab("segment length") + 
          geom_smooth(method="loess",se=TRUE) 
  ggExtra::ggMarginal(p, type = "histogram")
```

<img src="man/figures/README-plot segments-1.png" width="50%" height="50%" />

Calculating empirical and pseudotheoretical MLDs:

``` r
rB <- c(seq(-0.5,35.5,3),10^seq(log10(35.5)+0.1,4.8,0.1))
rV <- 0.5*(rB[-length(rB)]+rB[-1])
p <- hist(r,breaks=rB,plot=FALSE)
mE <- p$counts/diff(rB)

rA <- 0:max(Ls+1)
mPT <- rA*0

for (ir in 1:length(rA)) 
{
    mPT[ir] <- sum((Ks>rA[ir])*(2*taus+(Ks-rA[ir]-1)*taus^2)*exp(log(1-taus)*rA[ir]))
}
for (iK in 1:length(Ks)) 
{
  mPT[Ks[iK]+1] <- mPT[Ks[iK]+1] + exp(log(1-taus[iK])*Ks[iK])
}
p <- weighted.hist(rA,w=mPT,breaks=rB,plot=FALSE)
mPT <- p$counts/diff(rB)



Ktot <- sum(Ls)
mExp <- (2*divergence+(Ktot-rA)*divergence^2)*(1-divergence)^(rA)
p <- weighted.hist(rA,w=mExp,breaks=rB,plot=FALSE)
mExp <- p$counts/diff(rB)

dat <- data.frame(rV=rV,mE=mE,mPT=mPT, mExp=mExp)
```

Plotting empirical and pseudotheoretical MLDs:

``` r
ggplot(data=dat, aes(rV, mE)) + 
          geom_line(aes(rV, mExp),col="grey", size=0.5) +
          geom_line(aes(rV, mPT),col="blue", size=0.7) +
          scale_x_continuous(limits = c(1e0, 3e3),trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          scale_y_continuous(limits = c(1e-4, 1e5),trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          geom_line(aes(rV,1e8/rV^3),linetype = "dashed") + 
          geom_line(aes(rV,1e11/rV^4),linetype = "dotted") + 
          geom_point(size=4,shape='o') +
      labs(x="r", y = "m(r)",
       title="match length distribution",
       caption="Here the circles correspond to the empirical MLD, grey line to the theoretical precition based on the average genome diversity, blue line to the theoretical precition based on the diversity of calculated segments. Dashed and dotted lines represent the -3 and -4 power-law accordingly.", width = 40) +
            theme_bw() +
    theme(
    plot.caption =
      ggtext::element_textbox(
        width = unit(1.0, "npc"),
        hjust = 0, vjust = 0,
        halign = 0))
```

<img src="man/figures/README-mld-1.png" width="50%" height="50%" />

## Segmentation of the human genome (chr1)

Segmenting using segmut package:

``` r
muts <- sort(read.table(paste0("/home/misha/Documents/Development/segmut/data/HG02715_1_1_250000000.txt"), header = FALSE, sep = "\t",row.names=NULL)$V1)
muts <- muts[muts<1e7]
L <- max(muts) + 1 
divergence <- length(muts)/L
Ks <- c()
taus <- c()
nmutsS <- c()
start_time <- Sys.time()
breaks0L <- improve(muts, L, c(0,L))
for (j in 1:(length(breaks0L)-1))
{
  Ks <- c(Ks, breaks0L[j+1]-breaks0L[j]+1)
  nmutsS <- c(nmutsS, sum(muts<breaks0L[j+1] & muts>breaks0L[j]))
}
taus <- nmutsS/Ks
end_time <- Sys.time()
print(end_time - start_time)
#> Time difference of 2.686884 mins

colors <- rep(brewer.pal(name="Paired", n=8),round(length(breaks0L)/8+1))
par(mar=c(2,0,0,0))
plot(density(muts, bw = 100),main="")   
for (i in 1:(length(breaks0L)-1))
{
  lines(c(breaks0L[i],breaks0L[i+1]),c(4e-8,4e-8),col=colors[i], lwd=5)
}
```

<img src="man/figures/README-human breaking to segments-1.png" width="100%" />

Plotting divergences and lengths of the human segments:

``` r
p <- ggplot(data=data.frame(Ks=Ks[Ks>100],taus=taus[Ks>100]), aes(Ks, taus)) + 
          geom_point(size=0.2) +
          theme_bw() +
          scale_x_continuous(trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          scale_y_continuous(limits = c(0, 0.01)) +
          ylab("segment divergence") + xlab("segment length") + 
          geom_smooth(method="loess",se=TRUE) 
  ggExtra::ggMarginal(p, type = "histogram")
```

<img src="man/figures/README-plot human segments-1.png" width="50%" height="50%" />

Calculating empirical and pseudotheoretical MLDs:

``` r
r <- diff(muts)
rB <- c(seq(-0.5,35.5,3),10^seq(log10(35.5)+0.1,7.8,0.1))
rV <- 0.5*(rB[-length(rB)]+rB[-1])
p <- hist(r,breaks=rB,plot=FALSE)
mE <- p$counts/diff(rB)

rA <- 0:max(r+1)
mPT <- rA*0
for (ir in 1:length(rA)) 
{
    mPT[ir] <- sum((Ks>rA[ir])*(2*taus+(Ks-rA[ir]-1)*taus^2)*exp(log(1-taus)*rA[ir]))
}
for (iK in 1:length(Ks)) 
{
  mPT[Ks[iK]+1] <- mPT[Ks[iK]+1] + exp(log(1-taus[iK])*Ks[iK])
}
mPT[is.na(mPT)] <- 0
p <- weighted.hist(rA,w=mPT,breaks=rB,plot=FALSE)
mPT <- p$counts/diff(rB)


mExp <- (2*divergence+(L-rA)*divergence^2)*(1-divergence)^(rA)
p <- weighted.hist(rA,w=mExp,breaks=rB,plot=FALSE)
mExp <- p$counts/diff(rB)

dat <- data.frame(rV=rV,mE=mE,mPT=mPT,mExp=mExp)
```

Plotting empirical and pseudotheoretical MLDs:

``` r
ggplot(data=dat, aes(rV, mE)) + 
          geom_line(aes(rV, mExp),col="grey", size=0.7) +
          geom_line(aes(rV, mPT),col="blue", size=0.7) +
          scale_x_continuous(limits = c(1e0, 1e5),trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          scale_y_continuous(limits = c(2e-6, 1e2),trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) +
          geom_line(aes(rV,1e11/rV^3),linetype = "dashed") + 
          geom_point(size=4,shape='o') +
      labs(x="r", y = "m(r)",
       title="match length distribution",
       caption="Here the circles correspond to the empirical MLD, grey line to the theoretical precition based on the average genome diversity, blue line to the theoretical precition based on the diversity of calculated segments. Dashed line represent the -3  power-law.", width = 40) +
            theme_bw() +
    theme(
    plot.caption =
      ggtext::element_textbox(
        width = unit(1.0, "npc"),
        hjust = 0, vjust = 0,
        halign = 0))
```

<img src="man/figures/README-mld chr1-1.png" width="50%" />
