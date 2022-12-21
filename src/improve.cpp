// #include <Rcpp.h>
#include <iostream>
#include <random>
#include <chrono>
#include <math.h>
#include <stdint.h>

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>
//using namespace Rcpp;
using namespace arma;
using namespace std;

double getChiSquare(arma::rowvec par, arma::uvec& muts, double& L)
{
  double chiSquare = 0;
  double density = (double)muts.n_rows / (double)L;
  double expected;
  int imuts = 0;
  int ipar;
  int nmuts = 0;

  for (ipar = 1; ipar < par.n_cols; ipar++)
  {
    nmuts = 0;
    while (imuts<muts.n_rows && muts(imuts) < par(ipar))
    {
      nmuts++;
      imuts++;
    }
    expected = density*(par(ipar)-par(ipar-1)+1);
    chiSquare += (nmuts - expected)*(nmuts - expected)/expected;
  }
  return(chiSquare);
}

long int getBestSingleBreak(arma::uvec& muts, double& L, arma::rowvec& parPrev)
{
  double expected;
  double expected1;
  double expected2;
  double bestChiSquareDiff = log((double)muts.n_rows);
  double ChiSquare;
  double ChiSquare1;
  double ChiSquare2;
  long int xBest=0;
  double density = (double)muts.n_rows / (double)L;

  int nmuts = 0;
  int imuts = 0;
  int ipar = 0;
  int imutsx = 0;
  int nmutsx = 0;
  for (ipar = 1; ipar < parPrev.n_cols; ipar++)
  {
    nmuts = 0;
    while (imuts<muts.n_rows && muts(imuts) < parPrev(ipar))
    {
      nmuts++;
      imuts++;
    }
    expected = density*(parPrev(ipar)-parPrev(ipar-1)+1);
    ChiSquare = (nmuts - expected)*(nmuts - expected)/expected;
    nmutsx = 0;
    for (long int x = parPrev(ipar-1); x <= parPrev(ipar); x++)
    {
      while (imutsx<muts.n_rows && muts(imutsx) < x)
      {
        nmutsx++;
        imutsx++;
      }
      if (x > parPrev(ipar-1) + 10 && x < parPrev(ipar)-10 && nmuts>5)
      {
        expected1 = density*(x-parPrev(ipar-1)+1);
        ChiSquare1 = (nmutsx - expected1)*(nmutsx - expected1)/expected1;
        expected2 = density*(parPrev(ipar) - x + 1);
        ChiSquare2 = (nmuts - nmutsx - expected2)*(nmuts - nmutsx - expected2)/expected2;
        if (ChiSquare1 + ChiSquare2 - ChiSquare > bestChiSquareDiff)
        {
          bestChiSquareDiff = ChiSquare1 + ChiSquare2 - ChiSquare;
          xBest = x;
        }
      }
    }
  }
  return(xBest);
}

// [[Rcpp::depends(RcppArmadillo)]]


//' function to improve the breaks.
//'
//' @param muts vector of mutations
//' @param L genome length
//' @param par breaks (including 0 and L)
//' @return breaks
// [[Rcpp::export]]
arma::rowvec improve(arma::uvec& muts, double L=0, arma::rowvec par={})
{
  if (L==0) L = max(muts);
  if (par.n_cols==0) par = {0,L};
  double lognmuts = log((double)muts.n_rows);
  double ChiSquareBest = -1000.0;
  int runs=0;
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist(0,1e6); // distribution in range
  long int newPos;

  newPos = getBestSingleBreak(muts, L, par);
  if (newPos==0 && par.n_cols==2) //cannot add even a single break
  {
    return(par);
  }
  else
  {
    arma::rowvec parAdded(par.n_cols+1);
    for (int i=0; i<par.n_cols; i++) parAdded(i) = par(i);
    parAdded(par.n_cols)=newPos;
    par = sort(parAdded);
  }

  while (runs < par.n_cols)
  {
    newPos = getBestSingleBreak(muts, L, par);
    if (newPos!=0) //add one break
    {
      arma::rowvec parAdded(par.n_cols+1);
      for (int i=0; i<par.n_cols; i++) parAdded(i) = par(i);
      parAdded(par.n_cols)=newPos;
      par = sort(parAdded);
    }

    int Ind = 1 + (dist(rng) % (par.n_cols-2));
    arma::rowvec parRemoved(par.n_cols-1);
    for (int i=0;i<Ind;i++) parRemoved(i) = par(i);
    for (int i=Ind+1;i<par.n_cols;i++) parRemoved(i-1) = par(i);
    newPos = getBestSingleBreak(muts, L, parRemoved);
    if (newPos==0) //remove one break
    {
      par = parRemoved;
    }
    else if (newPos!=par(Ind))  //move one break
    {
      par(Ind) = newPos;
      par = sort(par);
    }
    double ChiSquare = getChiSquare(par, muts, L) - (par.n_cols-2.0)*lognmuts;
    if (ChiSquare > ChiSquareBest)
    {
      ChiSquareBest = ChiSquare;
      runs = 0;
    }
    runs++;
  }
  return(par);
}










