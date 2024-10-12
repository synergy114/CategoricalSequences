#ifndef RNG_INIT_H
#define RNG_INIT_H

#include <Rcpp.h>
#include "02_Initialize.h" 
#include "03_RandomWalkImpute.h"
#include "04_TransitionsCount.h"
#include "05_EM_Step.h"

using namespace Rcpp;

// Function declarations

NumericVector SumWeightedVectors(NumericVector weights, NumericMatrix matrix);

// [[Rcpp::export]]
NumericMatrix SumWeightedMatrices(NumericVector weights, List matrix_list);

// [[Rcpp::export]]
double RandomInitLikelihood(List sequences, int nCategories, int nClusters, int seed);

// Compute the Q-function (log-likelihood)
// [[Rcpp::export]]
double Q_function(int K, int nCategories, List current_params, List allImputes, List Pi, List XY1); 

List RNG_Init(List sequences, int nCategories, int nClusters);


#endif
