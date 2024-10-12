#ifndef EM_STEP_H
#define EM_STEP_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "04_TransitionsCount.h"

using namespace Rcpp;

// Function declarations
void debug_print(std::string msg, int value); 
// Calculate pi_i for each sequence

NumericVector SumWeightedVectors(NumericVector weights, NumericMatrix matrix);

NumericMatrix SumWeightedMatrices(NumericVector weights, List matrix_list);

NumericMatrix Pi_i(IntegerVector y_i1, int nCategories, NumericVector tau,  
                   NumericMatrix initial_prob, List transition_matrices, List x_i);

// Expectation step: Collect all Pi

List E_Step(int nCategories, Rcpp::List XY1, Rcpp::List current_params);

// Maximization step: Update tau, alpha, and transition matrices
List M_Step(int K,int nCategories, List allImputes, List Pi, List XY1);

double Q_function(int K, int nCategories, List current_params, List allImputes, List Pi, List XY1); 


#endif
