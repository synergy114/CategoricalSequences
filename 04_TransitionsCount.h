#ifndef TRANSITIONS_COUNT_H
#define TRANSITIONS_COUNT_H

#include <Rcpp.h>
using namespace Rcpp;

// Function declarations
int rowSum(IntegerMatrix mat, int row);

NumericMatrix Transitions_X_iu(NumericVector sequence, int nCategories, 
                               IntegerVector missing_positions, IntegerVector imputation);

List Transitions_i(NumericVector sequence, int nCategories, IntegerVector missing_positions_i, List imputations_i);

List TransitionsCount_Y1(List sequences,int nCategories, List all_imputes);

#endif
