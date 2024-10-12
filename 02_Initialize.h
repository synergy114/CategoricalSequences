#ifndef Initialize_h
#define Initialize_h

#include <Rcpp.h>
using namespace Rcpp;

// Function declarations
IntegerVector InitilizeClusters(Rcpp::List sequences, int nClusters, int seed);
List filterSequencesByCluster(Rcpp::List sequences, Rcpp::IntegerVector clusters, int k);
NumericVector InitialStateProb(Rcpp::List sequences, int categories);
NumericMatrix InitialTransitionMatrix(Rcpp::List sequence, int nCategories);
List InitializeParameters(Rcpp::List sequences, int nClusters, int nCategories, int seed);

#endif