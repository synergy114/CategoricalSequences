#ifndef SYNTHETIC_MIXTURE_DATA_H
#define SYNTHETIC_MIXTURE_DATA_H

#include <Rcpp.h>
using namespace Rcpp;

// Function to generate a single Markov sequence
IntegerVector generateMarkovSequence(int length, NumericMatrix transitionMatrix, int initialState);

// Function to generate synthetic mixture categorical sequence data
// [[Rcpp::export]]
List generateSyntheticData(int nSequences, int sequenceLengthMin, int sequenceLengthMax,
                           int nClusters, NumericVector mixingProbs,
                           List transitionMatrices, NumericMatrix initialProbabilities);

// Function to randomly impute NA values in a numeric vector
IntegerVector RandomImputeVector(IntegerVector x, double prob);

// Function to randomly impute NA values in a list of numeric vectors
// [[Rcpp::export]]
List RandomImputeList(List seq_list, double prob);

#endif // SYNTHETIC_MIXTURE_DATA_H
