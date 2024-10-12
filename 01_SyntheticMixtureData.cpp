#include "01_SyntheticMixtureData.h"

// Function to generate a single Markov sequence
IntegerVector generateMarkovSequence(int length, NumericMatrix transitionMatrix, int initialState) {
  IntegerVector sequence(length);
  int currentState = initialState;
  
  for (int i = 0; i < length; i++) {
    sequence[i] = currentState + 1; // Store 1-based state
    NumericVector prob = transitionMatrix(currentState, _);
    IntegerVector states = seq(0, prob.size() - 1);
    currentState = Rcpp::sample(states, 1, true, prob)[0];
  }
  return sequence;
}

// Function to generate synthetic mixture categorical sequence data
// [[Rcpp::export]]
List generateSyntheticData(int nSequences,int sequenceLengthMin, int sequenceLengthMax, int nClusters, NumericVector mixingProbs, List transitionMatrices, NumericMatrix initialProbabilities) {
  List data(nSequences);
  IntegerVector clusterAssignments(nSequences);
  
  // Check that the number of clusters matches the size of transitionMatrices and initialProbabilities
  if (transitionMatrices.size() != nClusters || initialProbabilities.nrow() != nClusters) {
    stop("Number of clusters does not match the size of transitionMatrices or initialProbabilities");
  }
  
  // Generate sequences
  for (int i = 0; i < nSequences; i++) {
    IntegerVector clusters = seq(0, nClusters - 1);
    int cluster = Rcpp::sample(clusters, 1, true, mixingProbs)[0];
    int sequenceLength;
    if (sequenceLengthMax!=sequenceLengthMin) 
      {
      IntegerVector seqLengths=seq(sequenceLengthMin, sequenceLengthMax);
      sequenceLength=Rcpp::sample(seqLengths, 1, true)[0];
      }
    else 
      {sequenceLength=sequenceLengthMax;}
    
    clusterAssignments[i] = cluster + 1; // Store 1-based cluster
    NumericMatrix transitionMatrix = transitionMatrices[cluster];
    NumericVector initialProb = initialProbabilities(cluster, _);
    IntegerVector states = seq(0, initialProb.size() - 1);
    int initialState = Rcpp::sample(states, 1, true, initialProb)[0];
    data[i] = generateMarkovSequence(sequenceLength, transitionMatrix, initialState);
  }
  
  // Return the generated data and cluster assignments
  return List::create(
    Named("sequences") = data,
    Named("clusterAssignments") = clusterAssignments
  );
}


// Function to randomly impute NA values in a numeric vector
IntegerVector RandomImputeVector(IntegerVector x, double prob) {
  int len = x.size();
  int imputation_size = static_cast<int>(len * prob);
  
  IntegerVector imputation_places = Rcpp::sample(len, imputation_size, false);
  // Sort the sampled positions to ensure correct checking for consecutive NAs
  std::sort(imputation_places.begin(), imputation_places.end());
  
  for (int i = 0; i < imputation_places.size(); i++) {
    x[imputation_places[i]-1] = NA_INTEGER;
  }
  // Rcpp::Rcout << "Imputation places: " << imputation_places << std::endl;
  
  return x;
}

// Function to randomly impute NA values in a list of numeric vectors
// [[Rcpp::export]]
List RandomImputeList(List seq_list, double prob) {
  List imputed_list = clone(seq_list);
  for (int i = 0; i < imputed_list.size(); ++i) {
    IntegerVector seq = as<IntegerVector>(imputed_list[i]);
    imputed_list[i] = RandomImputeVector(seq, prob);
  }
  
  return imputed_list;
}

