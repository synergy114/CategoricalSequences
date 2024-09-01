//#include <Rcpp.h>
#include "02_Initialize.h"  // Include the header file
//using namespace Rcpp;


////Initialize clusters randomly
// [[Rcpp::export]]
IntegerVector InitilizeClusters(List sequences, int nClusters){
  int nSequences=sequences.size();
  IntegerVector clusterAssignments(nSequences);
  IntegerVector clusters = seq(0, nClusters - 1);
  clusterAssignments = Rcpp::sample(clusters, nSequences, true);
  return clusterAssignments;
}

//// Extract each cluster's sequences
// [[Rcpp::export]]
List filterSequencesByCluster(List sequences, IntegerVector clusters, int k) {
  int n = sequences.size();
  List filtered_sequences;

  for (int i = 0; i < n; ++i) {
    if (clusters[i] == k) {
      filtered_sequences.push_back(sequences[i]);
    }
  }
  return filtered_sequences;
}

/// Initial State Prob
// [[Rcpp::export]]
NumericVector InitialStateProb(List sequences, int categories) { //probability for the 1st element only
  NumericVector probs(categories); // Initialize a vector to store counts of each category

  // Loop through each sequence in the list
  for (int i = 0; i < sequences.size(); i++) {
    IntegerVector seq = sequences[i]; // Convert the list element to IntegerVector
    if (!IntegerVector::is_na(seq[0])) {
      probs[seq[0] - 1] += 1; // Increment the count for the first element (adjusted for 0-based indexing)
    }
  }

  // Normalize the counts to get probabilities
  probs = probs / sum(probs);
  return probs;
}


// [[Rcpp::export]]
NumericMatrix InitialTransitionMatrix(List sequence, int nCategories) {
  // Initialize the transition matrix with zeros
  NumericMatrix transitionLocal(nCategories, nCategories);
  
  // Loop through the sequence and count transitions
  for (int seq_index = 0; seq_index < sequence.size(); seq_index++) {
    IntegerVector seq = sequence[seq_index];
    for (int i = 1; i < seq.size(); i++) {
      if (!IntegerVector::is_na(seq[i - 1]) && !IntegerVector::is_na(seq[i])) {
        int fromState = seq[i - 1] - 1; // Convert to 0-based index
        int toState = seq[i] - 1;       // Convert to 0-based index
        
        // Ensure indices are within bounds
        if (fromState >= 0 && fromState < nCategories && toState >= 0 && toState < nCategories) {
          transitionLocal(fromState, toState) += 1;
        }
      }
    }
  }
  
  // Convert frequencies to probabilities for each row
  for (int i = 0; i < nCategories; i++) {
    double total_transitions = sum(transitionLocal(i, _));
    
    // Normalize the row to convert frequencies to probabilities
    for (int j = 0; j < nCategories; j++) {
      if (total_transitions > 0) {
        transitionLocal(i, j) /= total_transitions;
      } else {
        transitionLocal(i, j) = 0.0; // Handle division by zero case
      }
    }
  }
  
  return transitionLocal;
}
// 
// //Initialize parameters
// [[Rcpp::export]]
List InitializeParameters(List sequences, int nClusters, int nCategories) {
  IntegerVector assigned_clusters = InitilizeClusters(sequences, nClusters);
  
  // Initialize mixing proportions to be equal for all clusters
  NumericVector mixing_prop(nClusters);
  for (int i = 0; i < nClusters; ++i) {
    mixing_prop[i] = 1.0 / nClusters;
  }
  
  // Initialize initial state probabilities and transition matrices
  NumericMatrix initial_values_alpha(nClusters, nCategories);
  List transition_matrices(nClusters);
  List grouped_sequences_by_clusters(nClusters);
  
  for (int cluster_i = 0; cluster_i < nClusters; ++cluster_i) {
    grouped_sequences_by_clusters[cluster_i] = filterSequencesByCluster(sequences, assigned_clusters, cluster_i);
    initial_values_alpha(cluster_i, _) = InitialStateProb(grouped_sequences_by_clusters[cluster_i], nCategories);
    transition_matrices[cluster_i] = InitialTransitionMatrix(grouped_sequences_by_clusters[cluster_i], nCategories);
  }
  
  return List::create(
    //Named("assigned_clusters") = assigned_clusters,
    Named("mixing_prop") = mixing_prop,
    Named("initial_values_alpha") = initial_values_alpha,
    Named("transition_matrices") = transition_matrices//,
    //Named("grouped_sequences_by_clusters") = grouped_sequences_by_clusters
  );
}





