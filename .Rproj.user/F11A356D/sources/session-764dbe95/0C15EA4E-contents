#include "06_RNG_Init.h"

const double doubleMax = 1e+10;
NumericVector SumWeightedVectors(NumericVector weights, NumericMatrix matrix){
  int n=weights.size();
  int m = matrix.ncol();  // Get the number of columns in the matrix
  
  NumericVector weighted_vector(m);
  for (int i=0; i<n; i++){
    weighted_vector += weights[i] * matrix(i, _);
  }
  return weighted_vector;
}

// [[Rcpp::export]]
NumericMatrix SumWeightedMatrices(NumericVector weights, List matrix_list) {
  int n = matrix_list.size(); // Number of matrices
  if (n == 0) return NumericMatrix(0, 0); // Return an empty matrix if the list is empty
  
  // Get the dimensions of the first matrix
  NumericMatrix first_matrix = matrix_list[0];
  int rows = first_matrix.nrow();
  int cols = first_matrix.ncol();
  
  // Initialize the result matrix with zeros
  NumericMatrix result(rows, cols);
  
  // Loop through each matrix and weight
  for (int i = 0; i < n; i++) {
    NumericMatrix current_matrix = matrix_list[i];
    double current_weight = weights[i];
    
    // Element-wise addition of the weighted matrices
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        result(r, c) += current_weight * current_matrix(r, c);
      }
    }
  }
  
  return result;
}


// [[Rcpp::export]]
double RandomInitLikelihood(List sequences, int nCategories, int nClusters, int seed){
  
  List init_params=InitializeParameters(sequences, nCategories, nClusters,seed);
  //int n=sequences.size();
  
  List current_imputes= AllSeqJointImputeRCPP(sequences,init_params, 100); 
  List XY1=TransitionsCount_Y1(sequences,nCategories,current_imputes);

  List Pi= E_Step(nCategories, XY1, init_params);
  List current_params=M_Step(nClusters, nCategories, current_imputes, Pi, XY1);
  
  double Q= Q_function(nClusters,nCategories, current_params, current_imputes,  Pi,  XY1);
  return Q;
  
}

// [[Rcpp::export]]
List RNG_Init(List sequences, int nClusters, int nCategories) {

  // Set the random seed
  Rcpp::RNGScope scope;
  std::srand(42);

  // Number of samples
  int num_samples = 100;
  int min_value = 0;
  int max_value = 10000;
  int top_n = 10;

  NumericVector Q_values(num_samples);
  // Sample 100 values from the vector without replacement
  IntegerVector range = seq(min_value, max_value);
  IntegerVector seeds = Rcpp::sample(range, num_samples, false); // Range from which to sample

  // Generate Q_values using RandomInitLikelihood for each seed
  for (int s = 0; s < num_samples; s++) {
    Q_values[s] = RandomInitLikelihood(sequences, nCategories, nClusters, seeds[s]);
  }

  // Create a DataFrame with seeds and Q_values
  DataFrame df = DataFrame::create(Named("Seed") = seeds,
                                   Named("Q_value") = Q_values);

  // Sort the data frame by Q-function results in descending order
  Rcpp::IntegerVector idx = Rcpp::seq(0, df.nrow() - 1);

  // Sort the index vector based on the Q_values in descending order
  std::sort(idx.begin(), idx.end(),
            [&Q_values](int i, int j) { return Q_values[i] > Q_values[j]; });

  // Store the top `top_n` seeds and Q-values
  List good_seeds;
  List top_Q;

  for (int k = 0; k < top_n; k++) {
    int sorted_idx = idx[k]; // Get the sorted index
    good_seeds.push_back(seeds[sorted_idx]); // Use sorted index to access the correct seed
    top_Q.push_back(Q_values[sorted_idx]);   // Use sorted index to access the correct Q-value

    // Debugging: Print the sorted likelihoods for verification
    debug_print("Highest Likelihood", Q_values[sorted_idx]);
  }

  return good_seeds;
}

