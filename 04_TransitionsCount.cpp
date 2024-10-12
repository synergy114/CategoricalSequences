#include "04_TransitionsCount.h"

// [[Rcpp::export]] 
int rowSum(IntegerMatrix mat, int row) {
  
  int nCols = mat.ncol();
  int rowSum= 0.0;
  for (int j = 0; j < nCols; j++) {
    rowSum += mat(row, j);
  }
  
  return rowSum;
}

/////////Transitions count for the given imputation
NumericMatrix Transitions_X_iu(NumericVector sequence, int nCategories, 
                               IntegerVector missing_positions, IntegerVector imputation) {
  
  // Make a copy of sequence to avoid modifying the original
  NumericVector y = clone(sequence); 
  
  // Replace missing positions with imputations
  for(int i = 0; i < missing_positions.size(); ++i) {
    y[missing_positions[i]] = imputation[i];
  }
  
  // Initialize the transition matrix
  NumericMatrix transitions(nCategories, nCategories);
  
  // Compute the transition matrix directly based on state indices
  for(int i = 1; i < y.size(); ++i) {
    int from = y[i - 1] - 1; // Assuming states are 1-indexed (1, 2, 3, ..., nCategories)
    int to = y[i] - 1;
    
    if(from >= 0 && from < nCategories && to >= 0 && to < nCategories) {
      transitions(from, to) += 1;
    }
  }
  return transitions;
}

/////////Transition matrices for one sequence, list of all transitions corresponding for each imputation uses Transitions_X_iu for each imputation
// // [[Rcpp::export]]
List Transitions_i(NumericVector sequence, int nCategories, IntegerVector missing_positions_i, List imputations_i){
  List x_i;
  IntegerVector y_i1;
  for (int u=0; u< imputations_i.size();u++){
    IntegerVector impute=imputations_i[u];
    NumericMatrix xiu=Transitions_X_iu(sequence, nCategories, 
                                       missing_positions_i,impute);
    x_i.push_back(xiu); 
    if (missing_positions_i[0]==0){
      y_i1.push_back(impute[0]);
    }
    if (!(missing_positions_i[0]==0)){
      y_i1.push_back(sequence[0]); 
    }
  }
  return List::create(
    Named("Y_i1") = y_i1,
    Named("x_i") = x_i);
} 

// [[Rcpp::export]] 
List TransitionsCount_Y1(List sequences,int nCategories, List all_imputes){
  List transitionX_y1;
  int n=sequences.size();
  for (int i=0; i<n; i++){
    List current_seq_impute=all_imputes[i]; 
    IntegerVector missing_positions=current_seq_impute["missing_positions"];
    List imputations=current_seq_impute["unique_imputes"]; 
    List xy1=Transitions_i(sequences[i], nCategories, missing_positions, imputations);
    transitionX_y1.push_back(xy1);
  }
  return transitionX_y1;
}


