#include <Rcpp.h>
using namespace Rcpp;

///////////////////////////////////////////////
// [[Rcpp::export]]
////Function to simulate 1-step random walk
int RMultinom(NumericVector probs) {
  int k = probs.size();
  IntegerVector ans(k);
  int impute_value = -1;  // Initialize impute_value to a sentinel value
  
  // Simulate a single draw from the multinomial distribution
  R::rmultinom(1, probs.begin(), k, ans.begin());
  
  // Find the index of the first non-zero element in the result
  for (int i = 0; i < k; i++) {
    if (ans[i] > 0) {
      impute_value = i;  // Store 0-based index of non-zero element
      break;  // Stop searching after finding the first non-zero element
    }
  }
  
  // Return the 1-based index of the first non-zero element
  return impute_value + 1;
}

// Function to perform group random walk imputation
// [[Rcpp::export]]
List GroupRandomWalk(IntegerVector x, IntegerVector missing_positions, NumericVector initial_prob,  NumericMatrix transition_matrix, int num_sim) {
  int i = 0;
  int prev_state;
  int next_state = -1;
  int original_prev_state = -1;
  int group_size = missing_positions.size();
  
  // Debugging: Print initial states and sizes
  Rcpp::Rcout << "Group size: " << group_size << std::endl;
  Rcpp::Rcout << "Missing positions: " << missing_positions << std::endl;
  Rcpp::Rcout << "Sequence length: " << x.size() << std::endl;
  
  List imputes_list;
  IntegerVector imputes_count;
  
  bool first_element_missing = (missing_positions[0] == 0);
  bool last_element_missing = (missing_positions[group_size - 1] == x.size() - 1);
  
  if (!first_element_missing) {
    original_prev_state = x[missing_positions[0] - 1];
  }
  if (!last_element_missing) {
    // Debugging: Check bounds before access
    if (missing_positions[group_size - 1] + 1 >= x.size()) {
      Rcpp::stop("Error: Index out of bounds when accessing next_state");
    }
    next_state = x[missing_positions[group_size - 1] + 1];
  }
  
  while (i < num_sim) {
    IntegerVector sim_impute(group_size);
    
    if (first_element_missing) {
      sim_impute[0] = RMultinom(initial_prob);
    } else {
      sim_impute[0] = RMultinom(transition_matrix(original_prev_state - 1, _));
    }
    prev_state = sim_impute[0];
    int current_transition_state = prev_state - 1;
    
    for (int j = 1; j < group_size; j++) {
      // Debugging: Check transition state before access
      if (current_transition_state < 0 || current_transition_state >= transition_matrix.nrow()) {
        Rcpp::Rcout << "Error: current_transition_state out of bounds" << std::endl;
        Rcpp::Rcout << "current_transition_state: " << current_transition_state << std::endl;
        Rcpp::Rcout << "transition_matrix.nrow(): " << transition_matrix.nrow() << std::endl;
        Rcpp::stop("Error in transition matrix access");
      }
      
      int RW_single_value = RMultinom(transition_matrix(current_transition_state, _));
      sim_impute[j] = RW_single_value;
      current_transition_state = RW_single_value - 1;
    }
    
    if (last_element_missing || (!last_element_missing && RMultinom(transition_matrix(current_transition_state, _)) == next_state)) {
      bool found_in_list = false;
      
      for (int k = 0; k < imputes_list.size(); k++) {
        if (Rcpp::is_true(all(as<IntegerVector>(imputes_list[k]) == sim_impute))) {
          imputes_count[k] += 1;
          found_in_list = true;
          break;
        }
      }
      
      if (!found_in_list) {
        imputes_list.push_back(clone(sim_impute));
        imputes_count.push_back(1);
      }
      
      i += 1;
    }
  }
  
  NumericVector probabilities = as<NumericVector>(imputes_count) / num_sim;
  
  Rcpp::List GroupRW = Rcpp::List::create(
    Rcpp::Named("missing_positions") = missing_positions,
    Rcpp::Named("unique_imputes") = imputes_list,
    Rcpp::Named("probabilities") = probabilities
  );
  
  return GroupRW;
}

// [[Rcpp::export]]
List TestGroupRandomWalk() {
  IntegerVector x = {NA_INTEGER,2,1, 1, 2, NA_INTEGER, 2}; // The sequence with missing values
  IntegerVector missing_positions = {0, 5};  // Positions where x has NA
  NumericVector initial_prob = {0.5, 0.5};  // Equal probability to start with state 1 or 2
  NumericMatrix transition_matrix(2, 2);  // Transition matrix for 2 states
  transition_matrix(0, 0) = 0.7; transition_matrix(0, 1) = 0.3; // From state 1 to 1 or 2
  transition_matrix(1, 0) = 0.4; transition_matrix(1, 1) = 0.6; // From state 2 to 1 or 2
  int num_sim = 10;  // Run 10 simulations
//   
  List result = GroupRandomWalk(x, missing_positions, initial_prob, transition_matrix, num_sim);
//   
 return result;
 }

/***R
print(TestGroupRandomWalk())
*/
