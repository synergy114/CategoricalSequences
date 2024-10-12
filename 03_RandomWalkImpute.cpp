#include "03_RandomWalkImpute.h"
// [[Rcpp::depends(RcppArmadillo)]]


/////////////////////////////////////////////////////////////////////////////////////
//Function to find missing positions for the given sequence which groups the sequential NAs returns list of vectors. each vector is for 
// the cluster of NAs

// [[Rcpp::export]]
arma::vec SumWeightedVectors(arma::vec weights, std::vector<arma::vec> vectors) {
  arma::vec result = weights[0] * vectors[0];  // Initialize with first weighted vector
  
  for (std::size_t k = 1; k < weights.n_elem; k++) {
    result += weights[k] * vectors[k];
  }
  
  return result;
}

List FindMissingGroups(IntegerVector y) {
  int n = y.size();
  List missing_positions;
  IntegerVector current_slots;
  
  // Finding missing positions
  for (int i = 0; i < n; i++) {
    if (!IntegerVector::is_na(y[i])) {
      if (current_slots.size() > 0) {
        missing_positions.push_back(current_slots);
        current_slots = IntegerVector(); 
      }
    } else {
      current_slots.push_back(i); // Store 0-based index of missing value
    }
  }
  
  return missing_positions;
}
///////////////////////////////////////////////


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

// ////////////////////////////////////////////////

// [[Rcpp::export]]
List GroupRandomWalk(IntegerVector y, IntegerVector missing_positions, NumericVector initial_prob,  NumericMatrix transition_matrix, int num_sim) {
  int i = 0;
  int prev_state;
  int next_state=-1;
  int original_prev_state=-1;
  int group_size = missing_positions.size();
  
  List imputes_list;
  IntegerVector imputes_count;
  
  bool first_element_missing = (missing_positions[0] == 0);
  bool last_element_missing = (missing_positions[group_size - 1] == y.size() - 1);
  
  // Initialize previous state for the group when 1st element of the seq is not missing
  if (!first_element_missing) {
    original_prev_state = y[missing_positions[0] - 1];
  } 
  if (!last_element_missing) {
    next_state = y[missing_positions[group_size - 1] + 1];
  }
  
  while (i < num_sim) {
    IntegerVector sim_impute(group_size);
    
    if (first_element_missing) {
      sim_impute[0] = RMultinom(initial_prob);
    } else {
      sim_impute[0] = RMultinom(transition_matrix(original_prev_state - 1, _));
    }
    prev_state = sim_impute[0];
    int current_transition_state = prev_state - 1; // Convert to 0-based index
    
    for (int j = 1; j < group_size; j++) {
      int RW_single_value = RMultinom(transition_matrix(current_transition_state, _));
      sim_impute[j] = RW_single_value;
      current_transition_state = RW_single_value - 1; // Update to 0-based index
    }
    
    if (last_element_missing || (!last_element_missing && RMultinom(transition_matrix(current_transition_state, _)) == next_state)) {
      bool found_in_list = false;
      
      // Check if the newly generated imputation (sim_impute) already exists in the imputes_list.
      for (int k = 0; k < imputes_list.size(); k++) {
        if (Rcpp::is_true(all(as<IntegerVector>(imputes_list[k]) == sim_impute))) {
          imputes_count[k] += 1;
          found_in_list = true;
          break;
        }
      }
      
      // If the newly generated imputation (sim_impute) does not exist in the imputes_list, add it to the list and initialize its count in imputes_count to 1.
      if (!found_in_list) {
        imputes_list.push_back(clone(sim_impute)); // Clone the IntegerVector to avoid modifying sim_impute
        imputes_count.push_back(1);
      }
      
      i += 1;
    }
  }
  
  // Convert counts to probabilities
  NumericVector probabilities = as<NumericVector>(imputes_count) / num_sim;
  
  // Prepare the output list
  Rcpp::List GroupRW = Rcpp::List::create(
    Rcpp::Named("missing_positions") = missing_positions,
    Rcpp::Named("unique_imputes") = imputes_list,
    Rcpp::Named("probabilities") = probabilities
  );
  
  return GroupRW;
}


List ImputeOneSequence(IntegerVector y, NumericVector initial_prob, NumericMatrix transition_matrix, int num_sim) {
  List one_seq_imputes;
  List missing_groups = FindMissingGroups(y);
  // Debugging print for missing groups
  //Rcpp::Rcout << "Missing groups: " << missing_groups << std::endl;
  for (int i = 0; i < missing_groups.size(); i++) {
    IntegerVector group = missing_groups[i];
    //int group_len = group.size();
    
    List impute_result = GroupRandomWalk(y, group, initial_prob,  transition_matrix, num_sim);
    one_seq_imputes.push_back(impute_result);
  }
  
  return one_seq_imputes;
}


List PairJoints(List group1, List group2) {
  // Validate input
  if (!group1.containsElementNamed("missing_positions") ||
      !group1.containsElementNamed("unique_imputes") ||
      !group1.containsElementNamed("probabilities") ||
      !group2.containsElementNamed("missing_positions") ||
      !group2.containsElementNamed("unique_imputes") ||
      !group2.containsElementNamed("probabilities")) {
      Rcpp::stop("Input Lists must contain 'missing_positions', 'unique_imputes', and 'probabilities' elements");
  }
  
  // Combine the missing positions
  IntegerVector missing_positions = group1["missing_positions"];
  IntegerVector group2_missing_positions = group2["missing_positions"];
  for (int i = 0; i < group2_missing_positions.size(); i++) {
    missing_positions.push_back(group2_missing_positions[i]);
  }
  
  // Extract unique imputes and probabilities from the groups
  List group1_unique_imputes = group1["unique_imputes"];
  NumericVector group1_probabilities = group1["probabilities"];
  List group2_unique_imputes = group2["unique_imputes"];
  NumericVector group2_probabilities = group2["probabilities"];
  
  int group1_size = group1_unique_imputes.size();
  int group2_size = group2_unique_imputes.size();
  
  // Containers for the joint unique imputes and joint probabilities
  List pair_unique_imputes;
  NumericVector pair_joint_probabilities;
  
  
  // Loop over all combinations of unique imputes from both groups
  for (int i = 0; i < group1_size; i++) {
    IntegerVector imputes1 = group1_unique_imputes[i];
    for (int j = 0; j < group2_size; j++) {
      // Combine unique imputes from group1 and group2
      IntegerVector imputes2 = group2_unique_imputes[j];
      IntegerVector current_pair = clone(imputes1);
      for (int k = 0; k < imputes2.size(); k++) {
        current_pair.push_back(imputes2[k]);
      }
      
      pair_unique_imputes.push_back(current_pair);
      
      // Compute the joint probability
      double prob = group1_probabilities[i] * group2_probabilities[j];
      pair_joint_probabilities.push_back(prob);
      
      // Debugging output
      
    }
  }
  
  // Create the resulting list
  List joint_pairs = List::create(
    Rcpp::Named("missing_positions") = missing_positions,
    Rcpp::Named("unique_imputes") = pair_unique_imputes,
    Rcpp::Named("probabilities") = pair_joint_probabilities
  );
  
  return joint_pairs;
}
///for each sequence iteratively finds the cartesian product of imputations for different groups

// [[Rcpp::export]]
List SeqJointImpute(IntegerVector y, NumericVector initial_prob, NumericMatrix transition_matrix, int num_sim=100) {
  List group_imputes = ImputeOneSequence(y, initial_prob, transition_matrix, num_sim);
  int groups_n = group_imputes.size();
  
  // Debugging output
  // Rcpp::Rcout << "Missing Groups: " << groups_n << std::endl;
  // Rcpp::Rcout << "Group Imputes: " << group_imputes << std::endl;
  
  List paired_joins = group_imputes[0];
  
  for (int i = 1; i < groups_n; i++) {
    List next_group_imputes = group_imputes[i];
    paired_joins = PairJoints(paired_joins, next_group_imputes);
  }
  
  return paired_joins;
}

// Function to apply SeqJointImpute to a list of sequences
// [[Rcpp::export]]
List AllSeqJointImputeRCPP(List seq_list, List current_params, int num_sim=100) {
  int n = seq_list.size();
  List results(n);
  NumericVector tau = current_params["tau"];
  NumericMatrix alpha = current_params["alpha"];
  List Gamma = current_params["gamma"];
  
  NumericVector initial_prob=SumWeightedVectors(tau,alpha );
  NumericMatrix transition_matrix=SumWeightedMatrices(tau, Gamma);

  for (int i = 0; i < n; i++) {
    // convert sequence to IntegerVector
    IntegerVector seq = seq_list[i];
    List result = SeqJointImpute(seq, initial_prob, transition_matrix, num_sim);
    results[i] = result;
  }
  return results;
}

// [[Rcpp::export]]
List AllSeqJointImputePi(List seq_list, List current_params,List all_imputes, List Pi, int num_sim=100) {
  int n = seq_list.size();
  List results(n);
  
  NumericMatrix alpha = current_params["alpha"];
  List Gamma = current_params["gamma"];
  
  int K=alpha.nrow();
  
  for (int i = 0; i < n; i++) {
    debug_print("i=", i);
    IntegerVector seq = seq_list[i];
    arma::mat pi_i=Pi[i];
    List imputes_i=all_imputes[i];
    arma::vec Pu=imputes_i["probabilities"];
    arma::vec weights=pi_i*Pu;
    
    
    NumericMatrix transition_matrix=SumWeightedMatrices(weights,Gamma);
    NumericVector initial_prob=SumWeightedVectors(weights,alpha); 
    
    List result = SeqJointImpute(seq, initial_prob, transition_matrix, num_sim);
    results[i] = result;
  }
  return results;
}


/***R

AllSeqJointImputeR <- function(seq_list, current_params, num_sim) {
  # Apply SeqJointImpute to each sequence in the list using lapply
  tau=current_params$tau
  alpha=current_params$alpha
  Gamma=current_params$gamma
  
  initial_prob=SumWeightedVectors(tau,alpha );  #initial_prob in the context pf imputation
  transition_matrix=SumWeightedMatrices(tau, Gamma);
  
  results <- lapply(seq_list, function(seq) {
    SeqJointImpute(seq, initial_prob, transition_matrix, num_sim)
  })
  
  # Return the list of results
  return(results)
}

*/
