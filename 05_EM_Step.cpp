#include "05_EM_Step.h"


// Debugging helper function
void debug_print(std::string msg, int value) {
  Rcpp::Rcout << msg << ": " << value << std::endl;
}

// // [[Rcpp::export]]
// NumericVector SumWeightedVectors(NumericVector weights, NumericMatrix matrix){
//   int n=weights.size();
//   int m = matrix.ncol();  // Get the number of columns in the matrix
//   
//   NumericVector weighted_vector(m);
//   for (int i=0; i<n; i++){
//     weighted_vector += weights[i] * matrix(i, _);
//   }
//   return weighted_vector;
// }

// [[Rcpp::export]]
arma::vec SumWeightedVectors(arma::vec weights, std::vector<arma::vec> vectors) {
  arma::vec result = weights[0] * vectors[0];  // Initialize with first weighted vector
  
  for (std::size_t k = 1; k < weights.n_elem; k++) {
    result += weights[k] * vectors[k];
  }
  
  return result;
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


// Constants
const double doubleMin = 1e-10;  // Small threshold to avoid zeros
const double doubleMax = 1e+10;  // Large threshold for log operation


NumericMatrix Pi_i(IntegerVector y_i1, int nCategories, List current_params, List x_i) {
  
  NumericVector tau=current_params["tau"];
  NumericMatrix initial_prob=current_params["alpha"]; 
  List transition_matrices=current_params["gamma"]; 
  
  int J = nCategories; // Number of states
  int K = tau.size();  // Number of clusters
  int Li = x_i.size(); // Number of imputations for the obs i
  NumericMatrix pi_i(K, Li);
  
  // Check that initial_prob dimensions match expected K and categories
  if (initial_prob.nrow() != K || initial_prob.ncol() != J) {
    Rcpp::stop("Error: initial_prob dimensions do not match expected K and J.");
  }
  
  for (int u = 0; u < Li; u++) {
    IntegerVector X = x_i[u];
    
    for (int k = 0; k < K; k++) {
      NumericMatrix gamma_k = transition_matrices[k];
      double inv_pi_iuk = 0.0;
      
      for (int k_prime = 0; k_prime < K; k_prime++) {
        NumericMatrix gamma_k_prime = transition_matrices[k_prime];
        
        // Check that y_i1[u] is a valid column index for initial_prob
        if (y_i1[u] < 0 || (y_i1[u]-1) >= initial_prob.ncol()) {
          Rcpp::stop("Error: y_i1[u] is out of bounds");
        }
        
        double inner_pi_iuk = std::max(-doubleMax, 
                                       log(std::max(doubleMin, tau[k_prime]) / 
                                         std::max(doubleMin, tau[k]))) +
                                         std::max(-doubleMax, 
                                                  log(std::max(doubleMin, initial_prob(k_prime, y_i1[u]-1)) / 
                                                    std::max(doubleMin, initial_prob(k, y_i1[u]-1))));
        
        for (int j = 0; j < J; j++) {
          for (int j_prime = 0; j_prime < J; j_prime++) {
            if (j >= X.size() || j_prime >= X.size()) {
              Rcpp::stop("Error: X index out of bounds");
            }
            
            inner_pi_iuk += X(j, j_prime) * 
              std::max(-doubleMax, 
                       log(std::max(doubleMin, gamma_k_prime(j, j_prime)) / 
                         std::max(doubleMin, gamma_k(j, j_prime))));
          }
        }
        
        inv_pi_iuk += exp(inner_pi_iuk);
      }
      
      pi_i(k, u) = std::max(1.0 / inv_pi_iuk, doubleMin);
    }
  }
  
  return pi_i;
}

// [[Rcpp::export]]
List E_Step(int nCategories, List XY1, List current_params){
  
  List Pi;
  //int nCategories=initial_prob.ncol();
    
  for (int i = 0; i < XY1.size(); i++) {
    NumericMatrix pi_i;
    List XY1_i = XY1[i];
    IntegerVector y_i1 = XY1_i["Y_i1"];
    List X_i = XY1_i["x_i"];
    pi_i = Pi_i(y_i1, nCategories,current_params, X_i);
    Pi.push_back(pi_i); 
  }
  return Pi;
}

// [[Rcpp::export]]
List M_Step(int K,int nCategories, List allImputes, List Pi, List XY1){
  int J=nCategories;
  int n=Pi.size();
  
  NumericVector tau(K);
  NumericMatrix alpha(K, J);
  List transition_matrces;
  
  
  for (int k=0; k<K; k++){
    tau[k]=0.0;
    NumericMatrix gamma(J,J);
    
    for (int j=0; j<J; j++) {
      alpha(k,j)=0;
      
      for (int j_prime=0; j_prime<J; j_prime++){
        //gamma(j, j_prime)=0;
        double gamma_denom=0.0;
        double gamma_num=0.0;
        
        for (int i=0; i<n;i++ ){
          List imputes_i=allImputes[i];
          NumericVector probs=imputes_i["probabilities"];
          NumericMatrix Pi_i=Pi[i];
          List XY1_i = XY1[i];
          IntegerVector y_i1=XY1_i["Y_i1"];
          List X_i=XY1_i["x_i"];
          
          int Li=probs.size();
          
          for (int u=0; u<Li; u++){
            IntegerMatrix X_iu=X_i[u];
            double PPiku=Pi_i(k,u)*probs[u];
            
            tau[k]+=PPiku;
            
            if((y_i1[u]-1)==j) { //I(y_iu1=j)
              alpha(k,j)+=PPiku;
            }
            gamma_num+= PPiku*X_iu(j, j_prime); //
            gamma_denom+= PPiku*rowSum(X_iu, j_prime);
          }
        }
        gamma(j,j_prime)=gamma_num/gamma_denom;
      }
      alpha(k,j)=alpha(k,j)/tau[k];
    }
    tau[k]=tau[k]/n;
    transition_matrces.push_back(gamma);
  }
  
  return List::create(
    Named("tau") = tau,
    Named("alpha") = alpha,
    Named("gamma") = transition_matrces);}

// [[Rcpp::export]]
double Q_function(int K, int nCategories, List current_params, List allImputes, List Pi, List XY1) {
  
  double doubleMin = 1e-10; // Small value to prevent underflow
  double doubleMax = 1e10;  // Large value to prevent overflow in logs
  
  NumericVector tau = current_params["tau"];
  NumericMatrix alpha = current_params["alpha"];
  List transition_matrices = current_params["gamma"];
  
  int n = Pi.size();
  int J = nCategories;
  double Q = 0.0;

  for (int i = 0; i < n; i++) {
    List XY1_i = XY1[i];
    IntegerVector y_i1 = XY1_i["Y_i1"];
    List X_i = XY1_i["x_i"];
    
    List imputes_i = allImputes[i];
    NumericVector Pui = imputes_i["probabilities"];
    
    int Li = Pui.size();
    NumericMatrix Pi_i = Pi[i];

    double InnerQu = 0.0;
    
    for (int u = 0; u < Li; u++) {
      IntegerMatrix X_iu = X_i[u];
      NumericVector Pi_iu = Pi_i(_,u);
      

      double InnerQk = 0.0;
      
      for (int k = 0; k < K; k++) {
        NumericMatrix gamma_k = transition_matrices[k];

        double InnerQ = log(std::max(doubleMin, tau[k])); // Prevent log(0)
        
        for (int j = 0; j < J; j++) {
          if ((y_i1[u]-1) == j) {
            InnerQ += log(std::max(doubleMin, alpha(k, j))); // Prevent log(0)
          }
          for (int j_prime = 0; j_prime < J; j_prime++) {
            InnerQ += X_iu(j, j_prime) * log(std::max(doubleMin, gamma_k(j, j_prime))); // Prevent log(0)
          }
        }
        InnerQk += std::max(doubleMin, std::min(doubleMax, Pi_iu[k])) * InnerQ; // Ensure Pi_iu[k] is not too small
      }
      InnerQu += std::max(doubleMin, std::min(doubleMax, Pui[u])) * InnerQk; // Prevent Pui[u] from being too small or too large
    }
    Q += InnerQu;
  }
  
  return Q;
}