#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List structureData(NumericVector values, List imputations, NumericVector probabilities, NumericMatrix transitions, List responsibilities) {
  
  // Create a list to hold the structured data
  List structuredData = List::create(
    Named("values") = values,
    Named("imputations") = imputations,
    Named("probabilities") = probabilities,
    Named("responsibilities") = responsibilities
  );
  
  return structuredData;
}

/*** R
# Example data
values <- c(1, 2, NA, 4)
imputations <- list(c(1, 2), c(2, 2), c(1, 1))
probabilities <- c(0.6, 0.3, 0.1)

# Responsibilities (k = 2 for example)
responsibilities <- list(
  c(0.4, 0.6, 0.5),
  c(0.6, 0.4, 0.5)
)

# Call the Rcpp function
structuredData <- structureData(values, imputations, probabilities, transitions, responsibilities)
print(structuredData)
*/
