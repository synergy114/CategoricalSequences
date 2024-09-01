#ifndef CompleteImputeV3_h
#define CompleteImputeV3_h
#include <Rcpp.h>
using namespace Rcpp;

// Function declarations
List FindMissingGroups(IntegerVector x);
int RMultinom(NumericVector probs);
NumericVector EmpiricProb(List sequences, int categories);
List GroupRandomWalk(IntegerVector x, NumericVector emp_prob, IntegerVector missing_positions, NumericMatrix transition_matrix, int num_sim);
List ImputeOneSequence(IntegerVector x, NumericVector emp_prob, NumericMatrix transition_matrix, int num_sim);
List PairJoints(List group1, List group2);
List SeqJointImpute(IntegerVector x, NumericVector emp_prob, NumericMatrix transition_matrix, int num_sim);
List AllSeqJointImputeRCPP(List seq_list, NumericVector emp_prob, NumericMatrix transition_matrix, int num_sim);

#endif