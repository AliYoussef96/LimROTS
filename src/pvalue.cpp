#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector pvalue(SEXP a, SEXP b) {
  NumericVector observed(a);
  NumericVector permuted(b);
  NumericVector pvalues(observed.length());
  
  int j = 0;
  for (int i=0; i<observed.length(); i++) {
    while(permuted[j]>=observed[i] && j<permuted.length()) {
      j++;
    }
    pvalues[i] = double(j) / double(permuted.length());
  }

  return pvalues;
}
