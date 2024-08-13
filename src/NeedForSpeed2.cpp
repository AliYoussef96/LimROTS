#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Sort array b based on the array a (decreasingly)
void sort2_2( double* a, double* b, int n ) {
  // Create vector of pairs.
  static vector< pair<double, double> > pairs;
  pairs.reserve( n );
  for( int i = 0; i < n; ++i)
    pairs.push_back( pair<double, double>(a[i], b[i]) );
  // Sort the pairs (inc). By default pairs are sorted by the first value and
  // in the case of a tie, the second values are used.
  sort( pairs.begin(), pairs.end());
  // Split the pairs back into the original vectors (dec).
  for( int i = 0; i < n; ++i ) {
    a[ n-1-i ] = pairs[i].first;
    b[ n-1-i ] = pairs[i].second;
  }
	// Empty the vector
	pairs.clear();
}

// Array printing routine for debugging.
template<class T> void printArray( T *a, int a_len ){
  for( int i = 0; i < a_len; ++i )
    cout << a[i] << " ";
  cout << endl;
}

// Calculate the overlap
void calculateOverlap_2( double *r1, double *r2, int r_len, IntegerVector N_ovlp, int N_len_ovlp,
		       int b, int B_ovlp, NumericVector &result ){
  // Copy r2 and sort the copy.
  double *r3 = new double[ r_len ];
	for( int i = 0; i < r_len; ++i )
    r3[i] = r2[i];
  sort(r3, r3 + r_len );
  reverse(r3, r3 + r_len );

  // Sort r2 by r1
  sort2_2( r1, r2, r_len );

  // Calculate the overlap
  for( int i = 0; i < N_len_ovlp; ++i ){
    static double sum = 0;
    for( int j = 0; j <= ( N_ovlp[i] - 1 ); ++j )
      sum += ( r2[j] >= r3[ N_ovlp[i] - 1 ] );
    result[ (b-1) + i*B_ovlp ] = sum / N_ovlp[i];
    sum = 0;
  }
	delete[] r3;
}

//Loop2

// [[Rcpp::export]]

List NeedForSpeed2( SEXP D, SEXP pD, SEXP nrow,
		     SEXP N,  SEXP N_len, SEXP B,
		     SEXP overlaps, SEXP overlaps_P){

	NumericVector D_ovlp(D);
	NumericVector pD_ovlp(pD);
	int nrow_ovlp = Rcpp::as<int>(nrow);
	IntegerVector N_ovlp(N);
	int N_len_ovlp = Rcpp::as<int>(N_len);
	int B_ovlp = Rcpp::as<int>(B);
	NumericVector overlaps_ovlp(overlaps);
	NumericVector overlaps_P_ovlp(overlaps_P);

    double *res1 = new double[ nrow_ovlp ];
    double *res2 = new double[ nrow_ovlp ];
    double *pres1 = new double[ nrow_ovlp ];
    double *pres2 = new double[ nrow_ovlp ];

    for( int b = 1; b <= B_ovlp; ++b ){
      for( int i = 0; i < nrow_ovlp; ++i ){
				res1[i] = fabs( D_ovlp[ ( b - 1 ) * (nrow_ovlp) + i ] );
				res2[i] = fabs( D_ovlp[ ( b + (B_ovlp) - 1 ) * (nrow_ovlp) + i ] );
				pres1[i] = fabs( pD_ovlp[ ( b - 1 ) * (nrow_ovlp) + i ] );
				pres2[i] = fabs( pD_ovlp[ ( b + (B_ovlp) - 1 ) * (nrow_ovlp) + i ] );
      }
     calculateOverlap_2( res1, res2, nrow_ovlp, N_ovlp, N_len_ovlp, b, B_ovlp, overlaps_ovlp );
     calculateOverlap_2( pres1, pres2, nrow_ovlp, N_ovlp, N_len_ovlp, b, B_ovlp, overlaps_P_ovlp );

    }
    delete[] res1; delete[] res2; delete[] pres1; delete[] pres2;
    
    return Rcpp::List::create(
        Rcpp::Named("overlaps") = overlaps_ovlp ,
        Rcpp::Named("overlaps_P") = overlaps_P_ovlp ) ;
	
  }
