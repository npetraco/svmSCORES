#include <Rcpp.h>
//#include <boost/foreach.hpp>
#include<iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP empiricalPvalues(NumericVector nullvec, NumericVector vec) { //This seems to do fine. ~5x faster than in R alone.
  
  int n = nullvec.size();
  int m = vec.size();
  Rcpp::NumericVector pvalVec(m);
  
  for(int i=0; i<m; i++) {
    int count = 0;
    for(int j=0; j<n; j++) {
      count = count + (nullvec[j]>=vec[i]);
      //cout<<count<<endl;
    }
    double val = ((double)count)/((double)n);
    //Rcout << "p-value #: " << i <<  " = " << val << endl;
    pvalVec[i] = val;
  }
  
  return pvalVec;
  
}

// [[Rcpp::export]]
SEXP compPvals2(NumericVector nullvec, NumericVector vec) { //Try above except with iterators. Runs about as fast.
  
  int n = nullvec.size();
  int m = vec.size();
  Rcpp::NumericVector pvalVec(m);
  
  typedef Rcpp::NumericVector::iterator vec_iterator;
  vec_iterator it_nv = nullvec.begin();
  vec_iterator it_v = vec.begin();
  vec_iterator it_pv = pvalVec.begin();
  
  for(int i=0; i<m; i++) {
    int count = 0;
    for(int j=0; j<n; j++) {
      count = count + (it_nv[j]>=it_v[i]);
      //cout<<count<<endl;
    }
    double val = ((double)count)/((double)n);
    //cout<<val<<endl;
    it_pv[i] = val;
  }
  
  return pvalVec;
  
}

// [[Rcpp::export]]
SEXP compPvals3(NumericVector nullvec, NumericVector vec) { //A fancyer version using iterators. Actually runs much slower!!!!
  
  int n = nullvec.size();
  int m = vec.size();
  Rcpp::NumericVector pvalVec(m);
  
  typedef Rcpp::NumericVector::iterator vec_iterator;
  //vec_iterator it_nv = nullvec.begin();
  //vec_iterator it_v = vec.begin();
  //vec_iterator it_pv = pvalVec.begin();
  
  for(vec_iterator it_v = vec.begin(), it_pv = pvalVec.begin(); it_v != vec.end(); ++it_v, ++it_pv) {
    int count = 0;
    for(vec_iterator it_nv = nullvec.begin(); it_nv != nullvec.end(); ++it_nv) {
      count = count + ( *it_nv >= *it_v );
      //cout<<count<<endl;
    }
    double val = ((double)count)/((double)n);
    //cout<<val<<endl;
    *it_pv = val;
  }
  
  return pvalVec;
  
}