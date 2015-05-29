#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector BootstrapIterKNMScores(NumericMatrix iterProbMat) {
  
  //"as" the pointer sent in into a Rcpp Xptr 
  Rcpp::NumericVector nullVecChunk( iterProbMat.nrow() );
  //Rcout << nullVecChunk << endl;
  
  Rcout << "With all the index dropping in the loop, relized that many new mallocs/reallocs will be required. Keep skeleton around for later though." << endl;
  
  //Index dropping can be done with std::set_difference, but two new index arrays need to be malloced/freed per index set....
  
  for(int i = 0; i<iterProbMat.nrow(); i++){
    //Rcout << nullVecChunk[i] << endl;
    //NumericVector zz1 = xx(_,1);
    NumericVector tmpProbVec = iterProbMat(i,_);
    //lbls[-mod.bsidx][j])] //DO THIS OUTSIDE THE FUNCTION?????? 
    //Still though, have to do the second index drop on what remains of lbls in here
    
    //R code. Is there a better way????????:
    //tmp.prob.vec<-sample(tmp.prob.vec[-as.numeric(lbls[-mod.bsidx][j])],1)  #Just randomly grab one of the KNM scores
  }
  //Rcout << endl;
  
  return nullVecChunk;
    
}
