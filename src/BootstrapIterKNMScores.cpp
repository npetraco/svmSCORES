#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector BootstrapIterKNMScores(NumericMatrix iterProbMat) {
  
  //"as" the pointer sent in into a Rcpp Xptr 
  Rcpp::NumericVector nullVecChunk( iterProbMat.nrow() );
  //Rcout << nullVecChunk << endl;
  
  for(int i = 0; i<iterProbMat.nrow(); i++){
    //Rcout << nullVecChunk[i] << endl;
    //NumericVector zz1 = xx(_,1);
    NumericVector tmpProbVec = iterProbMat(i,_);
    //lbls[-mod.bsidx][j])] //DO THIS OUTSIDE THE FUNCTION??????
    //tmp.prob.vec<-sample(tmp.prob.vec[-as.numeric(lbls[-mod.bsidx][j])],1)  #Just randomly grab one of the KNM scores
  }
  Rcout << endl;
  
  return nullVecChunk; //IMPLEMENT SMILE SIDE ERROR HANDLEING FOR PRODUCTION!!!!!!!!
    
}
