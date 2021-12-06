#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat rdist_cpp(mat X, Rcpp::Function crossproduct) {
  // Ref: https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
  colvec Xn =  sum(square(X),1); // dim = 1 means row sum
  mat C = -2 * (Rcpp::as<mat>(crossproduct(X.t()))); // fast C XX^T in R base
  C.each_col() += Xn;
  C.each_row() += Xn.t();
  return sqrt(C); 
}