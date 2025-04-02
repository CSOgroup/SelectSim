#include <RcppArmadillo.h>
using namespace Rcpp;


//' Compute  overlap stats
//' 
//' @param am The alteration matrix list
//' @return overlap the weight overlap between the pairs
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List rcpp_overlap(List v, int t=1){
    //Rcout << "The length of vec : " <<  v.size() << "\n";
    Rcpp::List x(v.size());
    for(int i=0; i<v.length(); ++i){
        arma::mat temp = as<arma::mat>(v[i]);
        x[i] = temp * temp.t(); 
    }
    return(x);
}

//' Compute weighted overlap stats
//'
//' @param am The alteration matrix list
//' @param W The weight matrix 
//' @return overlap the weight overlap between the pairs
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List rcpp_w_overlap(List v, NumericMatrix w, int t=1){
    //Rcout << "The length of vec : " <<  v.size() << "\n";
    Rcpp::List x(v.size());
    for(int i=0; i<v.length(); ++i){
        arma::mat temp = as<arma::mat>(v[i]);
        arma::mat W = as<arma::mat>(w);
        x[i] = (W%temp) * temp.t(); 
    }
    return(x);
}