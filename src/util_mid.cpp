#include <Rcpp.h>
  
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix generate_J(NumericMatrix R_solve) {
	int m = R_solve.nrow();
	NumericMatrix J(m * (m - 1) / 2, m * (m - 1) / 2);
	for(int j = 0; j < m - 1; ++j) {
	  for(int i = j + 1; i < m; ++i) {
	    for(int l = 0; l < m - 1; ++l) {
	      for(int k = l + 1; k < m; ++k) {
	        J((2 * m - j - 1) * j / 2 + i - j - 1, (2 * m - l - 1) * l / 2 + k - l - 1) = R_solve(i, k) * R_solve(j, l) + R_solve(i, l) * R_solve(j, k);
	      }
	    }
	  }
	}
	return J;
}

// [[Rcpp::export]]
NumericMatrix generate_J_extra(NumericMatrix R_solve, NumericMatrix R_sqrt_solve, NumericVector a_d4, NumericVector a_d2, double phi) {
	int m = R_solve.nrow();
	NumericMatrix J(m * (m - 1) / 2, m * (m - 1) / 2);
	for(int j = 0; j < m - 1; ++j) {
	  for(int i = j + 1; i < m; ++i) {
	    for(int l = 0; l < m - 1; ++l) {
	      for(int k = l + 1; k < m; ++k) {
	        double sum = 0;
	        for(int t = 0; t < m; ++t) {
	          sum += (a_d4[t] / (a_d2[t] * a_d2[t])) * R_sqrt_solve(t, i) * R_sqrt_solve(t, j) * R_sqrt_solve(t, l) * R_sqrt_solve(t, k);
	        }
	        sum *= phi;
	        J((2 * m - j - 1) * j / 2 + i - j - 1, (2 * m - l - 1) * l / 2 + k - l - 1) = R_solve(i, k) * R_solve(j, l) + R_solve(i, l) * R_solve(j, k) + sum;
	      }
	    }
	  }
	}
	return J;
}

// [[Rcpp::export]]
NumericMatrix generate_Mid(NumericVector lambda) {
	int m = lambda.size();
	int n = m * m;
	NumericMatrix Mid(n, n);
	const double EPSILON = 1e-10;
	for(int j = 0; j < m; ++j) {
	  for(int k = 0; k < m; ++k) {
	    int idx = j * m + k;
	    if(std::fabs(lambda[j] - lambda[k]) < EPSILON) {
	      Mid(idx, idx) = std::exp(lambda[j]);
	    } 
	    else {
	      double numerator = std::exp(lambda[j]) - std::exp(lambda[k]);
	      double denominator = lambda[j] - lambda[k];
	      Mid(idx, idx) = numerator / denominator;
	    }
	  }
	}
	return Mid;
}
