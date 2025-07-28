// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;

// [[Rcpp::export]]
Eigen::MatrixXd log_mat(const Eigen::MatrixXd& R) {
  Eigen::SelfAdjointEigenSolver<MatrixXd> solver(R);
  VectorXd eigvals = solver.eigenvalues();
  VectorXd log_eigvals = eigvals.array().log();
  return solver.eigenvectors() * log_eigvals.asDiagonal() * solver.eigenvectors().transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd exp_mat(const Eigen::MatrixXd& R) {
  Eigen::SelfAdjointEigenSolver<MatrixXd> solver(R);
  VectorXd eigvals = solver.eigenvalues();
  VectorXd exp_eigvals = eigvals.array().exp();
  return solver.eigenvectors() * exp_eigvals.asDiagonal() * solver.eigenvectors().transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd exp_solve_mat(const Eigen::MatrixXd& R) {
  Eigen::SelfAdjointEigenSolver<MatrixXd> solver(R);
  VectorXd eigvals = solver.eigenvalues();
  VectorXd exp_solve_eigvals = 1 / eigvals.array().exp();
  return solver.eigenvectors() * exp_solve_eigvals.asDiagonal() * solver.eigenvectors().transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd exp_sqrt_solve_mat(const Eigen::MatrixXd& R) {
  Eigen::SelfAdjointEigenSolver<MatrixXd> solver(R);
  VectorXd eigvals = solver.eigenvalues();
  VectorXd exp_sqrt_solve_eigvals = 1 / eigvals.array().exp().sqrt();
  return solver.eigenvectors() * exp_sqrt_solve_eigvals.asDiagonal() * solver.eigenvectors().transpose();
}
