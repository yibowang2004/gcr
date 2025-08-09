#' Generalized Correlation Regression
#'
#' Fits a Generalized Correlation Regression model,
#' supporting multiple distribution families and convergence criteria.
#'
#' @param Y Response variable as a list of length \code{n}, where each element is a numeric vector
#'          of responses for a cluster. For cluster \eqn{i}, \code{Y[[i]]} should be a vector of length \eqn{m_i}
#'          (number of observations in cluster \eqn{i}).
#'
#' @param X A list of length \code{n} containing moment-related covariate matrices for each cluster. Each element
#'          \code{X[[i]]} should be a numeric matrix of dimension \eqn{m_i \times p} where:
#'          \itemize{
#'            \item \eqn{m_i} is the number of observations in cluster \eqn{i}
#'            \item \eqn{p} is the number of covariates
#'          }
#'
#' @param W A list of length \code{n} containing correlation-related covariate matrices for each cluster. Each element
#'          \code{W[[i]]} should be a numeric matrix of dimension \eqn{m_i(1-m_i)/2 \times d} where:
#'          \itemize{
#'            \item \eqn{m_i} is the number of observations in cluster \eqn{i}
#'            \item \eqn{d} is the number of covariates
#'          }
#'
#' @param alpha_init Initial value for alpha
#' @param beta_init Initial value for beta
#' @param phi_init Initial value for phi
#' @param phi.include Logical indicating whether phi is known
#' @param family Response distribution family: "Gaussian", "Poisson", or "Binomial"
#' @param lambda Step size
#' @param max_iter_1 Maximum iterations for outer optimization loop
#' @param max_iter_2 Maximum iterations for inner optimization loop
#' @param tol Convergence tolerance threshold
#' @param criteria Convergence criteria: "sum" (sum of absolute changes) or
#'                "avg" (mean absolute change)
#' @param independent Logical indicating whether use independent correlation structure
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{alpha}{Estimated alpha}
#'   \item{pvalue_alpha}{p-value for estimated alpha}
#'   \item{beta}{Estimated beta}
#'   \item{pvalue_beta}{p-value for estimated beta}
#'   \item{phi}{Estimated phi}
#'   \item{H_1}{Estimated H_1}
#'   \item{H_2}{Estimated H_2}
#'   \item{Hessian}{Estimated Hessian matrix}
#' }
#'
#' @export
gcr <- function(Y, X, W,
                alpha_init, beta_init, phi_init, phi.include = TRUE,
                family, lambda = 1,
                max_iter_1 = 100, max_iter_2 = 100, tol = 1e-6, criteria = "sum",
                independent = FALSE) {

  n <- length(Y)
  p <- dim(X[[1]])[2]
  d <- dim(W[[1]])[2]

  beta_0 <- beta_init
  beta_1 <- beta_0 + tol

  alpha_0 <- alpha_init
  alpha_1 <- alpha_0 + tol

  if(family == "Poisson" || family == "Binomial") {
    phi.include <- FALSE
    phi_init <- 1
  }

  phi_0 <- phi_init
  if(phi.include) phi_1 <- phi_0 + tol
  else phi_1 <- phi_init

  step <- 0
  while(step < max_iter_1) {
    if(step != 0) {
      beta_0 <- beta_1
      alpha_0 <- alpha_1
      phi_0 <- phi_1
    }

    # update phi
    if(phi.include) phi_1 <- est_phi(Y, X, beta_0, family)
    else phi_1 <- phi_init

    # update alpha
    alpha_10 <- alpha_0
    alpha_11 <- alpha_10 + tol
    step_a <- 0
    while(step_a < max_iter_2) {
      if(independent) {
        alpha_11 <- rep(0, d)
        break
      }

      if(step_a != 0) alpha_10 <- alpha_11

      res_alpha <- calculate_alpha(Y, X, W, alpha_10, beta_0, phi_1, family)
      S_2 <- res_alpha$S
      H_2 <- res_alpha$H
      alpha_11 <- alpha_10 + lambda * solve(H_2) %*% S_2

      step_a <- step_a + 1

      # check convergence
      if(criteria == "sum" && sum(abs(alpha_11 - alpha_10)) <= tol) {
        break
      }
      else if(criteria == "avg" && mean(abs(alpha_11 - alpha_10)) <= tol) {
        break
      }
    }
    alpha_1 <- alpha_11

    # update beta
    res_beta <- calculate_beta(Y, X, W, alpha_1, beta_0, phi_1, family)
    S_1 <- res_beta$S
    H_1 <- res_beta$H
    beta_1 <- beta_0 + solve(H_1) %*% S_1

    step <- step + 1

    # check convergence
    if(criteria == "sum" && sum(abs(c(beta_1 - beta_0, alpha_1 - alpha_0, phi_1 - phi_0))) <= tol) {
      break
    }
    else if(criteria == "avg") {
      if(phi.include && mean(abs(c(beta_1 - beta_0, alpha_1 - alpha_0, phi_1 - phi_0))) <= tol) {
        break
      }
      if(!phi.include && mean(abs(c(beta_1 - beta_0, alpha_1 - alpha_0))) <= tol) {
        break
      }
    }
  }

  H_1 <- calculate_beta(Y, X, W, alpha_1, beta_1, phi_1, family)$H
  if(!independent) {
    H_2 <- calculate_alpha(Y, X, W, alpha_1, beta_1, phi_1, family)$H
    Hessian <- calculate_hessian(Y, X, W, alpha_1, beta_1, phi_1, family, 1e-14)
    Hessian_solve <- solve(Hessian)
    sigma_alpha <- sqrt(diag(Hessian_solve %*% H2 %*% Hessian_solve))
    std_alpha <- alpha_1 / sigma_alpha
    pvalue_alpha <- ifelse(std_alpha > 0, 2 - 2 * pnorm(std_alpha), 2 * pnorm(std_alpha))
  }
  else {
    H_2 <- NULL
    Hessian <- NULL
    pvalue_alpha <- NULL
  }

  sigma_beta <- sqrt(diag(solve(H_1)))
  std_beta <- beta_1 / sigma_beta
  pvalue_beta <- ifelse(std_beta > 0, 2 - 2 * pnorm(std_beta), 2 * pnorm(std_beta))

  return(list(alpha = alpha_1, pvalue_alpha = pvalue_alpha,
              beta = beta_1, pvalue_beta = pvalue_beta,
              phi = phi_1,
              H_1 = H_1, H_2 = H_2, Hessian = Hessian))
}
