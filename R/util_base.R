log_mat <- function(X) {
  eigen_X <- eigen(X, symmetric = TRUE)
  return(eigen_X$vectors %*% diag(log(eigen_X$values), nrow = dim(X)[1]) %*% t(eigen_X$vectors))
}

exp_mat <- function(X) {
  eigen_X <- eigen(X, symmetric = TRUE)
  return(eigen_X$vectors %*% diag(exp(eigen_X$values), nrow = dim(X)[1]) %*% t(eigen_X$vectors))
}

exp_solve_mat <- function(X) {
  eigen_X <- eigen(X, symmetric = TRUE)
  return(eigen_X$vectors %*% diag(1 / exp(eigen_X$values), nrow = dim(X)[1]) %*% t(eigen_X$vectors))
}

exp_sqrt_solve_mat <- function(X) {
  eigen_X <- eigen(X, symmetric = TRUE)
  return(eigen_X$vectors %*% diag(1 / sqrt(exp(eigen_X$values)), nrow = dim(X)[1]) %*% t(eigen_X$vectors))
}

vecl <- function(X) {
  return(X[lower.tri(X, diag = FALSE)])
}

lower_to_full <- function(x) {
  n <- (1 + sqrt(1 + 8 * length(x))) / 2
  G <- matrix(0, n, n)
  lower <- lower.tri(G, diag = FALSE)
  G[lower] <- x
  G <- t(G)
  G[lower] <- x
  G <- t(G)
  return(G)
}

GZT <- function(X) {
  return(vecl(log_mat(X)))
}

recover_GZT <- function(x, tol = 1e-14, max_iter = 100) {
  G <- lower_to_full(x)
  n <- dim(G)[1]

  x_0 <- numeric(n)
  diag(G) <- x_0
  x_1 <- x_0 - log(diag(exp_mat(G)))

  step <- 0
  while(sum(abs(x_1 - x_0)) > tol && step < max_iter) {
    x_0 <- x_1
    diag(G) <- x_0
    x_1 <- x_0 - log(diag(exp_mat(G)))
    step <- step + 1
  }

  diag(G) <- x_1
  return(G)
}

generate_El <- function(n){
  #n: dimension
  E <- matrix(0, n * (n - 1) / 2, n * n)
  k <- 1
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      E[k, n * (i - 1) + j] <- 1
      k <- k + 1
    }
  }
  return(E)
}

generate_Eu <- function(n){
  #n: dimension
  E <- matrix(0, n * (n - 1) / 2, n * n)
  k <- 1
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      E[k, n * (j - 1) + i] <- 1
      k <- k + 1
    }
  }
  return(E)
}

generate_Ed <- function(n){
  #n: dimension
  E=matrix(0, n, n * n)
  for(i in 1:n){
    E[i, n * (i - 1) + i] <- 1
  }
  return(E)
}
