make_family <- function(name, a1_fun, a2_fun, a4_fun, h_fun) {
  required_funs <- c(a1_fun, a2_fun, a4_fun, h_fun)
  if (!all(sapply(required_funs, is.function))) {
    stop("All components must be functions")
  }

  structure(
    list(
      family = name,
      a1_fun = a1_fun,
      a2_fun = a2_fun,
      a4_fun = a4_fun,
      h_fun = h_fun
    ),
    class = "gcr.family"
  )
}

gaussian_family <- function() {
  make_family(
    name = "gaussian",
    a1_fun = function(x) x,
    a2_fun = function(x) 1,
    a4_fun = function(x) 0,
    h_fun = function(x) x
  )
}

poisson_family <- function() {
  make_family(
    name = "Poisson",
    a1_fun = function(x) exp(x),
    a2_fun = function(x) exp(x),
    a4_fun = function(x) exp(x),
    h_fun = function(x) x
  )
}

binomial_family <- function() {
  a1_fun <- function(x) exp(x) / (1 + exp(x))

  make_family(
    name = "Binomial",
    a1_fun = a1_fun,
    a2_fun = function(x) exp(x) / (1 + exp(x)) ^ 2,
    a4_fun = function(x) {
      mu <- a1_fun(x)
      mu * (1 - mu) * ((1 - 2 * mu) ^ 2 - 2 * mu * (1 - mu))
    },
    h_fun = function(x) x
  )
}

calculate_beta <- function(Y, X, W, alpha, beta, phi, family) {

  if(family == "Gaussian") {
    family_type <- gaussian_family()
  }
  else if(family == "Poisson") {
    family_type <- poisson_family()
  }
  else if(family == "Binomial") {
    family_type <- binomial_family()
  }

  a1_fun <- family_type$a1_fun
  a2_fun <- family_type$a2_fun
  a4_fun <- family_type$a4_fun
  h_fun <- family_type$h_fun

  n <- length(Y)
  p <- length(beta)
  d <- length(alpha)

  S_1 <- numeric(p)
  H_1 <- matrix(0, p, p)
  for(i in 1:n) {
    m <- dim(X[[i]])[1]

    theta <- X[[i]] %*% beta
    mu <- a1_fun(h_fun((theta)))
    sigma <- phi * a2_fun(h_fun((theta)))
    A <- diag(as.vector(sigma), nrow = m)
    A_sqrt_solve <- diag(1 / sqrt(as.vector(sigma)), nrow = m)

    G <- recover_GZT(W[[i]] %*% alpha)
    R_solve <- exp_solve_mat(G)

    if(family == "Gaussian") {
      D <- X[[i]]
    }
    else if(family == "Poisson") {
      D <- diag(as.vector(mu), nrow = m) %*% X[[i]]
    }
    else if(family == "Binomial") {
      D <- diag(as.vector(mu * (1 - mu)), nrow = m) %*% X[[i]]
    }

    V_solve <- A_sqrt_solve %*% R_solve %*% A_sqrt_solve
    S_1 <- S_1 + t(D) %*% V_solve %*% (Y[[i]] - mu)
    H_1 <- H_1 + t(D) %*% V_solve %*% D
  }
  return(list(S = S_1, H = H_1))
}

calculate_alpha <- function(Y, X, W, alpha, beta, phi, accelerate, family) {

  if(family == "Gaussian") {
    family_type <- gaussian_family()
  }
  else if(family == "Poisson") {
    family_type <- poisson_family()
  }
  else if(family == "Binomial") {
    family_type <- binomial_family()
  }

  a1_fun <- family_type$a1_fun
  a2_fun <- family_type$a2_fun
  a4_fun <- family_type$a4_fun
  h_fun <- family_type$h_fun

  n <- length(Y)
  p <- length(beta)
  d <- length(alpha)

  S_2 <- numeric(d)
  H_2 <- matrix(0, d, d)
  for(i in 1:n) {
    m <- dim(X[[i]])[1]
    if(m == 1) next
    E_l <- generate_El(m)
    E_u <- generate_Eu(m)
    E_d <- generate_Ed(m)

    theta <- X[[i]] %*% beta
    mu <- a1_fun(h_fun((theta)))
    sigma <- phi * a2_fun(h_fun((theta)))
    A <- diag(as.vector(sigma), nrow = m)
    A_sqrt_solve <- diag(1 / sqrt(as.vector(sigma)), nrow = m)

    G <- recover_GZT(W[[i]] %*% alpha)
    eigen_G <- eigen(G, symmetric = TRUE)
    Q <- eigen_G$vectors
    lambda <- eigen_G$values
    Mid <- generate_Mid(lambda)
    Q_tensor <- kronecker(Q, Q)
    B <- Q_tensor %*% Mid %*% t(Q_tensor)
    # calculate partial
    partial <- E_l %*% (diag(m ^ 2) - B %*% t(E_d) %*% solve(E_d %*% B %*% t(E_d)) %*% E_d) %*% B %*% t(E_l + E_u)

    R <- exp_mat(G)
    R_solve <- exp_solve_mat(G)
    R_sqrt_solve <- exp_sqrt_solve_mat(G)

    R_hat <- A_sqrt_solve %*% (Y[[i]] - mu) %*% t(Y[[i]] - mu) %*% A_sqrt_solve
    eta <- vecl(R_solve %*% R_hat %*% R_solve - R_solve)
    S_2 <- S_2 + t(W[[i]]) %*% t(partial) %*% eta

    a_d4 <- a4_fun(theta)
    a_d2 <- a2_fun(theta)

    if(accelerate) {
      if(family == "Gaussian") {
        J <- generate_J(R_solve)
      }
      else if(family == "Poisson") {
        J <- generate_J_extra(R_solve, R_sqrt_solve, a_d4, a_d2, phi)
      }
      else if(family == "Binomial") {
        J <- generate_J_extra(R_solve, R_sqrt_solve, a_d4, a_d2, phi)
      }
    }
    else {
      J <- eta %*% t(eta)
    }

    H_2 <- H_2 + t(W[[i]]) %*% t(partial) %*% J %*% partial %*% W[[i]]
  }
  return(list(S = S_2, H = H_2))
}


calculate_alpha_S <- function(alpha, Y, X, W, beta, phi, family) {

  if(family == "Gaussian") {
    family_type <- gaussian_family()
  }
  else if(family == "Poisson") {
    family_type <- poisson_family()
  }
  else if(family == "Binomial") {
    family_type <- binomial_family()
  }

  a1_fun <- family_type$a1_fun
  a2_fun <- family_type$a2_fun
  a4_fun <- family_type$a4_fun
  h_fun <- family_type$h_fun

  n <- length(Y)
  p <- length(beta)
  d <- length(alpha)

  S_2 <- numeric(d)
  H_2 <- matrix(0, d, d)
  for(i in 1:n) {
    m <- dim(X[[i]])[1]
    if(m == 1) next
    E_l <- generate_El(m)
    E_u <- generate_Eu(m)
    E_d <- generate_Ed(m)

    theta <- X[[i]] %*% beta
    mu <- a1_fun(h_fun((theta)))
    sigma <- phi * a2_fun(h_fun((theta)))
    A <- diag(as.vector(sigma), nrow = m)
    A_sqrt_solve <- diag(1 / sqrt(as.vector(sigma)), nrow = m)

    G <- recover_GZT(W[[i]] %*% alpha)
    eigen_G <- eigen(G, symmetric = TRUE)
    Q <- eigen_G$vectors
    lambda <- eigen_G$values
    Mid <- generate_Mid(lambda)
    Q_tensor <- kronecker(Q, Q)
    B <- Q_tensor %*% Mid %*% t(Q_tensor)
    # calculate partial
    partial <- E_l %*% (diag(m ^ 2) - B %*% t(E_d) %*% solve(E_d %*% B %*% t(E_d)) %*% E_d) %*% B %*% t(E_l + E_u)

    R <- exp_mat(G)
    R_solve <- exp_solve_mat(G)
    R_sqrt_solve <- exp_sqrt_solve_mat(G)

    R_hat <- A_sqrt_solve %*% (Y[[i]] - mu) %*% t(Y[[i]] - mu) %*% A_sqrt_solve
    eta <- vecl(R_solve %*% R_hat %*% R_solve - R_solve)
    S_2 <- S_2 + t(W[[i]]) %*% t(partial) %*% eta
  }
  return(S_2)
}

calculate_alpha_H <- function(Y, X, W, alpha, beta, phi, family) {

  if(family == "Gaussian") {
    family_type <- gaussian_family()
  }
  else if(family == "Poisson") {
    family_type <- poisson_family()
  }
  else if(family == "Binomial") {
    family_type <- binomial_family()
  }

  a1_fun <- family_type$a1_fun
  a2_fun <- family_type$a2_fun
  a4_fun <- family_type$a4_fun
  h_fun <- family_type$h_fun

  n <- length(Y)
  p <- length(beta)
  d <- length(alpha)

  S_2 <- numeric(d)
  H_2 <- matrix(0, d, d)
  for(i in 1:n) {
    m <- dim(X[[i]])[1]
    if(m == 1) next
    E_l <- generate_El(m)
    E_u <- generate_Eu(m)
    E_d <- generate_Ed(m)

    theta <- X[[i]] %*% beta
    mu <- a1_fun(h_fun((theta)))
    sigma <- phi * a2_fun(h_fun((theta)))
    A <- diag(as.vector(sigma), nrow = m)
    A_sqrt_solve <- diag(1 / sqrt(as.vector(sigma)), nrow = m)

    G <- recover_GZT(W[[i]] %*% alpha)
    eigen_G <- eigen(G, symmetric = TRUE)
    Q <- eigen_G$vectors
    lambda <- eigen_G$values
    Mid <- generate_Mid(lambda)
    Q_tensor <- kronecker(Q, Q)
    B <- Q_tensor %*% Mid %*% t(Q_tensor)
    # calculate partial
    partial <- E_l %*% (diag(m ^ 2) - B %*% t(E_d) %*% solve(E_d %*% B %*% t(E_d)) %*% E_d) %*% B %*% t(E_l + E_u)

    R <- exp_mat(G)
    R_solve <- exp_solve_mat(G)
    R_sqrt_solve <- exp_sqrt_solve_mat(G)

    R_hat <- A_sqrt_solve %*% (Y[[i]] - mu) %*% t(Y[[i]] - mu) %*% A_sqrt_solve
    eta <- vecl(R_solve %*% R_hat %*% R_solve - R_solve)

    J <- eta %*% t(eta)

    H_2 <- H_2 + t(W[[i]]) %*% t(partial) %*% J %*% partial %*% W[[i]]
  }
  return(H_2)
}

calculate_hessian <- function(Y, X, W, alpha, beta, phi, family, eps = .Machine$double.eps ^ 0.5) {
  n <- length(Y)
  p <- length(beta)
  d <- length(alpha)

  I_d <- diag(d)

  perturbed_points <- matrix(0, nrow = d, ncol = 2 * d)
  for (j in 1:d) {
    perturbed_points[, j] <- alpha + eps * I_d[, j]
    perturbed_points[, j + d] <- alpha - eps * I_d[, j]
  }

  if (d == 1) {
    all_grads <- matrix(all_grads, nrow = 1)
  }

  grad_func <- function(theta) calculate_alpha_S(theta, Y = Y, X = X, W = W, beta = beta, phi = phi, family = family)
  all_grads <- apply(perturbed_points, 2, grad_func)

  grads_pos <- all_grads[, 1:d, drop = FALSE]
  grads_neg <- all_grads[, (d+1):(2*d), drop = FALSE]

  H <- (grads_pos - grads_neg) / (2 * eps)

  H_sym <- 0.5 * (H + t(H))

  return(H_sym)
}

est_phi <- function(Y, X, beta, family) {
  if(family == "Gaussian") {
    family_type <- gaussian_family()
  }
  else if(family == "Poisson") {
    family_type <- poisson_family()
  }
  else if(family == "Binomial") {
    family_type <- binomial_family()
  }

  a1_fun <- family_type$a1_fun
  a2_fun <- family_type$a2_fun
  a4_fun <- family_type$a4_fun
  h_fun <- family_type$h_fun

  y <- unlist(Y)
  x <- do.call(rbind, X)
  mu <- a1_fun(h_fun((x %*% beta)))
  var_fake <- a2_fun(h_fun(x %*% beta))
  res_fake <- (y - mu) / sqrt(var_fake)
  phi_hat <- sum(res_fake ^ 2) / (length(y) - dim(x)[2])
  return(phi_hat)
}
