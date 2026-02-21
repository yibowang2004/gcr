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

stacked_sandwich <- function(
    Y, X, W,
    alpha, beta,
    phi,                      # if phi_known=TRUE: known value; else phi_hat
    phi_known = FALSE,
    eps = 1e-6,               # used ONLY for A_{alpha,alpha} finite diff
    scale_eps = TRUE,
    force_sym = TRUE,
    use_minus = TRUE,         # TRUE: A_bb = -H1/n (recommended)
    family
) {
  if(family == "Gaussian") {
    family_type <- gaussian_family()
    # Gaussian: identity link, Var(Y)=phi*1
    dmu_deta <- function(eta, mu) {
      rep(1, length(mu))
    }

    Vprime_from_mu <- function(mu) {
      rep(0, length(mu))   # a2(mu)=1 => derivative 0
    }
  }
  else if(family == "Poisson") {
    family_type <- poisson_family()
    # Poisson: log link, Var(Y)=phi*mu  (phi=1 for true Poisson, can estimate for quasi-Poisson)
    dmu_deta <- function(eta, mu) {
      as.vector(mu)        # = exp(eta)
    }

    Vprime_from_mu <- function(mu) {
      rep(1, length(mu))   # a2(mu)=mu => derivative 1
    }
  }
  else if(family == "Binomial") {
    family_type <- binomial_family()
    # Binomial (Bernoulli): logit link, Var(Y)=phi*mu*(1-mu) (phi usually fixed at 1)
    dmu_deta <- function(eta, mu) {
      as.vector(mu * (1 - mu))    # derivative of logistic
    }

    Vprime_from_mu <- function(mu) {
      as.vector(1 - 2 * mu)       # derivative of mu*(1-mu)
    }
  }

  a1_fun <- family_type$a1_fun
  a2_fun <- family_type$a2_fun
  a4_fun <- family_type$a4_fun
  h_fun <- family_type$h_fun

  n <- length(Y)
  p <- ncol(X[[1]])
  d <- ncol(W[[1]])
  N <- sum(vapply(Y, length, 1L))
  if (is.null(phi) || length(phi) != 1) stop("Please provide scalar phi (known or estimated).")

  # ---------------- helper: psi3i Version A (df-corrected, cluster-split) ----------------
  psi3i_versionA <- function(Yi, Xi, beta, phi, p, N) {
    eta <- Xi %*% beta
    z <- h_fun(eta)
    mu <- a1_fun(z)
    var_base <- a2_fun(z)     # WITHOUT phi (Pearson-style)
    r <- (Yi - mu) / sqrt(as.vector(var_base))

    m <- length(Yi)
    ci <- m - (p * m / N)
    sum(r^2) - phi * ci
  }

  # ---------------- build partial and L = t(Wi)%*%t(partial) (depends on alpha only) ----------------
  build_partial_L <- function(Wi, G, m) {
    E_l <- generate_El(m)
    E_u <- generate_Eu(m)
    E_d <- generate_Ed(m)

    eigen_G <- eigen(G, symmetric = TRUE)
    Q <- eigen_G$vectors
    lambda <- eigen_G$values
    Mid <- generate_Mid(lambda)

    Q_tensor <- kronecker(Q, Q)
    Bmat <- Q_tensor %*% Mid %*% t(Q_tensor)

    partial <- E_l %*% (diag(m^2) - Bmat %*% t(E_d) %*% solve(E_d %*% Bmat %*% t(E_d)) %*% E_d) %*% Bmat %*% t(E_l + E_u)
    L <- t(Wi) %*% t(partial)   # (d x #pairs)

    list(partial = partial, L = L)
  }

  # ---------------- psi2 (alpha equation) ----------------
  psi2i_only <- function(Yi, Xi, Wi, alpha, beta, phi, d) {
    m <- length(Yi)
    psi2 <- rep(0, d)
    if (m <= 1) return(psi2)

    eta_lin <- Xi %*% beta
    z <- h_fun(eta_lin)
    mu <- a1_fun(z)

    var_base <- a2_fun(z)
    sigma <- phi * var_base
    S <- diag(1 / sqrt(as.vector(sigma)), nrow = m)  # A^{-1/2}

    G <- recover_GZT(Wi %*% alpha)
    M <- exp_solve_mat(G)  # R^{-1}

    tmp <- build_partial_L(Wi, G, m)
    L <- tmp$L

    r <- matrix(as.vector(Yi - mu), ncol = 1)
    R_hat <- S %*% r %*% t(r) %*% S
    eta2 <- vecl(M %*% R_hat %*% M - M)

    as.vector(L %*% eta2)
  }

  # ---------------- psi1 (beta equation) ----------------
  # IMPORTANT: keep D consistent with calculate_beta() used for H1 if you replace A_bb by H1.
  psi1i_only <- function(Yi, Xi, Wi, alpha, beta, phi) {
    m <- length(Yi)
    eta_lin <- Xi %*% beta
    z <- h_fun(eta_lin)
    mu <- a1_fun(z)

    var_base <- a2_fun(z)
    sigma <- phi * var_base
    S <- diag(1 / sqrt(as.vector(sigma)), nrow = m)

    G <- recover_GZT(Wi %*% alpha)
    M <- exp_solve_mat(G)
    V_solve <- S %*% M %*% S

    if(family == "Gaussian") {
      D <- Xi
    }
    else if(family == "Poisson") {
      D <- diag(as.vector(mu), nrow = m) %*% Xi
    }
    else if(family == "Binomial") {
      D <- diag(as.vector(mu * (1 - mu)), nrow = m) %*% Xi
    }
    as.vector(t(D) %*% V_solve %*% (Yi - mu))
  }

  # ============================================================
  # 1) Build Psi and B_hat
  # ============================================================
  if (phi_known) {
    q <- p + d
    theta_hat <- c(beta, alpha)

    Psi <- matrix(0, nrow = q, ncol = n)
    for (i in 1:n) {
      Yi <- Y[[i]]; Xi <- X[[i]]; Wi <- W[[i]]
      psi1 <- psi1i_only(Yi, Xi, Wi, alpha, beta, phi)
      psi2 <- psi2i_only(Yi, Xi, Wi, alpha, beta, phi, d)
      Psi[, i] <- c(psi1, psi2)
    }
    B_hat <- (Psi %*% t(Psi)) / n

  } else {
    q <- p + d + 1
    theta_hat <- c(beta, alpha, phi)

    Psi <- matrix(0, nrow = q, ncol = n)
    for (i in 1:n) {
      Yi <- Y[[i]]; Xi <- X[[i]]; Wi <- W[[i]]
      psi1 <- psi1i_only(Yi, Xi, Wi, alpha, beta, phi)
      psi2 <- psi2i_only(Yi, Xi, Wi, alpha, beta, phi, d)
      psi3 <- psi3i_versionA(Yi = Yi, Xi = Xi, beta = beta, phi = phi, p = p, N = N)
      Psi[, i] <- c(psi1, psi2, psi3)
    }
    B_hat <- (Psi %*% t(Psi)) / n
  }

  # ============================================================
  # 2) Build A_hat
  #   - beta-row: only A_bb from H1/n, others 0
  #   - alpha row: closed-form A_ab, A_aphi; diff only A_aa
  #   - phi row: closed-form A_phib, A_phiphi (if unknown), A_phia=0
  # ============================================================
  A_hat <- matrix(0, nrow = q, ncol = q)
  idx_beta  <- 1:p
  idx_alpha <- (p + 1):(p + d)

  # ---- beta-row (as requested) ----
  H1 <- calculate_beta(Y, X, W, alpha, beta, phi, family = family)$H
  sgn <- if (use_minus) -1 else 1
  A_hat[idx_beta, idx_beta] <- sgn * (H1 / n)
  # A_hat[idx_beta, idx_alpha] = 0 already
  # if phi unknown: A_hat[idx_beta, idx_phi] remains 0

  # ============================================================
  # alpha-row closed-form blocks: A_{alpha,beta}, A_{alpha,phi}
  # ============================================================
  A_ab <- matrix(0, nrow = d, ncol = p)
  A_aphi <- rep(0, d)

  dpsi2_dbeta_cluster <- function(Yi, Xi, Wi, alpha, beta, phi) {
    m <- length(Yi)
    if (m <= 1) return(matrix(0, d, p))

    eta_lin <- Xi %*% beta
    z <- h_fun(eta_lin)
    mu <- as.vector(a1_fun(z))
    dmu_eta <- dmu_deta(eta_lin, mu)

    var_base <- as.vector(a2_fun(z))
    Vp_mu <- Vprime_from_mu(mu)  # derivative of a2 w.r.t mu

    sig <- phi * var_base
    s <- 1 / sqrt(sig)
    S <- diag(s, m, m)

    # D = dmu/dbeta = diag(dmu/deta) %*% X
    Dmat <- diag(as.vector(dmu_eta), m, m) %*% Xi

    r <- matrix(as.vector(Yi - mu), ncol = 1)
    RR <- r %*% t(r)

    G <- recover_GZT(Wi %*% alpha)
    M <- exp_solve_mat(G)

    tmp <- build_partial_L(Wi, G, m)
    L <- tmp$L

    sig_pow <- sig^(-3/2)
    npair <- m * (m - 1) / 2
    dEta <- matrix(0, nrow = npair, ncol = p)

    for (j in 1:p) {
      dr <- -matrix(Dmat[, j], ncol = 1)

      # d sigma / d beta_j = phi * a2'(mu) * dmu/deta * X_{.,j}
      dsig <- phi * (Vp_mu * dmu_eta * Xi[, j])
      ds <- -0.5 * sig_pow * dsig
      dS <- diag(as.vector(ds), m, m)

      dRhat <- dS %*% RR %*% S +
        S %*% RR %*% dS +
        S %*% (dr %*% t(r) + r %*% t(dr)) %*% S

      dMmat <- M %*% dRhat %*% M
      dEta[, j] <- vecl(dMmat)
    }

    L %*% dEta
  }

  dpsi2_dphi_cluster <- function(Yi, Xi, Wi, alpha, beta, phi) {
    m <- length(Yi)
    if (m <= 1) return(rep(0, d))

    eta_lin <- Xi %*% beta
    z <- h_fun(eta_lin)
    mu <- a1_fun(z)

    var_base <- a2_fun(z)
    sig <- phi * var_base
    S <- diag(1 / sqrt(sig), nrow = m)

    G <- recover_GZT(Wi %*% alpha)
    M <- exp_solve_mat(G)

    tmp <- build_partial_L(Wi, G, m)
    L <- tmp$L

    r <- matrix(as.vector(Yi - mu), ncol = 1)
    R_hat <- S %*% r %*% t(r) %*% S

    dMmat <- M %*% R_hat %*% M
    dEta_phi <- (-1/phi) * vecl(dMmat)

    as.vector(L %*% dEta_phi)
  }

  for (i in 1:n) {
    Yi <- Y[[i]]; Xi <- X[[i]]; Wi <- W[[i]]
    A_ab <- A_ab + dpsi2_dbeta_cluster(Yi, Xi, Wi, alpha, beta, phi)
    if (!phi_known) A_aphi <- A_aphi + dpsi2_dphi_cluster(Yi, Xi, Wi, alpha, beta, phi)
  }
  A_hat[idx_alpha, idx_beta] <- (A_ab / n)

  if (!phi_known) {
    idx_phi <- p + d + 1
    A_hat[idx_alpha, idx_phi] <- (A_aphi / n)
  }

  # ============================================================
  # A_{alpha,alpha}: ONLY finite-difference (d x d)
  # ============================================================
  U_alpha_sum <- function(a) {
    out <- rep(0, d)
    for (i in 1:n) {
      out <- out + psi2i_only(Y[[i]], X[[i]], W[[i]], a, beta, phi, d)
    }
    out
  }

  for (k in 1:d) {
    step <- eps
    if (scale_eps) step <- eps * max(1, abs(alpha[k]))
    a_p <- alpha; a_m <- alpha
    a_p[k] <- a_p[k] + step
    a_m[k] <- a_m[k] - step
    Up <- U_alpha_sum(a_p)
    Um <- U_alpha_sum(a_m)
    A_hat[idx_alpha, idx_alpha[k]] <- (Up - Um) / (2 * step) / n
  }

  # ============================================================
  # phi-row (if phi unknown): closed-form
  #   A_{phi,alpha}=0, A_{phi,beta}, A_{phi,phi}
  # ============================================================
  if (!phi_known) {
    idx_phi <- p + d + 1

    c_sum <- 0
    A_phib <- rep(0, p)

    for (i in 1:n) {
      Yi <- as.vector(Y[[i]])
      Xi <- X[[i]]
      m <- length(Yi)
      ci <- m - (p * m / N)
      c_sum <- c_sum + ci

      eta_lin <- Xi %*% beta
      z <- h_fun(eta_lin)
      mu <- as.vector(a1_fun(z))
      dmu_eta <- dmu_deta(eta_lin, mu)

      Vmu <- as.vector(a2_fun(z))          # var without phi
      Vp_mu <- Vprime_from_mu(mu)          # derivative wrt mu

      r <- Yi - mu

      # derivative of sum (r^2 / V) wrt beta:  X^T * [(-2 r / V - r^2 V'/V^2) * dmu/deta]
      w <- (-2 * r / Vmu - (r^2) * Vp_mu / (Vmu^2)) * dmu_eta
      A_phib <- A_phib + as.vector(t(Xi) %*% w)
    }

    A_hat[idx_phi, idx_beta] <- (A_phib / n)
    A_hat[idx_phi, idx_phi]  <- (-c_sum / n)
    # A_hat[idx_phi, idx_alpha] remains 0
  }

  # ============================================================
  # 3) Sandwich variance
  # ============================================================
  Ainv <- solve(A_hat)
  Sigma_hat <- (Ainv %*% B_hat %*% t(Ainv)) / n
  if (force_sym) Sigma_hat <- 0.5 * (Sigma_hat + t(Sigma_hat))

  # ============================================================
  # 4) outputs
  # ============================================================
  if (phi_known) {
    list(
      phi_known = TRUE,
      theta_hat = c(beta, alpha),
      A_hat = A_hat,
      B_hat = B_hat,
      Sigma_hat = Sigma_hat,
      cov_beta  = Sigma_hat[idx_beta,  idx_beta,  drop = FALSE],
      cov_alpha = Sigma_hat[idx_alpha, idx_alpha, drop = FALSE],
      cov_beta_alpha = Sigma_hat[idx_beta, idx_alpha, drop = FALSE],
      Psi = Psi
    )
  } else {
    idx_phi <- p + d + 1
    list(
      phi_known = FALSE,
      theta_hat = c(beta, alpha, phi),
      A_hat = A_hat,
      B_hat = B_hat,
      Sigma_hat = Sigma_hat,
      cov_beta  = Sigma_hat[idx_beta,  idx_beta,  drop = FALSE],
      cov_alpha = Sigma_hat[idx_alpha, idx_alpha, drop = FALSE],
      var_phi   = Sigma_hat[idx_phi, idx_phi],
      cov_beta_phi  = Sigma_hat[idx_beta,  idx_phi, drop = FALSE],
      cov_alpha_phi = Sigma_hat[idx_alpha, idx_phi, drop = FALSE],
      Psi = Psi
    )
  }
}
