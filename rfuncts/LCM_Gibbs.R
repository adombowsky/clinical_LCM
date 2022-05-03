sample.theta <- function(y, type, pars, s, K){
  # samples full conditional of theta
  if (type == "c"){ # continuous
    theta <- matrix(0, nrow = K, ncol = 2)
    # key:
    for (h in 1:K){
      zeta <- pars[h, 1]
      lambda <- pars[h, 2]
      psi <- pars[h, 3]
      nu <- pars[h, 4]
      # update full conditional
      n_hj <- sum(s == h)
      y_bar_hj <- ifelse(n_hj > 0, mean(y[s==h]), 0)
      zeta_hj <- (zeta * lambda + n_hj * y_bar_hj)/(lambda + n_hj)
      lambda_hj <- lambda + n_hj
      psi_hj <- psi + 0.5 * n_hj
      nu_hj <- nu + 0.5 * sum( (y - y_bar_hj)^2 * (s == h)) +
        ((lambda * n_hj)/(lambda + n_hj)) * 0.5 * (y_bar_hj - zeta)^2
      # sample theta
      theta[h, 2] <- 1/rgamma(n = 1, shape = psi_hj, rate = nu_hj)
      theta[h, 1] <- rnorm(n = 1, mean = zeta_hj, sd = sqrt(theta[h, 2]/lambda_hj))
    }
  }
  else{ # binary
    theta <- rep(0, K)
    # key:
    a <- pars[1]
    b <- pars[2]
    for (h in 1:K){
      # update full conditional
      a_hj <- a + sum(y[s == h])
      b_hj <- b + sum(s == h) - sum(y[s == h])
      # sample theta
      theta[h] <- rbeta(n = 1, shape1 = a_hj, shape2 = b_hj)
    }
  }
  return(theta)
}

sample.y_ij <- function(theta_j, type, s_i){
  if (type == "c") {
    y_ij <- rnorm(n = 1, 
                  mean = theta_j[s_i, 1], 
                  sd = sqrt(theta_j[s_i, 2]))
    
  }
  else if (type == "b") {
    y_ij <- rbinom(n = 1,
                   size = 1,
                   prob = theta_j[s_i])
  }
  return(y_ij)
}



sample.marginal.s_i <- function(y_i, eta, priors, K, comb_theta, p){
  sigma2 <- comb_theta[, seq(2, 2*p, by = 2)]
  probs <- rep(0,K)
  for (h in 1:K){
    zeta <- priors[h,1,]
    lambda <- priors[h, 2,]
    probs[h] <- sum(dnorm(x = y_i, 
                          mean = zeta, 
                          sd = sqrt(sigma2[h,]*(1 + 1/lambda)), 
                          log = TRUE)) + log(eta[h])
  }
  probs <- unlist(probs)
  probs <- exp(probs - max(probs))
  s_i <- sample(x = 1:K, size = 1, replace = T, prob = probs)
  return(s_i)
}

sample.marginal.s_i.tdistib <- function(y_i, eta, priors, K){
  probs <- rep(0, K)
  for (h in 1:K){
    zeta <- priors[h, 1, ]
    lambda <- priors[h, 2, ]
    psi <- priors[h, 3, ]
    nu <- priors[h, 4, ]
    tau <- sqrt( (nu/psi) * (1 + (1/lambda))) 
    y_icent <- (y_i - zeta)/tau
    probs[h] <- sum(dt(x = y_icent, df = 2 * psi, log = T) - log(tau)) + log(eta[h])
  }
  probs <- exp(probs - max(probs))
  s_i <- sample(x = 1:K, size = 1, replace = T, prob = probs)
  return(s_i)
  
}

apply.marginal.s_i <- function(x, eta, priors, K, comb_theta, p){
  s <- sample.marginal.s_i(y_i = x, eta = eta, 
                           priors = priors, K = K, comb_theta = comb_theta, p = p)
  return(s)
}

apply.marginal.s_i.tdistrib <- function(x, eta, priors, K){
  s <- sample.marginal.s_i.tdistib(y_i = x, eta = eta,  priors = priors, K = K)
  return(s)
}

sample.eta <- function(alpha, s, K){
  n_i <- table(factor(s, levels = 1:K))
  alpha_star <- n_i + alpha
  eta <- rdirichlet_arma(n = 1, a = alpha_star)
  return(eta)
}


theta.list <- function(type, K, R){
  # constructs theta array
  p <- length(type)
  theta <- list()
  for (j in 1:p){
    if (type[j] == "c"){
      theta[[j]] <- array(1, dim = c(K, 2, R))
    }
    else if (type[j] == "b") {
      theta[[j]] <- matrix(0.5, nrow = R, ncol = K)
    }
  }
  return(theta)
}

inf_LCM_gibbs <- function(R, y, type, K,
                          priors,
                          alpha,
                          stops = 1){
  # Gibbs sampler for the GoM model w/ mixed data
  # packages
  require(Rcpp)
  require(RcppArmadillo)
  require(abind)
  
  # C++ functions
  print("Compiling C++ Functions")
  sourceCpp("rcppfuncts/rdirichlet_arma.cpp")
  
  # preliminaries
  n <- nrow(y)
  p <- ncol(y)
  
  # storage
  theta <- theta.list(type = type, K = K, R = R)
  s <- matrix(sample(1:K, size = n * R, replace = T, prob = rep(1/K,K)), nrow = R, ncol = n)
  eta <- matrix(1/K, nrow = R, ncol = K)
  y_array <- array(0, dim = c(R, n, p))
  y_array[1, , ] <- as.matrix(y)
  
  # configuring stops
  stops <- (1:(R/stops)) * stops
  
  # configure missing data
  na_matrix <- is.na(y)
  for (i in 1:n) {
    for (j in 1:p) {
      if (na_matrix[i, j]) {
        if (type[j] == "c") {
          y[i, j] = 0
        }
        else {
          y[i, j] = rbinom(n = 1, size = 1, prob = 1/2)
        }
      }
    }
  }
  na_logical = apply(na_matrix, 1, function(x) sum(x) > 0)
  na_index = (1:n)[na_logical]
  
  # sampling
  print("Sampling")
  for (r in 2:R){
    for (j in 1:p){
      # sample theta
      #if (type[j] == "c"){
        samp <- sample.theta(y = y[,j], type = "c", pars = priors[,,j], s = s[r-1,], K = K)
        theta[[j]][, 1, r] <- samp[, 1]
        theta[[j]][, 2, r] <- samp[, 2]
        # now missing data
        for (i in na_index){
          if (na_matrix[i, j]) {
            y[i,j] <- as.numeric(sample.y_ij(theta_j = theta[[j]][, ,r], 
                                             type = "c", s_i = s[r-1,i]))
          }
          }
        #}
      #else if (type[j] == "b"){
        #theta[[j]][r, ] <- sample.theta(y = y[,j], type = "b", pars = b_pars, s = s[r-1,], K = K)
        # now update labels
        #for (i in 1:n){
          #if (na_matrix[i, j]) {
            #y_ij <- sample.y_ij(theta_j = theta[[j]][r, ], 
                                #type = "b", s_ij = s[i,j])
          #}
          #}
          #s[i, j] <- sample.s_ij(y_ij = y_ij, type = "b", eta_i = eta[r-1, , i], 
                                 #theta_j = theta[[j]][r, ], K = K)
        #}
      }
    # now, update labels and weights
    s[r, ] <- apply(y, 1, apply.marginal.s_i.tdistrib, eta = eta[r-1,], priors = priors, 
                    K = K)
    eta[r, ] <- sample.eta(alpha = alpha, s = s[r,], K = K)
    #for (i in 1:n){
      #s[r, i] <- sample.marginal.s_i(y_i = y[i,], eta_i = eta[r-1, , i], priors = priors, K = K, theta = theta, r = r, p = p)
      #eta[r, , i] <- sample.eta_i(alpha = alpha, s_i = s[r,i], K = K)
    #}
    # save missing data
    y_array[r, , ] = as.matrix(y)
    
    # print stops
    if (r %in% stops){
      print(r)
    }
  }
  # returning
  return(list(s = s,theta = theta, eta = eta, y_array = y_array))
}