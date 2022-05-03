prior_predict <- function(n, p, alpha, priors, names){
  # sample from prior predictive, mainly for purpose of comparing to posterior predictive
  
  # packages
  require(Rcpp)
  require(RcppArmadillo)
  
  # compile C++ functions
  print("Compiling C++ Functions")
  sourceCpp("rcppfuncts/rdirichlet_arma.cpp")
  
  # preliminaries
  K <- length(alpha)
  
  # sample eta
  eta <- rdirichlet_arma(n = 1, a = alpha)
  
  # sample cluster labels
  s <- sample(1:K, size = n, replace = T, prob = eta)
  
  # sample prior parameters
  theta <- array(0, dim = c(K, 2, p))
  for (j in 1:p){
    for (h in 1:K) {
      zeta <- priors[h, 1, j]
      lambda <- priors[h, 2, j]
      psi <- priors[h, 3, j]
      nu <- priors[h, 4, j]
      theta[h, 2, j] <- 1/rgamma(n = 1, shape = psi, rate = nu)
      theta[h, 1, j] <- rnorm(n = 1, mean = zeta, sd = sqrt(theta[h,2,j]/lambda))
    }
  }
  
  # sample data
  y <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    for (j in 1:p){
      y[i,j] <- rnorm(n = 1, mean = theta[s[i], 1, j], sd = sqrt(theta[s[i], 2, j]))
    }
  }
  y <- as.data.frame(y)
  colnames(y) <- names
  return(list(y = y, s = s, theta = theta, eta = eta))
}