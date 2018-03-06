library(Rcpp)
library(microbenchmark)


# basic example -----------------------------------------------------------

# inefficient R implementation for squaring a vector
square_R <- function(vec) {
  n <- length(vec)
  squared_vec <- numeric(n)
  for (i in 1:n) {
    squared_vec[i] <- vec[i] * vec[i]
  }
  squared_vec
}

# create c++ function that squares a vector.
# cppFunction() handles compilation and makes
# it available 
cppFunction("NumericVector square_rcpp(NumericVector vec) {
              int n = vec.size();
              NumericVector x(n);
              for (int i=0; i < n; i++ ) {
                x[i] = vec[i]*vec[i];
              }
              return x;
            }")

# make sure it gives the correct values
vec <- runif(1000)
all.equal(square_rcpp(vec), square_R(vec), vec^2)

# even though ^2 is highly optimized, square_rcpp is pretty close
microbenchmark(vec^2, square_rcpp(vec), square_R(vec))


# two sample t-test simulation --------------------------------------------

# First, we define the simulation in R
t_test <- function(y1, y2, alpha = 0.05, method = c("pooled", "welch", "lln")) {
  m1 <- mean(y1)
  m2 <- mean(y2)
  v1 <- var(y1)
  v2 <- var(y2)
  n1 <- length(y1)
  n2 <- length(y2)
  if (method == "pooled") {
    df <- n1 + n2 - 2
    v_pool <- ((n1-1) * v1 + (n2-1) * v2) / df
    test_statistic <- (m1 - m2) / sqrt(v_pool * (1 / n1 + 1 / n2))
    critical_value <- qt(1 - alpha / 2, df)
    result <- ifelse(abs(test_statistic) >= critical_value, 1, 0) 
  } else if (method == "welch") {
    k <- (v1 / n1) / (v1 / n1 + v2 / n2)
    df <- (n1 - 1) * (n2 - 1) / ((1 - k) * (1 - k) * (n1 - 1) + k * k * (n2-1))
    test_statistic <- (m1 - m2) / sqrt(v1 / n1 + v2 / n2)
    critical_value <- qt(1 - alpha / 2, df)
    result <- ifelse(abs(test_statistic) >= critical_value, 1, 0)
  } else if (method == "lln") {
    test_statistic <- (m1 -  m2) / sqrt(v1 / n1 + v2 / n2)
    critical_value <- qnorm(1 - alpha / 2)
    result <- ifelse(abs(test_statistic) >= critical_value, 1, 0)
  }
  result
}

simulation_R <- function(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha) {
  engine <- function(n1, n2, mu1, mu2, sigma1, sigma2, alpha) {
    y1 <- rnorm(n1, mu1, sigma1)
    y2 <- rnorm(n2, mu2, sigma2)
    result <- c(t_test(y1, y2, alpha, "pooled"),
                t_test(y1, y2, alpha, "lln"),
                t_test(y1, y2, alpha, "welch"))
    result
  }
  
  output <- matrix(nrow = length(sigma2) * length(mu2), ncol = 5)
  colnames(output) <- c("mu2", "sigma2", "pool", "lln", "welch")
  counter <- 1
  for (i in 1:length(sigma2)) {
    for (j in 1:length(mu2)) {
      result <- rowSums(sapply(1:nreps, function(x) { 
          engine(n1, n2, mu1, mu2[j], sigma1, sigma2[i], alpha) 
        }))
      output[counter, ] <- c(mu2[j], sigma2[i], result / nreps)
      counter <- counter + 1
    }
  }
  output
}

# Next we define the simulation in C 
simulation_C <- function(n1,n2,mu1,sigma1,mu2,sigma2,nreps,alpha=0.05) {
  if ( !is.loaded("welch1938") ) dyn.load("welch1938.so")
  if ( length(n1) != 1 ) stop("n1 must be a vector of length 1.")
  if ( length(n2) != 1 ) stop("n2 must be a vector of length 1.")
  if ( length(mu1) == 0 ) stop("mu1 must have non zero length.")
  if ( length(sigma1) == 0 ) stop("sigma1 must have non zero length.")
  if ( length(mu2) == 0 ) stop("mu2 must have non zero length.")
  if ( length(sigma2) == 0 ) stop("sigma2 must have non zero length.")
  if ( length(nreps) != 1 ) stop("nreps must be a vector of length 1.")
  if ( nreps <= 0 ) stop("nreps must be a positive integer.")
  if ( length(alpha) != 1 ) stop("alpha must be a vector of length 1.")
  if ( alpha <= 0.0 || alpha >= 1.0 ) stop("alpha must be between 0 and 1.")
  result <- matrix(0,ncol=length(mu2)*length(sigma2),nrow=5)
  storage.mode(result) <- "double"
  all <- .C("simulation", n1=as.integer(n1), n2=as.integer(n2),
            mu1=as.double(mu1), sigma1=as.double(sigma1),
            mu2=as.double(mu2), mu2Length=as.integer(length(mu2)),
            sigma2=as.double(sigma2), sigma2Length=as.integer(length(sigma2)),
            nreps=as.integer(nreps), alpha=as.double(alpha), result=result)
  result <- t(all$result)
  colnames(result) <- c("mu2","sigma2","pool","lln","welch")
  result
}

# Finally, we load our c++ implementation
sourceCpp("./welch1938.cpp")

n1 <- 6
n2 <- 10
mu1 <- 0
sigma1 <- 1
mu2 <- c(0.0,0.5,1.0,2.0)
sigma2 <- c(0.25,1.0,2.0,5.0)
nreps <- 1000
alpha <- 0.1

(sim_results_R <- simulation_R(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha))
(sim_results_C <- simulation_C(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha))
(sim_results_rcpp <- simulation_rcpp(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha))

# compare computation time
microbenchmark(simulation_R(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha), 
               simulation_C(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha),
               simulation_rcpp(n1, n2, mu1, sigma1, mu2, sigma2, nreps, alpha),
               times=10)


# Gibbs sampler -----------------------------------------------------------

n <- 100
n_beta <- 6
X <- cbind(1, matrix(rnorm(n * (n_beta - 1), sd = 4), nrow = n))
beta <- rnorm(n_beta, runif(n_beta, -5, 5), sd = 1)
sigma2 <- 2.5
y <- X %*% beta + rnorm(n, sd = sqrt(sigma2))

gibbs_sampler <- function(X, y, ndraws = 1000, beta_s2 = 10) {
  beta_draws <- matrix(0, nrow = ndraws, ncol = ncol(X))
  sigma2_draws <- rep(1, ndraws)
  for (i in 2:ndraws) {
    Lambda <- (1/sigma2_draws[i-1] * t(X) %*% X + 1/beta_s2)
    Sigma <- solve(Lambda)
    mu <- 1 / sigma2_draws[i-1] * Sigma %*% t(X) %*% y
    beta_draws[i, ] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
    shape <- 2.001 + length(y) / 2
    rate <- 1.001 + 0.5 * t(y - X %*% beta_draws[i, ]) %*% (y - X %*% beta_draws[i, ])
    sigma2_draws[i] <- 1 / rgamma(1, shape = shape, rate = rate)
  }
  list(beta_draws = beta_draws, sigma2_draws = sigma2_draws)
}

mcmc_out <- gibbs_sampler(X = X, y = y)

sourceCpp("./gibbs_sampler.cpp")
mcmc_out_cpp <- gibbs_sampler_rcpp(X = X, y = as.numeric(y), ndraws = 1000, beta_s2 = 10)

            