partA <- function() {
  par(mfrow = c(2,2))
  
  x <- list()
  
  set.seed(2)
  lambda <- 2
  x[[1]] <- rpois(100,lambda)
  barplot(table(x[[1]]), xlab = "Values", ylab = "Frequency", main = "Poison with lambda=2, 100 samples")
  
  set.seed(2)
  lambda <- 2
  x[[2]] <- rpois(10000,lambda)
  barplot(table(x[[2]]), xlab = "Values", ylab = "Frequency", main = "Poison with lambda=2, 10000 samples")
  
  y <- list()
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[1]] <- rgamma(100, shape = alpha, scale = beta)
  hist(
    y[[1]], breaks = 10, xlab = "Values", ylab = "Frequency", main = "Gamma with a=2 and b=3, 100 samples"
  )
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[2]] <- rgamma(10000, shape = alpha, scale = beta)
  hist(
    y[[2]], breaks = 100, xlab = "Values", ylab = "Frequency", main = "Gamma with a=2 and b=3, 10000 samples"
  )
}

partB <- function() {
  x <- list()
  
  set.seed(2)
  lambda <- 2
  x[[1]] <- rpois(100,lambda)
  
  set.seed(2)
  lambda <- 2
  x[[2]] <- rpois(10000,lambda)
  
  y <- list()
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[1]] <- rgamma(100, shape = alpha, scale = beta)
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[2]] <- rgamma(10000, shape = alpha, scale = beta)
  
  means <- lapply(X = c(x[1],x[2],y[1],y[2]), FUN = mean)
  sds <- lapply(X = c(x[1],x[2],y[1],y[2]), FUN = sd)
  medians <- lapply(X = c(x[1],x[2],y[1],y[2]), FUN = median)
  iqrs <- lapply(X = c(x[1],x[2],y[1],y[2]), FUN = IQR)
  
  real_means <- c(lambda, lambda, alpha * beta, alpha * beta)
  real_sds <-
    c(sqrt(lambda), sqrt(lambda), sqrt(alpha * beta ^ 2), sqrt(alpha * beta ^
                                                                 2))
  real_medians <-
    c(
      qpois(0.5,lambda = lambda), qpois(0.5,lambda = lambda), qgamma(0.5, shape =
                                                                       alpha,scale = beta), qgamma(0.5, shape = alpha,scale = beta)
    )
  real_iqrs <-
    c(
      qpois(0.75,lambda = lambda) - qpois(0.25,lambda = lambda), qpois(0.75,lambda = lambda) -
        qpois(0.25,lambda = lambda), qgamma(0.75, shape = alpha, scale = beta) -
        qgamma(0.25, shape = alpha,scale = beta),qgamma(0.75, shape = alpha, scale = beta) -
        qgamma(0.25, shape = alpha,scale = beta)
    )
  
  data <-
    cbind(means, real_means, medians,real_medians, sds, real_sds, iqrs, real_iqrs)
  data <- as.data.frame(data)
  rownames(data) <-
    c("Poison-100", "Poisson-10000", "Gamma-100", "Gamma-10000")
  colnames(data) <-
    c(
      "Data Mean", "Distribution Mean", "Data Median", "Distribution Median", "Data SD", "Distribution SD", "Data IQR", "Distribution IQR"
    )
  
  print(data)
}

partC <- function() {
  x <- list()
  
  set.seed(2)
  lambda <- 2
  x[[1]] <- rpois(100,lambda)
  
  set.seed(2)
  lambda <- 2
  x[[2]] <- rpois(10000,lambda)
  
  y <- list()
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[1]] <- rgamma(100, shape = alpha, scale = beta)
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[2]] <- rgamma(10000, shape = alpha, scale = beta)
  
  means <- lapply(X = c(x[1],x[2],y[1],y[2]), FUN = mean)
  medians <- lapply(X = c(x[1],x[2],y[1],y[2]), FUN = median)
  real_means <- c(lambda, lambda, alpha * beta, alpha * beta)
  real_medians <-
    c(
      qpois(0.5,lambda = lambda), qpois(0.5,lambda = lambda), qgamma(0.5, shape =
                                                                       alpha,scale = beta), qgamma(0.5, shape = alpha,scale = beta)
    )
  proportion <- function(x, mean, median) {
    length(x[x >= min(mean,median) & x <= max(mean,median)]) / length(x)
  }
  data_proportions <- mapply(proportion, c(x,y), means, medians)
  real_proportions <-
    c(sum(dpois(seq(
      min(real_means[1],real_medians[1]), max(real_means[1],real_medians[1])
    ),lambda = lambda)),sum(dpois(seq(
      min(real_means[2],real_medians[2]), max(real_means[2],real_medians[2])
    ),lambda = lambda)), diff(pgamma(
      c(
        min(real_means[4],real_medians[4]),max(real_means[3],real_medians[3])
      ),shape = alpha, scale = beta
    )), diff(pgamma(
      c(
        min(real_means[4],real_medians[4]),max(real_means[4],real_medians[4])
      ),shape = alpha, scale = beta
    )))
  data <- cbind(data_proportions, real_proportions)
  data <- as.data.frame(data)
  rownames(data) <-
    c("Poison-100", "Poisson-10000", "Gamma-100", "Gamma-10000")
  colnames(data) <- c("Data", "Distribution")
  data
}

partD <- function() {
  x <- list()
  
  set.seed(2)
  lambda <- 2
  x[[1]] <- rpois(100,lambda)
  
  set.seed(2)
  lambda <- 2
  x[[2]] <- rpois(10000,lambda)
  
  y <- list()
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[1]] <- rgamma(100, shape = alpha, scale = beta)
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[2]] <- rgamma(10000, shape = alpha, scale = beta)
  
  data_quantiles <- lapply(c(x,y), quantile, probs = 0.01)
  real_quantiles <-
    c(
      qpois(0.01, lambda = lambda), qpois(0.01, lambda = lambda), qgamma(0.01, shape = alpha, scale = beta), qgamma(0.01, shape = alpha, scale = beta)
    )
  
  data <- cbind(data_quantiles,real_quantiles)
  data <- as.data.frame(data)
  rownames(data) <-
    c("Poison-100", "Poisson-10000", "Gamma-100", "Gamma-10000")
  colnames(data) <- c("Data 1%", "Distribution 1%")
  data
}

partE <- function() {
  x <- list()
  
  set.seed(2)
  lambda <- 2
  x[[1]] <- rpois(100,lambda)
  
  set.seed(2)
  lambda <- 2
  x[[2]] <- rpois(10000,lambda)
  
  y <- list()
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[1]] <- rgamma(100, shape = alpha, scale = beta)
  
  set.seed(2)
  alpha <- 2
  beta <- 3
  y[[2]] <- rgamma(10000, shape = alpha, scale = beta)
  
  data_quantiles <- lapply(c(x,y), quantile, probs = 0.99)
  real_quantiles <-
    c(
      qpois(0.99, lambda = lambda), qpois(0.99, lambda = lambda), qgamma(0.99, shape = alpha, scale = beta), qgamma(0.99, shape = alpha, scale = beta)
    )
  
  data <- cbind(data_quantiles,real_quantiles)
  data <- as.data.frame(data)
  rownames(data) <-
    c("Poison-100", "Poisson-10000", "Gamma-100", "Gamma-10000")
  colnames(data) <- c("Data 99%", "Distribution 99%")
  data
}

partF <- function() {
  print(
    "In parts B to E above we observer that there are differences between the values computed from the samples and the values calculated from the theoretic distributions. This is to be expected as the randomization process creates these differences. Also we can observe that as the sample increases the estimates align more and more with the theoretical parameters of the distributions."
  )
}

partG< function(){
  print("The best prediction we can make is the mean of the sample. See part B for the values of each sample.")
}