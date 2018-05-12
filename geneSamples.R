# Generate Simulation data
library(extraDistr)   #rsign()
library(sigmoid)      #sigmoid()
geneSamples <- function(n, p=2){
  X <- matrix(nrow=n, ncol=p)
  Y <- matrix(nrow=n, ncol=1)
  Theta = matrix(rep(1, p)) 
  for (i in 1:p){
    X[,i] = rsign(n)  #draw samples from i.i.d. Rademacher distribution 
    #library(extraDistr)
  }
  for (i in 1:n){
    Y[i] = rbinom(1, 1, sigmoid(X[i, ] %*% Theta))}
  return (list("X"=X, "Y"=Y))
}

# samples = geneSamples(n=n.samples, p=2)
# X = matrix(samples$X, ncol=2)
# Y = matrix(samples$Y)

# # [python]
# # data generation  
# X = np.random.binomial(1, 0.5, (n, p)) * 2 -1
# X = X/np.sqrt((X * X).mean(axis=1))[:, None]
# theta_true = np.ones(d)
# h_true = X.dot(theta_true)
# r_true = 1./(1.+np.exp(-h_true))
# Y = np.random.binomial(1, r_true)