source("samplers.R")

######################### 
# Multivariate Gaussian # You can un-comment this part to run simulation on "Multivariate Gaussian" model
#########################______________________________________________
# Mu        = c(3, 6)
# Sigma     = matrix(c(2, .5, .5, 1), nrow=2, ncol=2)
# Sigma.inv = solve(Sigma)
# Params    = list('mu'=Mu, 'sigma'=Sigma, 'sigma.inv'=Sigma.inv)
# 
# eval  = eigen(Sigma)$values
# evec  = eigen(Sigma)$vectors
# L     = 1/min(eval)
# m     = 1/max(eval) 
# 
# # Calculates the log-likelihood
# loglike <- function(x, params=Params){
#   mu        = params$mu
#   sigma.inv = params$sigma.inv
#   
#   return(as.numeric(- t((x - mu)) %*% sigma.inv %*% (x - mu)/2) )
# }
# 
# # Calculate the gradient of the log-likelihood
# loglike.grad <- function(x, params=Params){
#   mu        = params$mu
#   sigma.inv = params$sigma.inv
#   
#   return(as.numeric(sigma.inv %*% (x - mu)))
# }

#############################
# Bayes Posterior Inference #
#############################__________________________________________
# General params settings
dim.theta = 2
n.exps    = 100 #for "mainSingle.R" this parameter is useless
n.samples = 50
n.iters   = 1e3*2
pre.cond  = FALSE
lazy.step = FALSE

Alpha     = 0.9  #user-specified (for model) 
#__________________________________________________
# Data Generation
source("geneSamples.R")
samples = geneSamples(n=n.samples, p=dim.theta)
X = matrix(samples$X, ncol=dim.theta)
Y = matrix(samples$Y)
#__________________________________________________
# Data Info
n     = dim(X)[1]  #number of samples
D     = dim(X)[2]  #dimension of theta
Theta = matrix(rep(1, D)) #real theta
covx  = t(X) %*% X / n
eval  = eigen(covx)$values
evec  = eigen(covx)$vectors
#__________________________________________________
# Model (Logistic Regression) params 
if (pre.cond==TRUE){
  covxp = evec%*%diag((eval)^( 1/2))%*%t(evec)  #covx^( 1/2)
  covxn = evec%*%diag((eval)^(-1/2))%*%t(evec)  #covx^(-1/2)
  L     = Alpha+n/4
  m     = Alpha
} else{
  covxp = diag(D)
  covxn = diag(D)
  L     = max(eval) * (Alpha+n/4)
  m     = min(eval) *  Alpha
}
#__________________________________________________
# Model (Logistic Regression) functions

# The netative of log likelihood
f <- function(theta, alpha = Alpha){
  # pre-condition
  theta = covxn %*% theta
  
  return( - t(Y) %*% X %*% theta 
          + sum(log(1 + exp(X %*% theta))) 
          + alpha * t(theta) %*% covx %*% theta)
}

# Log Likelihood
loglike <- function(theta){return(-f(theta))}

# The gradient of f 
fgrad <- function(theta, alpha = Alpha){
  # pre-condition
  theta = covxn %*% theta
  
  theta = matrix(theta)
  temp  = 0
  for(i in 1:n){
    temp = temp + X[i,] / as.vector(1 + exp(- X[i,] %*% theta))
  }
  return( - t(X) %*% Y 
          + temp 
          + alpha * covx %*% theta * 2)
}

# The gradient of Log Likelihood
loglike.grad <- function(theta){return(-fgrad(theta))}


####################### The above is the same #########################
####################### for two main.R files. #########################

#################
#  Experiments  # without init_distr for staring points (all sampling starts at c(0, 0))
#################
#______________________________________________________________________
# Sampling

# ULA large
out.ula0 <- ula(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, delta=0.1,
                   loglike=loglike, loglike.grad=loglike.grad)
out.ula0$accept

# ULA
out.ula <- ula(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, delta=1,
               loglike=loglike, loglike.grad=loglike.grad)
out.ula$accept
# MALA
out.mala <- mala(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step,
                 loglike=loglike, loglike.grad=loglike.grad)
out.mala$accept

out.mrw <- mrw(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, #sigma=2.5,
                 loglike=loglike)
out.mrw$accept

# out.mala$accept    for: accept rate
# out.mala$chain[,1] for: theta 1
# out.mala$chain[,2] for: theta 2

#______________________________________________________________________
# Visualization
source("plotsSingle.R")
chain.list = list(out.ula0, out.ula, out.mala, out.mrw)

plt.density(chain.list)
plt.acfccf(chain.list)
plt.trace(chain.list)

#______________________________________________________________________
# Benchmarking
# library(rbenchmark)
# benchmark(ula(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, delta=0.1,
#                loglike=loglike, loglike.grad=loglike.grad, init=c(0, 0)), 
#            replications=10)
# 
# benchmark(mala(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step,
#                 loglike=loglike, loglike.grad=loglike.grad, init=c(0, 0)), 
#            replications=10)
# benchmark(mrw(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, #sigma=2.5,
#                  loglike=loglike, init=c(0, 0)), 
#            replications=10)



