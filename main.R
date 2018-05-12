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
n.exps    = 500
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
#  Experiments  # with init_distr for staring points 
#################
#______________________________________________________________________
# Sampling
# Distribution of the (sampling) starting point
init_distr_L = 1./sqrt(L) * matrix(rnorm(n.exps*dim.theta), ncol = dim.theta)
init_distr_m = 1./sqrt(m) * matrix(rnorm(n.exps*dim.theta), ncol = dim.theta)

out.ula0 = list()
acc.ula0 = list()
out.ula  = list()
acc.ula  = list()
out.mala = list()
acc.mala = list()
out.mrw  = list()
acc.mrw  = list()

for (i in 1:n.exps) {
  
  cat("iter: ", i, "\n")
  # ULA large
  chain.ula0 <- ula(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, delta=0.1,
                     loglike=loglike, loglike.grad=loglike.grad, init=init_distr_m[i, ])
  out.ula0[[i]] = chain.ula0$chain
  acc.ula0[[i]] = chain.ula0$accept

  # ULA
  chain.ula <- ula(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, delta=1,
                 loglike=loglike, loglike.grad=loglike.grad, init=init_distr_L[i, ])
  out.ula[[i]] = chain.ula$chain
  acc.ula[[i]] = chain.ula$accept
  
  # MALA
  chain.mala <- mala(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step,
                   loglike=loglike, loglike.grad=loglike.grad, init=init_distr_L[i, ])
  out.mala[[i]] = chain.mala$chain
  acc.mala[[i]] = chain.mala$accept
  
  # MRW
  chain.mrw <- mrw(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, #sigma=2.5,
                   loglike=loglike, init=init_distr_L[i, ])
  out.mrw[[i]] = chain.mrw$chain
  acc.mrw[[i]] = chain.mrw$accept
}
#______________________________________________________________________
# Visualization
params = list("n.samples"=n.samples, "n.exps"=n.exps, "n.iters"=n.iters,
              "pre.cond"=pre.cond, "alpha"=Alpha, "lazy.step"=lazy.step)
chain.list = list(out.ula0, out.ula, out.mala, out.mrw)
acc.list = list(mean(as.numeric(matrix(acc.ula0))),
                mean(as.numeric(matrix(acc.ula))),
                mean(as.numeric(matrix(acc.mala))),
                mean(as.numeric(matrix(acc.mrw))))

source("plots.R")
plt.err(chain.list, acc.list, params, maxiter = 2000, logScale=FALSE)

#______________________________________________________________________
# Benchmarking
library(rbenchmark)
benchmark(ula(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, delta=0.1,
               loglike=loglike, loglike.grad=loglike.grad, init=c(0, 0)), 
           replications=10)

benchmark(mala(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step,
                loglike=loglike, loglike.grad=loglike.grad, init=c(0, 0)), 
           replications=10)
benchmark(mrw(N=n.iters, d=dim.theta, m=m, L=L, lazy=lazy.step, #sigma=2.5,
                 loglike=loglike, init=c(0, 0)), 
           replications=10)



