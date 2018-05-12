source("helpingFun.R")  #source custom functions
plot.num = get.plot.number()

#__________________________________________________
# General params settings
# set.seed(1)
dim.theta = 2
n.exps    = 3
n.samples = 50
n.iters   = 1e3
pre.cond  = TRUE
lazy.step = FALSE

Alpha     = 0.9  #user-specified (for model) 

#__________________________________________________
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

samples = geneSamples(n=n.samples)
X = matrix(samples$X, ncol = 2)
# X = matrix(scale(samples$X), ncol = dim.theta)
Y = matrix(samples$Y)

#__________________________________________________
# Test if data was correctly generated
model <- glm(Y ~ ., family = binomial(link = 'logit'),
             data = data.frame(X))
# summary(model): coefficients should be [1, 1]
print(model$coefficients)

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

# L=20
# m=0.1

##exp42
#L=L/10

##exp43
# L=L/3
# m=m*3

# ##exp46
# geo_m = sqrt(L*m)
# L = sqrt(L*geo_m)
# m = sqrt(m*geo_m)
  
# ##exp47
# L = sqrt(L*m)

k       = L/m  #condition number
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

##############
##   MALA   ## 
##############

# MALA Params
H.mala = 1/L * min(1/sqrt(D*(L/m)), 1/D)  #stepsize (ref to Table 1 on page 8
# ## ULA step size
# delta  = 2
# H.mala = delta^2 / (D*k*L)

# Transition Probability between states x and z (newly proposed)
tprob <- function(z, x, h = H.mala){
  return(exp(- f(z) - norm(x - z + h*fgrad(z), type="2")^2/4/h))
}

# Mala sampler
mala <- function(N = 1e5, d = D, h = H.mala, lazy = TRUE, verbose = FALSE)  #h=stepsize
{
  mala.chain <- matrix(0, nrow = d, ncol = N)
  accept <- 0
  mala.chain[, 1] <- matrix(0, nrow = d, ncol = 1)  #starting point
  
  for(i in 2:N)
  {
    coin = runif(1) * lazy
    if(coin>0.5){
      mala.chain[, i] <- mala.chain[, i-1]} #lazy step
    else{
      grad <- fgrad(mala.chain[, i-1])
      mu   <- mala.chain[, i-1] - h * fgrad(mala.chain[, i-1])
      prop <- rnorm(d, mean = mu, sd = sqrt(2*h))
      
      x = mala.chain[, i-1]
      z = prop
      #print(c(tprob(z, x), tprob(x, z)))
      ratio <- min(1, tprob(z, x) / tprob(x, z))
      if(verbose==TRUE){cat("Ratio of accepting the proposed point:", ratio, "\n")}  #display only
      if(runif(1) >= ratio){
        # Reject
        mala.chain[, i] <- mala.chain[, i-1]
      }
      else{
        # Accept
        mala.chain[, i] <- prop
        accept <- accept+1
        if(verbose==TRUE){cat("[New move] ", mala.chain[, i], "\n")}  #for display only
      }
    }
  }
  return(list("chain" = mala.chain, "accept" = accept/N))
}

##############
##   MRW    ## 
##############

#MRW Params
H.mrw = 1 / (D^2*k*L)  #stepsize

#MRW sampler, with symmetric proposal.
mrw <- function(N = 1e5, d = D, h = H.mrw, lazy = TRUE)
{
  mrw.chain <- matrix(0, nrow = d, ncol = N)
  accept <- 0
  mrw.chain[, 1] <- matrix(0, nrow = d, ncol = 1)  #starting point
  
  for(i in 2:N)
  {
    coin = runif(1) * lazy
    if(coin>0.5){
      mrw.chain[, i] <- mrw.chain[, i-1]} #lazy step
    else{
      prop <- rnorm(2, mean = mrw.chain[, i-1], sd = sqrt(2*h))
      log.ratio <- loglike(prop) - loglike(mrw.chain[, i-1])
      if(log(runif(1)) < log.ratio)
      {
        mrw.chain[, i] <- prop
        accept <- accept+1
      }
      else{
        mrw.chain[, i] <- mrw.chain[, i-1]
      }
    }
  }
  return(list("chain" = mrw.chain, "accept" = accept/N))
}

##############
##   ULA    ## 
##############

ula <- function(N = 1e5, delta=0.1, d = D, h = H.ula, lazy = TRUE)  #manually choose delta
{
  # ULA Params
  H.ula = delta^2 / (D*k*L)  #stepsize
  
  ula.chain <- matrix(0, nrow = d, ncol = N)
  accept <- 0
  ula.chain[, 1] <- matrix(0, nrow = d, ncol = 1)  #starting point
  for(i in 2:N)
  {
    coin = runif(1) * lazy
    if(coin>0.5){
      ula.chain[, i] <- ula.chain[, i-1]} #lazy step
    else{
      grad <- fgrad(ula.chain[, i-1])
      mu   <- ula.chain[, i-1] - h * fgrad(ula.chain[, i-1])
      prop <- rnorm(d, mean = mu, sd = sqrt(2*h))
      
      log.ratio <- loglike(prop) - loglike(ula.chain[, i-1])
      if(log(runif(1)) < log.ratio)
      {
        ula.chain[, i] <- prop
        accept <- accept+1
      }
      else{
        ula.chain[, i] <- ula.chain[, i-1]
      }
    }
  }
  return(list("chain" = ula.chain, "accept" = accept/N))
}

###################################################
#__________________________________________________
# choose a method 
# chain.mrw  <- mrw(N = 1e4*5)  #each column is a point sampled
# chain.mala <- mala(N = 1e3)   #each column is a point sampled
# chain.ula  <- ula(N = 1e4)    #each column is a point sampled

#__________________________________________________
# A function for getting the chian
library(base) #colMeans()

get.chain <- function(method.chose, N=4*1e3, del=FALSE, lazy=lazy.step){
  if (del == FALSE){
    chain   = method.chose(N, lazy=lazy)}
  else{
    chain   = method.chose(N, del, lazy=lazy)}
  cat("\nAccept rate: ", chain$accept)
  return(list('chain'=covxp %*% chain$chain, #for pre-condition reverse
              'rate'=chain$accept)) 
}

#__________________________________________________
# A function for calculating L1 err of mean estimation over iterations
cal.l1.err <- function(chain){
  N.chain = dim(chain)[2]
  err     = numeric(N.chain)
  
  for (i in 1:N.chain){  
    err[i] = norm(data.matrix(chain[, i] - Theta), type = "1")
    # err[i] = norm(data.matrix(rowMeans(
    #   matrix(chain[, 1:i], nrow = dim(chain)[1])) - Theta), type = "1")
  }
  return(err / dim(chain)[1])
}

#__________________________________________________
# Summarize Repeated Experiment

repeat.exps <- function(method.chose, N=n.iters, del=FALSE, n.exps=5, d.theta=dim.theta){
  #initalization
  err.matrix = matrix(nrow = n.exps, ncol = N)
  chain.list = list()
  avg_rate   = 0
  #calculate L1 error & summarize theta
  for (n in seq(1, n.exps)){
    out = get.chain(method.chose=method.chose, N=N, del=del)
    err.matrix[n, ] = cal.l1.err(chain=out$chain)
    chain.list[[n]] = out$chain
    avg_rate = avg_rate+out$rate
  }
  return(list(err.matrix, chain.list, avg_rate/n.exps))
}

res1 = repeat.exps(method.chose=ula, del=0.1, n.exps=n.exps) #ULA
res2 = repeat.exps(method.chose=ula, del=1, n.exps=n.exps)   #ULA large
res3 = repeat.exps(method.chose=mala, n.exps=n.exps)
res4 = repeat.exps(method.chose=mrw, n.exps=n.exps)

# res1 = out1[1]
# res2 = out2[1]
# res3 = out3[1]
# res4 = out4[1]

#__________________________________________________
# Plotting L1 err v.s. iters () 

err1 = res1[[1]]
err2 = res2[[1]]
err3 = res3[[1]]
err4 = res4[[1]]

dev.copy(png, paste("Rplot", plot.num, "-", "err.png", sep=""), width = 600, height = 500)
par(mfrow=c(1, 1))
color = c("red", "darkorange", "blue", "darkseagreen3")
plot( log(colMeans(err1))/log(10), type="l", lwd=2, col=color[1], 
      xlim=c(0, n.iters), ylim=c(-1.4, 0), ann=FALSE)
lines(log(colMeans(err2))/log(10), type="l", lw=2, col=color[2])
lines(log(colMeans(err3))/log(10), type="l", lw=2, col=color[3])
lines(log(colMeans(err4))/log(10), type="l", lw=2, col=color[4])
legend(0.75 * n.iters, -.75, col = color, lty = c(1,1,1,1), bty = "n",
       legend = c(paste(format(res1[[3]], digits=4), "ULA"), 
                  paste(format(res2[[3]], digits=4), "ULA large"),
                  paste(format(res3[[3]], digits=4), "MALA"), 
                  paste(format(res4[[3]], digits=4), "MRW")))
#Use option bty = "n" in legend to remove the box around the legend
title(paste("n.samples=", n.samples, ", n.exps=", n.exps, 
            ", pre-cond=", pre.cond, ", alpha=", Alpha, sep=''),
      xlab="iterations", ylab="L1 error log_10 scale")
dev.off()

#__________________________________________________
# Traceplot of theta

plot.theta.trace <- function(res, method){
  plot( res[[2]][[1]][1, ], type="l", col="blue",
        main=method, xlab="iterations", ylab="value")
  lines(res[[2]][[1]][2, ], type="l", col="deeppink")
  legend("bottomright", col=c("blue", "deeppink"), lty=c(1,1), bty = "n",
         legend=c("theta[1]", "theta[2]"))
}

dev.copy(png, paste("Rplot", plot.num, "-", "traceplot.png", sep=""), width = 870, height = 524)
par(mfrow=c(2, 2))
plot.theta.trace(res1, "ULA")
plot.theta.trace(res2, "ULA large")
plot.theta.trace(res3, "MALA")
plot.theta.trace(res4, "MRW")
title(paste("Traceplots,", "n.samples =", n.samples, ", n.exps =", n.exps, 
            ", pre-cond =", pre.cond), outer=TRUE, line = -1)
dev.off()

#__________________________________________________
# Plotting theta summary

theta.summary <- function(theta.list, method, k.chose=1){
  l.theta = dim(theta.list[[1]])[1]
  l.iters = dim(theta.list[[1]])[2]
  l.exps  = length(theta.list)
  
  temp = list()
  for (k in seq(1, l.theta)){
    temp[[k]] <- matrix(0, nrow=l.exps, ncol=l.iters)
  }
  ## Minus real Theta from each sampled theta & statistics
  for (i in seq(1, l.exps)){
    for (j in seq(1, l.iters)){
      theta.list[[i]][, j] = theta.list[[i]][, j] - Theta
      for (k in seq(1, l.theta)){
        temp[[k]][i, j] = theta.list[[i]][k, j]
      }
    }
  }
  ## Statistic Summary
  
  # par(mfrow=c(l.theta, 1))
  # for (k in seq(1, l.theta)){
  #   plot(colMeans(temp[[k]]), type="l")
  # }
  par(mfrow=c(1, 1))
  k = k.chose
  
  dev.copy(png, paste("Rplot", plot.num, "-", tolower(method), "-", k.chose, ".png", sep=""), 
           width = 600, height = 500)
  theta.mean = colMeans(temp[[k]])
  plot(theta.mean, type="l", main=method, xlab="iterations", ylab="error")
  #add error bars
  x = c(c(10, 100, 300, 500), 
        seq(1000, l.iters, length.out=6))
  y = matrix(0, nrow=1, ncol=length(x))
  q = matrix(0, nrow=2, ncol=length(x))
  i = 0
  for (ix in x){
    i      = i + 1
    y[i]   = theta.mean[ix]
    q[, i] = matrix(quantile(temp[[k]][, ix], c(.25, .75), name=FALSE))
    #q[1, ] is y.upper: y 75%
    #q[2, ] is y.lower: y 25%
  }
  #style of errorbars:
  myangle  = 90   #Angle from shaft to edge of arrow head.
  mylength = 0.1 #Length of arrow head in inches (!!!)
  #two--sided errorbars:
  par(new=TRUE)
  # plot(x, y, ylim=range(q[1, ], q[2, ]), ann=FALSE, axes=FALSE)
  arrows(x, q[1, ], x, q[2, ], code = 3, angle=myangle, length=mylength)
  dev.off()
}


theta.summary(res3[[2]], method="MALA", k.chose=1)
theta.summary(res3[[2]], method="MALA", k.chose=2)
 
theta.summary(res1[[2]], method="ULA", k.chose=1)
theta.summary(res2[[2]], method="ULA-large", k.chose=1)
theta.summary(res4[[2]], method="MRW", k.chose=1)