##############
##   RWM    ## 
##############
#MRW sampler, with symmetric proposal.
mrw <- function(N, d, m, L, loglike, init=c(0, 0), sigma=FALSE, lazy=FALSE)
{
  # Calculate step size
  if (sigma==FALSE){ 
    h = 1/(d^2*L*L/m)
    sigma = sqrt(2*h); cat("sigma:", sigma, "\n")
  } else {
    h = sigma^2/2
  } 
  # Calculates the log of transition probability from y to x
  proplike <- function(x, mu)
  {
    return(as.numeric(- t((x - mu)) %*% (x - mu)/(2*sigma^2)))
  }
  #__________________________________________________
  mrw.chain <- matrix(0, nrow = d, ncol = N)
  accept <- 0
  mrw.chain[, 1] <- init #c(0, 0) #starting point
  
  for(i in 2:N)
  {
    coin = runif(1) * lazy
    if(coin>0.5){
      mrw.chain[, i] <- mrw.chain[, i-1]} #lazy step
    else{
      #__________________________________________________
      prop <- rnorm(2, mean = mrw.chain[, i-1], sd = sigma)
      log.ratio <- loglike(prop) - loglike(mrw.chain[, i-1]) +
                   proplike(mrw.chain[, i-1], prop) - proplike(prop, mrw.chain[, i-1])
      if(log(runif(1)) < log.ratio)
      {
        mrw.chain[, i] <- prop
        accept <- accept+1
      }
      else{
        mrw.chain[, i] <- mrw.chain[, i-1]
      }
      #__________________________________________________
    }
  }
  return(list("chain" = mrw.chain, "accept" = accept/N))
}


##############
##   MALA   ## 
##############
# Mala sampler
mala <- function(N, d, m, L, loglike, loglike.grad, init=c(0, 0), sigma=FALSE, lazy=FALSE)
{
  # Calculate step size
  if (sigma==FALSE){ 
    h = 1/L*min(1/sqrt(d*(L/m)), 1/d) #stepsize
    sigma = sqrt(2*h); cat("sigma:", sigma, "\n")
  } else {
    h = sigma^2/2
  } 
  # Calculates the log of transition probability from y to x
  proplike <- function(x, y, mu.u=FALSE)
  {
    if (mu.m==FALSE){
      grad <- loglike.grad(y) 
      mu.m <- y + sigma^2 * grad/2
    } 
    return(as.numeric(- t((x - mu.m)) %*% (x - mu.m)/(2*sigma^2)))
  }
  #__________________________________________________
  chain.mala <- matrix(0, nrow = d, ncol = N)
  accept <- 0
  chain.mala[, 1] <- c(0, 0)
  for(i in 2:N)
  {      
    coin = runif(1) * lazy
    if(coin>0.5){
      chain.mala[, i] <- chain.mala[, i-1]} #lazy step
    else{
      #__________________________________________________
      grad <- loglike.grad(chain.mala[, i-1]) #note: loglike = -f
      mu.m <- chain.mala[, i-1] + h * grad
      prop <- rnorm(d, mean = mu.m, sd = sigma)
      log.ratio <- loglike(prop) - loglike(chain.mala[, i-1]) +
                   proplike(chain.mala[, i-1], prop, mu.m) - proplike(prop, chain.mala[, i-1])
      #log.ratio <- min(0, log.ratio)
      if(log(runif(1)) < log.ratio)
      {
        chain.mala[, i] <- prop
        accept <- accept+1
      }
      else{
        chain.mala[, i] <- chain.mala[, i-1]
      }
      #__________________________________________________
    }
  }
  return(list("chain" = chain.mala, "accept" = accept/N))
}

##############
##   ULA    ## 
##############

ula <- function(N, d, m, L, loglike, loglike.grad, delta=1, #manually choose delta
                init=c(0, 0), sigma=FALSE, lazy=FALSE, verbose=FALSE) 
{
  # Calculate step size
  if (sigma==FALSE){ 
    h = delta^2/(d*L*L/m)  #stepsize
    sigma = sqrt(2*h); cat("sigma:", sigma, "\n")
  } else{
    h = sigma^2/2} 
  #__________________________________________________
  ula.chain <- matrix(0, nrow = d, ncol = N)
  accept <- 0
  ula.chain[, 1] <- c(0, 0) #starting point
  for(i in 2:N)
  {
    coin = runif(1) * lazy
    if(coin>0.5){
      ula.chain[, i] <- ula.chain[, i-1]} #lazy step
    else{
      #__________________________________________________
      grad <- loglike.grad(ula.chain[, i-1])
      mu.m  <- ula.chain[, i-1] + h * grad
      prop <- rnorm(d, mean = mu.m, sd = sigma)
      ula.chain[, i] <- prop
      
      # #Reply: ULA stands for unadjusted and hence there is no accept-reject step.
      # #While it is motivated well as a discretization of Langevin diffusion, the 
      # #point of our paper is to point out precisely the fact that one should use 
      # #accept-reject (at least for log-concave cases) to speed up the algorithm.
      #
      # log.ratio <- loglike(prop) - loglike(ula.chain[, i-1])
      # if(log(runif(1)) < log.ratio)
      # {
      #   ula.chain[, i] <- prop
      #   accept <- accept+1
      # }
      # else{
      #   ula.chain[, i] <- ula.chain[, i-1]
      # }
      if (verbose==TRUE){
        cat("Ratio:", exp(log.ratio), "\n")  #display only
        cat("[New move] ", ula.chain[, i], "\n")
      }
      #__________________________________________________
    }
  }
  return(list("chain" = ula.chain, "accept" = accept/N))
}
