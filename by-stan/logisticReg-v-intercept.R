# Ref: https://github.com/danilofreire/r-scripts/blob/master/stan-logistic-regression.R

# intercept : alpha
# parameters: theta_1, theta_2

library(coda)
library(ggmcmc)
#__________________________________________________
# Generate Simulation data
library(extraDistr)   #rsign()
library(sigmoid)      #sigmoid()

# Real theta
Theta = matrix(c(1, 1))

geneSamples <- function(n, p=2){
  X <- matrix(nrow=n, ncol=p)
  Y <- matrix(nrow=n, ncol=1)
  set.seed(5)
  for (i in 1:p){
    X[,i] = rsign(n)  #draw samples from i.i.d. Rademacher distribution 
    #library(extraDistr)
  }
  for (i in 1:n){
    Y[i] = rbinom(1, 1, sigmoid(X[i, ] %*% Theta))}
  return (list("X"=X, "Y"=Y))
}

samples = geneSamples(n=100)
data.X = matrix(samples$X, ncol = 2)
# data.X = matrix(scale(samples$X), ncol = 2)
data.Y = matrix(samples$Y)[, 1]  #dim(matrix(samples$Y))=(50,1)

n = dim(data.X)[1]

######################################
### Logistic Regression with RStan ###
######################################
# Main ref: https://github.com/danilofreire/r-scripts/blob/master/stan-logistic-regression.R

# Load necessary packages
library(rstan)  # Stan interface for R
library(Zelig)  # Simulations for frequentist inference, data

# intercept : alpha
# parameters: theta_1, theta_2
m1 <- '
data {                          
  int<lower=0> N;                // number of observations
  int<lower=0,upper=1> Y[N];     // response Y
  vector[N] x_1;                 // independent variable X1
  vector[N] x_2;                 // independent variable X2
}
parameters {
  real alpha;                    // intercept
  real theta_1;                  // coef for x1
  real theta_2;                  // coef for x2
}
model {
  Y ~ bernoulli_logit(alpha + theta_1 * x_1 + theta_2 * x_2); 
}'

# Create a list with the chosen variables
data.list <- list(N = n, Y = data.Y, x_1 = data.X[, 1], x_2 = data.X[, 2])
str(data.list)

# Estimate the model
fit <- stan(model_code = m1, data = data.list, iter = 1000, chains = 4)
print(fit, digits = 3)

#______________________________________________________
# Compare with frequentist estimation (with Zelig) 
z1 <- zelig(Y ~ x_1 + x_2, model = "logit", data = data.frame(data.list))
summary(z1)

s1 <- setx(z1, x_1 =1 )   
print(sim(z1, x = s1))    
#______________________________________________________
# Using coda and ggmcmc to plot graphs
# library(coda)    #these two may be loaded before "rstan", see the top of this script
# library(ggmcmc)  #these two may be loaded before "rstan", see the top of this script

# Plotting some graphs and assessing convergence with the coda package.
# We can convert a fit object to coda with the following function
# (http://jeromyanglim.tumblr.com/post/91434443911/how-to-convert-a-stan-fit-object-to-work-with-coda-and)
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

fit.mcmc <- stan2coda(fit)

# Remove lp__ columns
fit.mcmc <- fit.mcmc[,1:3]

# codamenu() has all tests and plots one may want. Just type:
codamenu()
# Use an mcmc object
2
fit.mcmc  
# Then follow the menu instructions. 

#______________________________________________________
# Covergence Analysis
intercept_trace_1 = as.matrix(fit.mcmc[,1][1])
theta_1_trace_1   = as.matrix(fit.mcmc[,2][1])
theta_2_trace_1   = as.matrix(fit.mcmc[,3][1])

chain = matrix(c(theta_1_trace_1, theta_2_trace_1), nrow=2)
N     = dim(chain)[2]
err   = numeric(100)
ii    = 0

for (i in 1:N){  
  err[i] = norm(data.matrix(chain[, i] - Theta), type = "1")
  # err[i] = norm(data.matrix(rowMeans(
  #   matrix(chain[, 1:i], nrow = dim(chain)[1])) - Theta), type = "1")
}
err = err / dim(chain)[1]

plot.new()
plot(err, type='l', xlab="iterations", ylab="L1 error")

#______________________________________________________
# Summary on acceptance rate (accept_stat__)
sampler_params = get_sampler_params(fit, inc_warmup = FALSE)
sampler_params[[1]][1:5,]
# accept_stat__ stepsize__ treedepth__ n_leapfrog__ divergent__ energy__
# [1,]     0.8563656  0.8457639           2            7           0 48.22885
# [2,]     0.9253854  0.8457639           2            3           0 48.62060
# [3,]     1.0000000  0.8457639           3            7           0 48.36535
# [4,]     1.0000000  0.8457639           2            3           0 47.04900
# [5,]     0.9980155  0.8457639           2            7           0 46.36486

#______________________________________________________
# A Stan fit object can also be transformed into a ggmcmc object with ggs(). 

# # Remove lp__ columns
# fit.mcmc <- fit.mcmc[,1:3]

# Change the parameters' labels 
P <- data.frame(Parameter = c("beta[1]", "beta[2]", "beta[3]"),
                Label = c("Alpha", "Theta_1", "Theta_2"))
fit.ggmcmc <- ggs(fit.mcmc, par_labels = P)

ggs_traceplot(fit.ggmcmc) + ggtitle("Trace Plots") + theme_bw()
ggs_density(fit.ggmcmc) + ggtitle("Logistic Estimations for Voter Turnout") +
  xlab("Estimate") + ylab("Density") + theme_bw()
ggs_caterpillar(fit.ggmcmc) + ggtitle("Coefficient Plot") +
  xlab("HPD") + ylab("Parameter") + theme_bw()

# You may also check Stan's Github repository, it has many examples:
# https://github.com/stan-dev/example-models
# Several other ggmcmc() options here: http://xavier-fim.net/packages/ggmcmc/

# Lastly, if we want to estimate a full marginal distribution for a given
# predictor in the model, we employ the workhorse apply() function as follows:
# (from: http://bit.ly/1y0OMyC)
margins <- apply(as.matrix(fit), MARGIN = 2, FUN = quantile, probs = (1:100) / 100)
head(margins, 10)

# Plot the marginal distribution of educate (2nd column)
par(mfrow=c(1,2))
plot(jitter(margins[,2]), pch=20, xlab = "Theta 1 - Marginal Distribution (%)",
     ylab = "Probability of Y", main = "Predicted Values", axes=FALSE)
axis(1) # adds x axis
axis(2) # adds y axis
plot(jitter(margins[,3]), pch=20, xlab = "Theta 2 - Marginal Distribution (%)",
     ylab = "Probability of Y", main = "Predicted Values", axes=FALSE)
axis(1)
axis(2)
