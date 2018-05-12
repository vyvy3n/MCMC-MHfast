plt.legend <- function(loc="topleft"){
  color = c("red", "darkorange", "blue", "darkseagreen3")
  legend(loc, col = color, lty = c(1,1,1,1), bty = "n",
         legend = c(paste(format(chain.list[[1]]$accept, digits=3), "ULA"), 
                    paste(format(chain.list[[2]]$accept, digits=3), "ULA large"),
                    paste(format(chain.list[[3]]$accept, digits=3), "MALA"), 
                    paste(format(chain.list[[4]]$accept, digits=3), "MRW")))
}

########### 
# Density #
###########__________________________________________________
plt.density <- function(chain.list, save=TRUE, num=0){
  color = c("red", "darkorange", "blue", "darkseagreen3")
  par(mfrow = c(1,2))
  if (save==TRUE){
    dev.copy(png, paste("density-", num, ".png", sep=""), width =900, height = 500)
  }
  plot( density(chain.list[[1]]$chain[1,]), col=color[1], 
        xlab="", ylab="Density", main="theta 1", xlim=c(-0.2, 1.5))
  lines(density(chain.list[[2]]$chain[1,]), col=color[2])
  lines(density(chain.list[[3]]$chain[1,]), col=color[3])
  lines(density(chain.list[[4]]$chain[1,]), col=color[4])
  plt.legend()

  plot( density(chain.list[[1]]$chain[2,]), col=color[1], 
        xlab="", ylab="Density", main="theta 2", xlim=c(-0.2, 1.5))
  lines(density(chain.list[[2]]$chain[2,]), col=color[2])
  lines(density(chain.list[[3]]$chain[2,]), col=color[3])
  lines(density(chain.list[[4]]$chain[2,]), col=color[4])
  plt.legend()
  
  if (save==TRUE){dev.off()}
}

########### 
# ACF CCF #
###########__________________________________________________
plt.acfccf <- function(chain.list, save=TRUE, num=0){
  color = c("red", "darkorange", "blue", "darkseagreen3")
  par(mfrow = c(1,2))
  if (save==TRUE){
    dev.copy(png, paste("acfccf-", num, ".png", sep=""), width =900, height = 500)
  }
  par(mfrow = c(4,3))
  acf(chain.list[[1]]$chain[1,], main = "ULA: First")
  acf(chain.list[[1]]$chain[2,], main = "ULA: Second")
  ccf(chain.list[[1]]$chain[1,], chain.list[[1]]$chain[2,], main = "ULA: CCF")
  acf(chain.list[[2]]$chain[1,], main = "ULA Large: First")
  acf(chain.list[[2]]$chain[2,], main = "ULA Large: Second")
  ccf(chain.list[[2]]$chain[1,], chain.list[[2]]$chain[2,], main = "ULA Large: CCF")
  acf(chain.list[[3]]$chain[1,], main = "MALA: First")
  acf(chain.list[[3]]$chain[2,], main = "MALA: Second")
  ccf(chain.list[[3]]$chain[1,], chain.list[[3]]$chain[2,], main = "MALA: CCF")
  acf(chain.list[[4]]$chain[1,], main = "MRW: First")
  acf(chain.list[[4]]$chain[2,], main = "MRW: Second")
  ccf(chain.list[[4]]$chain[1,], chain.list[[4]]$chain[2,], main = "MRW: CCF")
  if (save==TRUE){dev.off()}
}

########### 
#  Trace  #
###########__________________________________________________
plt.trace <- function(chain.list, save=TRUE, num=0){
  color = c("red", "darkorange", "blue", "darkseagreen3")
  par(mfrow = c(1,2))
  if (save==TRUE){
    dev.copy(png, paste("trace-", num, ".png", sep=""), width =900, height = 500)
  }
  #__________________________________________________________
  N = length(chain.list[[1]]$chain[1,])
  n = as.integer(0.9*N)
  par(mfrow = c(2, 1))
  plot( tail(1:N, n), tail(chain.list[[1]]$chain[1,], n), col=color[1], type='l',
        xlab="iters", ylab="theta 1", main = "trace plot", ylim=c(0.2, 1.2))
  lines(tail(1:N, n), tail(chain.list[[2]]$chain[1,], n), col=color[2])
  lines(tail(1:N, n), tail(chain.list[[3]]$chain[1,], n), col=color[3])
  lines(tail(1:N, n), tail(chain.list[[4]]$chain[1,], n), col=color[4])
  plt.legend()
  
  plot( tail(1:N, n), tail(chain.list[[1]]$chain[2,], n), col=color[1], type='l',
        xlab="iters", ylab="theta 2", main = "trace plot", ylim=c(0.2, 1.2))
  lines(tail(1:N, n), tail(chain.list[[2]]$chain[2,], n), col=color[2])
  lines(tail(1:N, n), tail(chain.list[[3]]$chain[2,], n), col=color[3])
  lines(tail(1:N, n), tail(chain.list[[4]]$chain[2,], n), col=color[4])
  plt.legend()

  #__________________________________________________________          
  if (save==TRUE){dev.off()}
}